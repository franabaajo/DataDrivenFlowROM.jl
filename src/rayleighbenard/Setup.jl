import IncompressibleNavierStokes as NS

include(joinpath(@__DIR__, "..", "workflows", "Backend.jl"))

getbackend() = workflow_getbackend(; env_var = "RAYLEIGH_BACKEND")
ensure_cuda_psolver!(backend) = ensure_workflow_cuda_psolver!(backend)
backend_summary() = workflow_backend_summary(; env_var = "RAYLEIGH_BACKEND")

function grid_axis(lims::NTuple{2,Float64}, n::Int, stretch::Union{Nothing,Float64})
    if isnothing(stretch)
        return range(lims..., n + 1)
    end
    return NS.tanh_grid(lims..., n, stretch)
end

function getsetup(cfg::FOMConfig; backend = getbackend())
    ax = grid_axis(cfg.xlims, cfg.nx, cfg.stretch_x)
    ay = grid_axis(cfg.ylims, cfg.ny, cfg.stretch_y)

    return NS.Setup(;
        x = (ax, ay),
        boundary_conditions = (;
            u = (
                (NS.DirichletBC(), NS.DirichletBC()),
                (NS.DirichletBC(), NS.DirichletBC()),
            ),
            temp = (
                (NS.SymmetricBC(), NS.SymmetricBC()),
                (NS.DirichletBC(cfg.temp_hot), NS.DirichletBC(cfg.temp_cold)),
            ),
        ),
        backend = backend,
    )
end

boussinesq_params(cfg::FOMConfig) = (
    viscosity = cfg.viscosity,
    conductivity = cfg.conductivity,
    gravity = cfg.gravity,
    gdir = cfg.gdir,
    dodissipation = cfg.dodissipation,
)

const _RAYLEIGH_PROFILE_CONDUCTION = UInt8(0)
const _RAYLEIGH_PROFILE_UNIFORM = UInt8(1)

function _rayleigh_profile_kind(profile::AbstractString)
    prof = lowercase(strip(profile))
    if prof in ("conduction", "linear", "plates")
        return _RAYLEIGH_PROFILE_CONDUCTION
    elseif prof in ("uniform", "mean", "flat")
        return _RAYLEIGH_PROFILE_UNIFORM
    else
        error("Unknown initial_temperature_profile: $(profile). Use \"conduction\" or \"uniform\".")
    end
end

struct InitialTemperatureFunctor{T}
    temp_hot::T
    temp_cold::T
    y0::T
    yrange::T
    perturbation_amplitude::T
    perturbation_wavenumber::T
    profile_kind::UInt8
end

@inline function (f::InitialTemperatureFunctor)(x, y)
    yfrac = (y - f.y0) / f.yrange
    base = if f.profile_kind == _RAYLEIGH_PROFILE_CONDUCTION
        f.temp_hot + (f.temp_cold - f.temp_hot) * yfrac
    else
        one(y) * ((f.temp_hot + f.temp_cold) / 2)
    end
    # Match the RayleighBenard2D.jl example:
    # temp(x,y) = 0.5 + max(sinpi(20*x)/100, 0)
    # Here `perturbation_wavenumber` is in physical x-coordinates (not normalized).
    perturb = f.perturbation_amplitude * max(sinpi(f.perturbation_wavenumber * x), zero(x))
    return base + perturb
end

function initial_temperature_functor(cfg::FOMConfig, ::Type{T}) where {T}
    y0 = T(cfg.ylims[1])
    yrange = T(cfg.ylims[2] - cfg.ylims[1])
    kind = _rayleigh_profile_kind(cfg.initial_temperature_profile)
    return InitialTemperatureFunctor(
        T(cfg.temp_hot),
        T(cfg.temp_cold),
        y0,
        yrange,
        T(cfg.perturbation_amplitude),
        T(cfg.perturbation_wavenumber),
        kind,
    )
end

function initial_temperature(cfg::FOMConfig, x, y)
    return initial_temperature_functor(cfg, typeof(x))(x, y)
end

function initial_state(cfg::FOMConfig, setup, psolver)
    u = NS.velocityfield(setup, (dim, x, y) -> zero(x); psolver = psolver)
    tempfunc = initial_temperature_functor(cfg, eltype(setup.x[1]))
    temp = NS.temperaturefield(setup, tempfunc)
    return (; u, temp)
end

state_vector_length(state) = length(state.u) + length(state.temp)

function copy_state_to_vector!(dest::AbstractVector, state)
    nu = length(state.u)
    nt = length(state.temp)
    length(dest) == nu + nt || error("State vector length mismatch")

    copyto!(view(dest, 1:nu), vec(Array(state.u)))
    copyto!(view(dest, (nu + 1):(nu + nt)), vec(Array(state.temp)))
    return dest
end

function copy_vector_to_state!(state, x::AbstractVector)
    nu = length(state.u)
    nt = length(state.temp)
    length(x) == nu + nt || error("State vector length mismatch")

    # NOTE: Avoid passing `view`/`reshape(view(...))` into `copyto!` when the
    # destination lives on the GPU. CUDA.jl does not support bulk host->device
    # copies from `SubArray` sources, and will fall back to scalar indexing
    # (which is disallowed and extremely slow).
    #
    # Using the linear-index `copyto!(dest, doffs, src, soffs, n)` path keeps
    # the source as a plain vector and triggers the efficient bulk copy.
    copyto!(vec(state.u), 1, x, 1, nu)
    copyto!(vec(state.temp), 1, x, nu + 1, nt)
    return state
end

function snapshot_memory_gb(len_state::Integer, nsnapshot::Integer; T::Type = Float64)
    numbersize = sizeof(T)
    size_bytes = len_state * 2 * (nsnapshot + 1) * numbersize
    return size_bytes / 1e9
end

"""
Return a nonnegative weight vector `w` such that `x' * Diagonal(w) * x / 2`
approximates the discrete kinetic energy of a staggered velocity field `x = vec(u)`.

Weights are the control-volume sizes associated with each velocity DOF and are
zero outside the physical velocity index sets `setup.Iu`.
"""
function velocity_mass_vector(setup, utemplate)
    D = setup.dimension()
    D in (2, 3) || error("velocity_mass_vector only supports D=2 or D=3, got D=$D")

    T = eltype(utemplate)
    wfield = similar(utemplate, T, size(utemplate))
    wfield .= zero(T)

    w_lastdim = ndims(wfield)
    w = ntuple(D) do α
        selectdim(wfield, w_lastdim, α)
    end

    if D == 2
        Δ1x = Array(setup.Δu[1])
        Δ1y = Array(setup.Δ[2]')
        Δ2x = Array(setup.Δ[1])
        Δ2y = Array(setup.Δu[2]')

        area1 = @. Δ1x * Δ1y
        area2 = @. Δ2x * Δ2y

        view(w[1], setup.Iu[1]) .= area1[setup.Iu[1]]
        view(w[2], setup.Iu[2]) .= area2[setup.Iu[2]]
    else
        Δx = Array(setup.Δ[1])
        Δy = Array(setup.Δ[2])
        Δz = Array(setup.Δ[3])
        Δux = Array(setup.Δu[1])
        Δuy = Array(setup.Δu[2])
        Δuz = Array(setup.Δu[3])

        vol1 = @. reshape(Δux, :, 1, 1) * reshape(Δy, 1, :, 1) * reshape(Δz, 1, 1, :)
        vol2 = @. reshape(Δx, :, 1, 1) * reshape(Δuy, 1, :, 1) * reshape(Δz, 1, 1, :)
        vol3 = @. reshape(Δx, :, 1, 1) * reshape(Δy, 1, :, 1) * reshape(Δuz, 1, 1, :)

        view(w[1], setup.Iu[1]) .= vol1[setup.Iu[1]]
        view(w[2], setup.Iu[2]) .= vol2[setup.Iu[2]]
        view(w[3], setup.Iu[3]) .= vol3[setup.Iu[3]]
    end

    return vec(wfield)
end

"""
Return a nonnegative weight vector `w` such that `x' * Diagonal(w) * x / 2`
approximates the discrete L2 energy of a scalar field `x = vec(p)` at cell centers.

Weights are the control-volume sizes associated with each scalar DOF and are
zero outside the physical pressure index set `setup.Ip`.
"""
function scalar_mass_vector(setup, ptemplate)
    D = setup.dimension()
    D in (2, 3) || error("scalar_mass_vector only supports D=2 or D=3, got D=$D")

    T = eltype(ptemplate)
    wfield = similar(ptemplate, T, size(ptemplate))
    wfield .= zero(T)

    if D == 2
        Δx = Array(setup.Δ[1])
        Δy = Array(setup.Δ[2]')
        area = @. Δx * Δy
        view(wfield, setup.Ip) .= area[setup.Ip]
    else
        Δx = Array(setup.Δ[1])
        Δy = Array(setup.Δ[2])
        Δz = Array(setup.Δ[3])
        vol = @. reshape(Δx, :, 1, 1) * reshape(Δy, 1, :, 1) * reshape(Δz, 1, 1, :)
        view(wfield, setup.Ip) .= vol[setup.Ip]
    end

    return vec(wfield)
end

"""
Build row weights for a concatenated state vector `x = [vec(u); vec(temp)]`.
"""
function state_mass_vector(setup, utemplate, temp_template)
    wu = velocity_mass_vector(setup, utemplate)
    wt = scalar_mass_vector(setup, temp_template)
    return vcat(wu, wt)
end
