import IncompressibleNavierStokes as NS

include(joinpath(@__DIR__, "..", "workflows", "Backend.jl"))

getbackend() = workflow_getbackend(; env_var = "DIPOLE_BACKEND")
ensure_cuda_psolver!(backend) = ensure_workflow_cuda_psolver!(backend)
backend_summary() = workflow_backend_summary(; env_var = "DIPOLE_BACKEND")

function grid_axis(cfg::FOMConfig)
    if isnothing(cfg.stretch)
        return range(cfg.lims..., cfg.n + 1)
    else
        return NS.tanh_grid(cfg.lims..., cfg.n, cfg.stretch)
    end
end

function getsetup(cfg::FOMConfig; backend = getbackend())
    ax = grid_axis(cfg)
    return NS.Setup(;
        x = (ax, ax),
        boundary_conditions = (;
            u = (
                (NS.DirichletBC(), NS.DirichletBC()),
                (NS.DirichletBC(), NS.DirichletBC()),
            ),
        ),
        backend = backend,
    )
end

function dipole_velocity(cfg::FOMConfig, dim, x, y)
    x1, y1 = cfg.p1
    x2, y2 = cfg.p2
    r1 = (x - x1)^2 + (y - y1)^2
    r2 = (x - x2)^2 + (y - y2)^2
    d1 = ifelse(dim == 1, -(y - y1), x - x1)
    d2 = ifelse(dim == 1,  (y - y2), -(x - x2))
    return d1 * exp(-r1 / cfg.r0^2) + d2 * exp(-r2 / cfg.r0^2)
end

function initial_velocity(cfg::FOMConfig, setup, psolver)
    x1, y1 = cfg.p1
    x2, y2 = cfg.p2
    r0 = cfg.r0

    Ufun(dim, x, y) = begin
        r1 = (x - x1)^2 + (y - y1)^2
        r2 = (x - x2)^2 + (y - y2)^2
        d1 = ifelse(dim == 1, -(y - y1), x - x1)
        d2 = ifelse(dim == 1,  (y - y2), -(x - x2))
        d1 * exp(-r1 / r0^2) + d2 * exp(-r2 / r0^2)
    end

    return NS.velocityfield(setup, Ufun; psolver = psolver)
end

function scalevelocity!(u, setup, cfg::FOMConfig)
    u1 = selectdim(u, 3, 1)
    u2 = selectdim(u, 3, 2)

    Δ1x = setup.Δu[1]
    Δ1y = setup.Δ[2]'
    Δ2x = setup.Δ[1]
    Δ2y = setup.Δu[2]'

    eu = @. Δ1x * Δ1y * u1^2 / 2
    ev = @. Δ2x * Δ2y * u2^2 / 2

    kin = sum(view(eu, setup.Iu[1])) + sum(view(ev, setup.Iu[2]))
    @. u = sqrt(cfg.initialenergy / kin) * u
    return u
end

function snapshot_memory_gb(cfg::FOMConfig; T::Type = Float64)
    numbersize = sizeof(T)
    ncomponent = 2
    nfield = 2
    size_bytes = cfg.n^2 * ncomponent * nfield * (cfg.nsnapshot + 1) * numbersize
    return size_bytes / 1e9
end

function velocity_mass_vector(setup, utemplate)
    T = eltype(utemplate)
    wfield = similar(utemplate, T, size(utemplate))
    wfield .= zero(T)

    w1 = selectdim(wfield, 3, 1)
    w2 = selectdim(wfield, 3, 2)

    Δ1x = Array(setup.Δu[1])
    Δ1y = Array(setup.Δ[2]')
    Δ2x = Array(setup.Δ[1])
    Δ2y = Array(setup.Δu[2]')

    area1 = @. Δ1x * Δ1y
    area2 = @. Δ2x * Δ2y

    view(w1, setup.Iu[1]) .= area1[setup.Iu[1]]
    view(w2, setup.Iu[2]) .= area2[setup.Iu[2]]

    return vec(wfield)
end
