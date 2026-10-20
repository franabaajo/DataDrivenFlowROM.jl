if !isdefined(@__MODULE__, :_WORKFLOW_FOM_LOADED)
const _WORKFLOW_FOM_LOADED = true

using JLD2
using Mmap
using ProgressMeter
import IncompressibleNavierStokes as NS

const WORKFLOW_JLD2_IOTYPE = IOStream

function workflow_load_object_iostream(path::AbstractString; iotype::Type = WORKFLOW_JLD2_IOTYPE)
    return jldopen(path, "r"; iotype = iotype) do file
        all_keys = collect(keys(file))
        isempty(all_keys) && error("File $path does not contain any object")
        length(all_keys) > 1 && error("File $path contains more than one object; use jldopen directly.")
        return file[all_keys[1]]
    end
end

function _workflow_checked_mmap_matrix(path::AbstractString, ::Type{T}, dims::Tuple{Int,Int}; label::AbstractString) where {T}
    expected = Int(prod(dims)) * sizeof(T)
    actual = filesize(path)
    if actual != expected
        error(
            "Snapshot stream file size mismatch for $(label):\n" *
            "  file = $(path)\n" *
            "  expected_bytes = $(expected)\n" *
            "  actual_bytes = $(actual)\n" *
            "This typically means the streamed snapshots are incomplete/corrupted. " *
            "Re-run the FOM (or delete the run directory and re-run).",
        )
    end
    return Mmap.mmap(path, Matrix{T}, dims)
end

function snapshot_stream_paths(dir::AbstractString)
    return (
        u = joinpath(dir, "snapshots_u.bin"),
        dudt = joinpath(dir, "snapshots_dudt.bin"),
        meta = joinpath(dir, "snapshots_meta.jld2"),
    )
end

function has_streamed_snapshots(dir::AbstractString)
    paths = snapshot_stream_paths(dir)
    return isfile(paths.u) && isfile(paths.dudt) && isfile(paths.meta)
end

function fullrun(cfg, u, setup, psolver)
    state = (; u = copy(u))
    method = NS.LMWray3()
    ode_cache = NS.get_cache(method, state, setup)
    force_cache = NS.get_cache(NS.navierstokes!, setup)
    stepper = NS.create_stepper(method; setup, psolver, state, t = 0.0)

    dudt = similar(u)
    snaps = (; u = zeros(Float64, length(u), cfg.nsnapshot + 1),
        dudt = zeros(Float64, length(u), cfg.nsnapshot + 1))
    cpubuf = Array(u)

    prog = Progress(cfg.nsnapshot + 1; desc = "Running full-order model")

    for isnap in 0:cfg.nsnapshot
        tstop = isnap / cfg.nsnapshot * cfg.tsim
        isubstep = 0

        while stepper.t < prevfloat(tstop)
            dt = cfg.cfl * NS.propose_timestep(NS.navierstokes!, stepper.state, setup, (; viscosity = cfg.viscosity))
            dt = min(dt, tstop - stepper.t)

            stepper = NS.timestep!(
                method,
                NS.navierstokes!,
                stepper,
                dt;
                params = (; viscosity = cfg.viscosity),
                ode_cache,
                force_cache,
            )
            isubstep += 1
        end

        NS.navierstokes!(
            (; u = dudt),
            state,
            stepper.t;
            setup,
            cache = force_cache,
            viscosity = cfg.viscosity,
        )
        NS.project!(dudt, setup; psolver, ode_cache.p)

        copyto!(cpubuf, state.u)
        copyto!(view(snaps.u, :, isnap + 1), cpubuf)

        copyto!(cpubuf, dudt)
        copyto!(view(snaps.dudt, :, isnap + 1), cpubuf)

        next!(prog; showvalues = [("snapshot", isnap), ("substeps", isubstep), ("time", stepper.t)])
    end

    return snaps
end

function fullrun_streaming(cfg, u, setup, psolver, paths)
    state = (; u = copy(u))
    method = NS.LMWray3()
    ode_cache = NS.get_cache(method, state, setup)
    force_cache = NS.get_cache(NS.navierstokes!, setup)
    stepper = NS.create_stepper(method; setup, psolver, state, t = 0.0)

    dudt = similar(u)
    cpubuf = Array(u)

    prog = Progress(cfg.nsnapshot + 1; desc = "Running full-order model")

    open(paths.u, "w") do io_u
        open(paths.dudt, "w") do io_d
            for isnap in 0:cfg.nsnapshot
                tstop = isnap / cfg.nsnapshot * cfg.tsim
                isubstep = 0

                while stepper.t < prevfloat(tstop)
                    dt = cfg.cfl * NS.propose_timestep(NS.navierstokes!, stepper.state, setup, (; viscosity = cfg.viscosity))
                    dt = min(dt, tstop - stepper.t)

                    stepper = NS.timestep!(
                        method,
                        NS.navierstokes!,
                        stepper,
                        dt;
                        params = (; viscosity = cfg.viscosity),
                        ode_cache,
                        force_cache,
                    )
                    isubstep += 1
                end

                NS.navierstokes!(
                    (; u = dudt),
                    state,
                    stepper.t;
                    setup,
                    cache = force_cache,
                    viscosity = cfg.viscosity,
                )
                NS.project!(dudt, setup; psolver, ode_cache.p)

                copyto!(cpubuf, state.u)
                write(io_u, cpubuf)

                copyto!(cpubuf, dudt)
                write(io_d, cpubuf)

                next!(prog; showvalues = [("snapshot", isnap), ("substeps", isubstep), ("time", stepper.t)])
            end
        end
    end

    return nothing
end

function one_step_fom_from_state(
    u0_vec::AbstractVector,
    setup,
    psolver,
    cfg;
    dt::Real,
    method::Symbol = :explicit_euler,
)
    T = eltype(u0_vec)
    ufield = NS.vectorfield(setup)
    length(u0_vec) == length(ufield) ||
        error("Length mismatch: u0 has length $(length(u0_vec)) but setup expects $(length(ufield))")

    copyto!(ufield, reshape(u0_vec, size(ufield)))
    state = (; u = ufield)

    if method == :explicit_euler
        rk = NS.RKMethods.FE11(; T = T)
        ode_cache = NS.get_cache(rk, state, setup)
        force_cache = NS.get_cache(NS.navierstokes!, setup)
        stepper = NS.create_stepper(rk; setup, psolver, state, t = 0.0)
        stepper = NS.timestep!(
            rk,
            NS.navierstokes!,
            stepper,
            dt;
            params = (; viscosity = cfg.viscosity),
            ode_cache,
            force_cache,
        )
        u1 = vec(Array(stepper.state.u))
        return u1, (method = "explicit_euler", is_explicit_euler = true)
    elseif method == :lmwray3
        rk = NS.LMWray3(; T = T)
        ode_cache = NS.get_cache(rk, state, setup)
        force_cache = NS.get_cache(NS.navierstokes!, setup)
        stepper = NS.create_stepper(rk; setup, psolver, state, t = 0.0)
        stepper = NS.timestep!(
            rk,
            NS.navierstokes!,
            stepper,
            dt;
            params = (; viscosity = cfg.viscosity),
            ode_cache,
            force_cache,
        )
        u1 = vec(Array(stepper.state.u))
        return u1, (method = "lmwray3", is_explicit_euler = false)
    else
        error("Unsupported one-step FOM method: $(method). Use :explicit_euler or :lmwray3.")
    end
end

function load_streamed_snapshots(dir::AbstractString)
    stream_paths = snapshot_stream_paths(dir)
    has_streamed_snapshots(dir) || error("Streamed snapshots not found in: $dir")

    meta = workflow_load_object_iostream(stream_paths.meta)
    len = Int(meta.len)
    nsnapshot = Int(meta.nsnapshot)
    K = nsnapshot + 1

    T = meta.eltype == "Float64" ? Float64 : error("Unsupported snapshot eltype: $(meta.eltype)")
    dims = (len, K)
    u = _workflow_checked_mmap_matrix(stream_paths.u, T, dims; label = "snapshots_u")
    dudt = _workflow_checked_mmap_matrix(stream_paths.dudt, T, dims; label = "snapshots_dudt")
    return (; u, dudt)
end

end
