using Dates
using JLD2
using KernelAbstractions
using Mmap
using ProgressMeter
import IncompressibleNavierStokes as NS

function snapshot_stream_paths(dir::AbstractString)
    return (
        state = joinpath(dir, "snapshots_state.bin"),
        dstate = joinpath(dir, "snapshots_dstate.bin"),
        meta = joinpath(dir, "snapshots_meta.jld2"),
    )
end

function has_streamed_snapshots(dir::AbstractString)
    paths = snapshot_stream_paths(dir)
    return isfile(paths.state) && isfile(paths.dstate) && isfile(paths.meta)
end

function snapshot_time(cfg::FOMConfig, isnap::Integer)
    if cfg.nsnapshot == 0
        return 0.0
    end
    return isnap / cfg.nsnapshot * cfg.tsim
end

function fullrun_streaming(cfg::FOMConfig, start, setup, psolver, paths)
    state = (; u = copy(start.u), temp = copy(start.temp))
    method = NS.LMWray3()
    ode_cache = NS.get_cache(method, state, setup)
    force_cache = NS.get_cache(NS.boussinesq!, setup)
    stepper = NS.create_stepper(method; setup, psolver, state, t = 0.0)

    force = (; u = similar(start.u), temp = similar(start.temp))
    params = boussinesq_params(cfg)

    u_cpu = Array(start.u)
    temp_cpu = Array(start.temp)
    du_cpu = similar(u_cpu)
    dtemp_cpu = similar(temp_cpu)

    prog = Progress(cfg.nsnapshot + 1; desc = "Running full-order model")

    open(paths.state, "w") do io_state
        open(paths.dstate, "w") do io_dstate
            for isnap in 0:cfg.nsnapshot
                tstop = snapshot_time(cfg, isnap)
                isubstep = 0

                while stepper.t < prevfloat(tstop)
                    dt = cfg.cfl * NS.propose_timestep(NS.boussinesq!, stepper.state, setup, params)
                    dt = min(dt, tstop - stepper.t)

                    stepper = NS.timestep!(
                        method,
                        NS.boussinesq!,
                        stepper,
                        dt;
                        params,
                        ode_cache,
                        force_cache,
                    )
                    isubstep += 1
                end

                NS.boussinesq!(force, stepper.state, stepper.t; setup, cache = force_cache, params...)
                NS.project!(force.u, setup; psolver, ode_cache.p)

                copyto!(u_cpu, stepper.state.u)
                copyto!(temp_cpu, stepper.state.temp)
                write(io_state, u_cpu)
                write(io_state, temp_cpu)

                copyto!(du_cpu, force.u)
                copyto!(dtemp_cpu, force.temp)
                write(io_dstate, du_cpu)
                write(io_dstate, dtemp_cpu)

                next!(prog; showvalues = [("snapshot", isnap), ("substeps", isubstep), ("time", stepper.t)])
            end
        end
    end

    return (len_u = length(start.u), len_temp = length(start.temp))
end

function run_fom(
    cfg::FOMConfig;
    force::Bool = false,
    skip_if_exists::Bool = true,
    run_tag::Union{Nothing,AbstractString} = nothing,
)
    backend = getbackend()
    cpu_marker = backend isa KernelAbstractions.CPU
    run_id = fom_id(cfg; tag = run_tag, cpu_marker = cpu_marker)
    dir = fom_dir(run_id)
    stream_paths = snapshot_stream_paths(dir)
    forced_backend = lowercase(get(ENV, "RAYLEIGH_BACKEND", ""))

    if isdir(dir)
        if skip_if_exists
            if !has_streamed_snapshots(dir)
                @warn "FOM run directory exists but snapshots missing; rerunning" run_id dir
            else
                @info "FOM run already exists, reusing $run_id"
                return run_id
            end
        elseif !force
            error("FOM directory already exists: $dir")
        else
            rm(dir; recursive = true, force = true)
        end
    end

    ensure_run_dirs!(dir)
    save_fom_config(joinpath(dir, "config.toml"), cfg)

    ensure_cuda_psolver!(backend)
    setup = getsetup(cfg; backend = backend)
    psolver = rayleigh_psolver(cfg, setup)
    start = initial_state(cfg, setup, psolver)
    len_state = state_vector_length(start)

    @info "Estimated snapshot memory = $(round(snapshot_memory_gb(len_state, cfg.nsnapshot); digits = 3)) GB"

    binfo = backend_summary()
    @info "Selected backend" (backend = binfo.backend, device = binfo.device, selection = isempty(forced_backend) ? "auto" : forced_backend)
    if binfo.backend == "cpu" && !(forced_backend in ("cpu", "cpu-only", "cpu_only"))
        cuda_ok = false
        try
            cuda_ok = CUDA.functional()
        catch
            cuda_ok = false
        end
        if !cuda_ok
            @error "CUDA backend not available; running on CPU" hint = "Check `nvidia-smi` and `julia --project=. -e 'using CUDA; CUDA.versioninfo()'`. You can force CPU with `RAYLEIGH_BACKEND=cpu`."
        end
    end
    elapsed = @elapsed lens = fullrun_streaming(cfg, start, setup, psolver, stream_paths)

    snapshot_meta = (
        len_u = lens.len_u,
        len_temp = lens.len_temp,
        len_state = lens.len_u + lens.len_temp,
        nsnapshot = cfg.nsnapshot,
        eltype = string(Float64),
        layout = "column_major",
    )
    save_object(stream_paths.meta, snapshot_meta)

    metadata = (
        created_at = string(Dates.now()),
        runtime_seconds = elapsed,
        backend = binfo.backend,
        device = binfo.device,
        snapshot_memory_gb = snapshot_memory_gb(len_state, cfg.nsnapshot),
        snapshot_storage = "streamed",
        run_tag = run_tag,
    )
    save_object(joinpath(dir, "metadata.jld2"), metadata)

    register_fom!(run_id, cfg)
    @info "Saved FOM run to $dir"

    return run_id
end

function rayleigh_psolver(cfg::FOMConfig, setup)
    solver_override = lowercase(strip(get(ENV, "RAYLEIGH_PSOLVER", "")))
    if solver_override == "default"
        @info "Using default Poisson solver" source = "RAYLEIGH_PSOLVER"
        return NS.default_psolver(setup)
    elseif solver_override == "transform"
        @info "Using transform Poisson solver" source = "RAYLEIGH_PSOLVER"
        return NS.psolver_transform(setup)
    elseif !isempty(solver_override)
        error("Unsupported RAYLEIGH_PSOLVER=$(solver_override). Use \"default\" or \"transform\".")
    end

    if isnothing(cfg.stretch_x) && isnothing(cfg.stretch_y)
        # Match RayleighBenard2D.jl: use FFT/DCT Poisson solver on uniform grids.
        try
            @info "Using transform Poisson solver (uniform grid)"
            return NS.psolver_transform(setup)
        catch err
            @warn "Transform Poisson solver unavailable; falling back to default" err
            return NS.default_psolver(setup)
        end
    end
    return NS.default_psolver(setup)
end

function load_streamed_snapshots(dir::AbstractString)
    stream_paths = snapshot_stream_paths(dir)
    has_streamed_snapshots(dir) || error("Streamed snapshots not found in: $dir")

    meta = load_object(stream_paths.meta)
    len_u = Int(meta.len_u)
    len_temp = Int(meta.len_temp)
    len_state = Int(meta.len_state)
    nsnapshot = Int(meta.nsnapshot)
    K = nsnapshot + 1

    T = meta.eltype == "Float64" ? Float64 : error("Unsupported snapshot eltype: $(meta.eltype)")
    state = Mmap.mmap(stream_paths.state, Matrix{T}, (len_state, K))
    dstate = Mmap.mmap(stream_paths.dstate, Matrix{T}, (len_state, K))

    u = @view state[1:len_u, :]
    temp = @view state[(len_u + 1):(len_u + len_temp), :]
    du = @view dstate[1:len_u, :]
    dtemp = @view dstate[(len_u + 1):(len_u + len_temp), :]

    return (; state, dstate, u, temp, du, dtemp, len_u, len_temp)
end

function one_step_fom_from_state(
    x0::AbstractVector,
    setup,
    psolver,
    cfg::FOMConfig;
    dt::Real,
    method::Symbol = :explicit_euler,
)
    T = eltype(x0)
    state = initial_state(cfg, setup, psolver)
    copy_vector_to_state!(state, x0)
    # `timestep!` for explicit Runge-Kutta methods evaluates the force before
    # (re-)applying boundary conditions, so the state must already satisfy BCs.
    NS.apply_bc_u!(state.u, 0.0, setup)
    NS.apply_bc_temp!(state.temp, 0.0, setup)

    params = boussinesq_params(cfg)

    if method == :explicit_euler
        rk = NS.RKMethods.FE11(; T = T)
        ode_cache = NS.get_cache(rk, state, setup)
        force_cache = NS.get_cache(NS.boussinesq!, setup)
        stepper = NS.create_stepper(rk; setup, psolver, state, t = 0.0)
        stepper = NS.timestep!(
            rk,
            NS.boussinesq!,
            stepper,
            dt;
            params,
            ode_cache,
            force_cache,
        )
        x1 = similar(x0)
        copy_state_to_vector!(x1, stepper.state)
        return x1, (method = "explicit_euler", is_explicit_euler = true)
    elseif method == :lmwray3
        rk = NS.LMWray3(; T = T)
        ode_cache = NS.get_cache(rk, state, setup)
        force_cache = NS.get_cache(NS.boussinesq!, setup)
        stepper = NS.create_stepper(rk; setup, psolver, state, t = 0.0)
        stepper = NS.timestep!(
            rk,
            NS.boussinesq!,
            stepper,
            dt;
            params,
            ode_cache,
            force_cache,
        )
        x1 = similar(x0)
        copy_state_to_vector!(x1, stepper.state)
        return x1, (method = "lmwray3", is_explicit_euler = false)
    else
        error("Unsupported one-step FOM method: $(method). Use :explicit_euler or :lmwray3.")
    end
end

function load_fom_snapshots(run_id::String)
    dir = fom_dir(run_id)
    if has_streamed_snapshots(dir)
        return load_streamed_snapshots(dir)
    end
    error("Streamed snapshots not found for $run_id")
end
load_fom_metadata(run_id::String) = load_object(joinpath(fom_dir(run_id), "metadata.jld2"))
load_fom_config_from_run(run_id::String) = load_fom_config(joinpath(fom_dir(run_id), "config.toml"))
