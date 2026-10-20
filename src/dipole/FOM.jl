using Dates
using JLD2
using KernelAbstractions
import IncompressibleNavierStokes as NS

include(joinpath(@__DIR__, "..", "workflows", "FOM.jl"))

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

    @info "Estimated snapshot memory = $(round(snapshot_memory_gb(cfg); digits = 3)) GB"

    ensure_cuda_psolver!(backend)
    setup = getsetup(cfg; backend = backend)
    psolver = NS.default_psolver(setup)
    u = initial_velocity(cfg, setup, psolver)
    scalevelocity!(u, setup, cfg)

    binfo = backend_summary()
    elapsed = @elapsed fullrun_streaming(cfg, u, setup, psolver, stream_paths)

    snapshot_meta = (
        len = length(u),
        nsnapshot = cfg.nsnapshot,
        eltype = string(Float64),
        layout = "column_major",
    )
    save_object(stream_paths.meta, snapshot_meta; iotype = WORKFLOW_JLD2_IOTYPE)

    metadata = (
        created_at = string(Dates.now()),
        runtime_seconds = elapsed,
        backend = binfo.backend,
        device = binfo.device,
        snapshot_memory_gb = snapshot_memory_gb(cfg),
        snapshot_storage = "streamed",
        run_tag = run_tag,
    )
    save_object(joinpath(dir, "metadata.jld2"), metadata; iotype = WORKFLOW_JLD2_IOTYPE)

    register_fom!(run_id, cfg)
    @info "Saved FOM run to $dir"

    return run_id
end

function load_fom_snapshots(run_id::String)
    dir = fom_dir(run_id)
    if has_streamed_snapshots(dir)
        return load_streamed_snapshots(dir)
    end
    error("Streamed snapshots not found for $run_id")
end
load_fom_metadata(run_id::String) = workflow_load_object_iostream(joinpath(fom_dir(run_id), "metadata.jld2"))
load_fom_config_from_run(run_id::String) = load_fom_config(joinpath(fom_dir(run_id), "config.toml"))
