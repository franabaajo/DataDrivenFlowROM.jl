include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Config.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Paths.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Registry.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Setup.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "FOM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "POD.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "ROM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Plotting.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Stats.jl"))

function usage()
    println("Usage:")
    println("  julia --project=. experiments/DipoleWallCollision/scripts/dipole_tasks.jl run-all <fom.toml> <rom.toml>")
    println("  julia --project=. experiments/DipoleWallCollision/scripts/dipole_tasks.jl run-fom <fom.toml>")
    println("  julia --project=. experiments/DipoleWallCollision/scripts/dipole_tasks.jl sweep-fom <grid_study.toml> <fom.toml>")
    println("  julia --project=. experiments/DipoleWallCollision/scripts/dipole_tasks.jl diagnostics-fom [fom_run_id ...]")
    println("  julia --project=. experiments/DipoleWallCollision/scripts/dipole_tasks.jl postprocess-fom [flags]")
    println("  julia --project=. experiments/DipoleWallCollision/scripts/dipole_tasks.jl postprocess-one <fom_run_id> [flags]")
    println("  julia --project=. experiments/DipoleWallCollision/scripts/dipole_tasks.jl plot-10-modes <fom_run_id>")
    println("  julia --project=. experiments/DipoleWallCollision/scripts/dipole_tasks.jl rerun-fom <fom_run_id> <tag>")
    println("  julia --project=. experiments/DipoleWallCollision/scripts/dipole_tasks.jl run-rom <fom_run_id> <rom.toml> [--no-plots] [--no-movie]")
    println("  julia --project=. experiments/DipoleWallCollision/scripts/dipole_tasks.jl rerun-rom <rom_run_id>")
    println("")
    println("Postprocess flags:")
    println("  --force-all")
    println("  --force-grid")
    println("  --force-vorticity")
    println("  --force-vorticity-zoom")
    println("  --force-pressure-zoom")
    println("  --force-singular-values")
    println("  --force-pod-modes")
    println("  --force-movie")
end

function parse_run_rom_flags(args)
    length(args) >= 2 || error("Missing <fom_run_id> <rom_base.toml>")
    parent_fom_id = args[1]
    rom_path = args[2]
    make_plots = true
    make_movie = true
    for arg in args[3:end]
        if arg == "--no-plots"
            make_plots = false
            make_movie = false
        elseif arg == "--no-movie"
            make_movie = false
        else
            error("Unexpected run-rom argument: $arg")
        end
    end
    return parent_fom_id, rom_path, make_plots, make_movie
end

function list_from_toml(d::AbstractDict, key::AbstractString, default)
    if !haskey(d, key)
        return default
    end
    v = d[key]
    return v isa AbstractVector ? v : [v]
end

function stretch_from_toml_list(vals)
    return [_stretch_from_toml(v) for v in vals]
end

function parse_postprocess_flags(args)
    force = (
        grid = false,
        vorticity = false,
        vorticity_zoom = false,
        pressure_zoom = false,
        singular_values = false,
        pod_modes = false,
        movie = false,
    )
    remaining = String[]
    for arg in args
        if arg == "--force-all"
            force = (
                grid = true,
                vorticity = true,
                vorticity_zoom = true,
                pressure_zoom = true,
                singular_values = true,
                pod_modes = true,
                movie = true,
            )
        elseif arg == "--force-grid"
            force = (; force..., grid = true)
        elseif arg == "--force-vorticity"
            force = (; force..., vorticity = true)
        elseif arg == "--force-vorticity-zoom"
            force = (; force..., vorticity_zoom = true)
        elseif arg == "--force-pressure-zoom"
            force = (; force..., pressure_zoom = true)
        elseif arg == "--force-singular-values"
            force = (; force..., singular_values = true)
        elseif arg == "--force-pod-modes"
            force = (; force..., pod_modes = true)
        elseif arg == "--force-movie"
            force = (; force..., movie = true)
        elseif startswith(arg, "--")
            error("Unknown flag: $arg")
        else
            push!(remaining, arg)
        end
    end
    return force, remaining
end

function run_fom_from_config(fom_path::AbstractString)
    cfg = load_fom_config(fom_path)
    run_id = run_fom(cfg; skip_if_exists = true)
    println("FOM run id: ", run_id)
    println("Saved in: ", fom_dir(run_id))
end

function run_all_from_configs(fom_path::AbstractString, rom_path::AbstractString)
    fomcfg = load_fom_config(fom_path)
    romcfg = load_rom_config(rom_path)

    fom_run_id = run_fom(fomcfg; skip_if_exists = true)
    plot_fom_outputs(fom_run_id; make_movie = true)

    rom_run_id = run_rom(fom_run_id, romcfg; skip_if_exists = true)
    plot_rom_outputs(rom_run_id; make_movie = true)

    println("FOM run id: ", fom_run_id)
    println("ROM run id: ", rom_run_id)
end

function list_fom_run_ids()
    root = fom_root()
    if !isdir(root)
        return String[]
    end
    runs = String[]
    for entry in readdir(root; join = true)
        if isdir(entry) && has_streamed_snapshots(entry)
            push!(runs, basename(entry))
        end
    end
    sort!(runs)
    return runs
end

function run_fom_diagnostics(run_id::String)
    outdir = save_fom_diagnostic_plots(run_id)
    println("Diagnostics saved in: ", outdir)
end

function run_all_fom_diagnostics()
    run_ids = list_fom_run_ids()
    if isempty(run_ids)
        println("No FOM runs found in: ", fom_root())
        return
    end
    outdir = save_fom_diagnostics_compare(run_ids)
    println("Diagnostics saved in: ", outdir)
    for run_id in run_ids
        @info "Postprocessing FOM (movie)" run_id
        postprocess_one_fom(run_id)
    end
end

function run_compare_fom_diagnostics(run_ids::AbstractVector{String})
    outdir = save_fom_diagnostics_compare(run_ids)
    println("Diagnostics saved in: ", outdir)
    for run_id in run_ids
        @info "Postprocessing FOM (movie)" run_id
        postprocess_one_fom(run_id)
    end
end

function run_sweep(grid_path::AbstractString, base_fom_path::AbstractString)
    basecfg = load_fom_config(base_fom_path)
    sweep = TOML.parsefile(grid_path)

    n_list = Int.(list_from_toml(sweep, "n", [basecfg.n]))
    cfl_list = Float64.(list_from_toml(sweep, "cfl", [basecfg.cfl]))
    stretch_list = stretch_from_toml_list(list_from_toml(sweep, "stretch", [basecfg.stretch]))
    nsnapshot_list = Int.(list_from_toml(sweep, "nsnapshot", [basecfg.nsnapshot]))
    viscosity_list = Float64.(list_from_toml(sweep, "viscosity", [basecfg.viscosity]))

    for n in n_list, cfl in cfl_list, stretch in stretch_list, ns in nsnapshot_list, nu in viscosity_list
        cfg = FOMConfig(
            viscosity = nu,
            r0 = basecfg.r0,
            p1 = basecfg.p1,
            p2 = basecfg.p2,
            initialenergy = basecfg.initialenergy,
            n = n,
            stretch = stretch,
            lims = basecfg.lims,
            tsim = basecfg.tsim,
            cfl = cfl,
            nsnapshot = ns,
        )

        backend = getbackend()
        existing = find_existing_fom(cfg; cpu_marker = backend isa KernelAbstractions.CPU)
        if !isnothing(existing)
            @warn "FOM already exists, skipping" run_id = existing
            continue
        end

        run_id = run_fom(cfg; skip_if_exists = true)
        println("Finished: ", run_id)
    end
end

function is_fom_outputs_complete(dir::AbstractString)
    needed = [
        joinpath(dir, "plots", "grid.png"),
        joinpath(dir, "plots", "final_vorticity.png"),
        joinpath(dir, "plots", "final_vorticity_zoom.png"),
        joinpath(dir, "plots", "final_pressure_zoom.png"),
        joinpath(dir, "plots", "singular_values.pdf"),
        joinpath(dir, "plots", "pod_modes.png"),
        joinpath(dir, "movies", "fom_evolution.mp4"),
    ]
    return all(isfile, needed)
end

function postprocess_fom_outputs(;
    force = (
        grid = false,
        vorticity = false,
        vorticity_zoom = false,
        pressure_zoom = false,
        singular_values = false,
        pod_modes = false,
        movie = false,
    ),
)
    root = fom_root()
    if !isdir(root)
        println("No FOM output directory found at: ", root)
        return
    end

    for entry in readdir(root; join = true)
        if !isdir(entry)
            continue
        end
        run_id = basename(entry)
        if !has_streamed_snapshots(entry)
            @warn "Missing snapshots, skipping" run_id entry
            continue
        end

        force_any = any(values(force))
        if is_fom_outputs_complete(entry) && !force_any
            @info "Plots and movie already present, skipping" run_id
            continue
        end

        ensure_run_dirs!(entry)
        movie_path = joinpath(entry, "movies", "fom_evolution.mp4")
        make_movie = !isfile(movie_path) || force.movie

        try
            plot_fom_outputs(run_id; make_movie = make_movie, force = force)
            @info "Postprocessed FOM outputs" run_id
        catch err
            @warn "Failed to postprocess FOM outputs, continuing" run_id err
        end
    end
end

function postprocess_one_fom(
    run_id::String;
    force = (
        grid = false,
        vorticity = false,
        vorticity_zoom = false,
        pressure_zoom = false,
        singular_values = false,
        pod_modes = false,
        movie = false,
    ),
)
    dir = fom_dir(run_id)
    if !isdir(dir)
        error("FOM run directory not found: $dir")
    end
    if !has_streamed_snapshots(dir)
        error("Missing streamed snapshots for: $run_id")
    end

    force_any = any(values(force))
    if is_fom_outputs_complete(dir) && !force_any
        @info "Plots and movie already present, skipping" run_id
        return
    end

    ensure_run_dirs!(dir)
    movie_path = joinpath(dir, "movies", "fom_evolution.mp4")
    make_movie = !isfile(movie_path) || force.movie
    plot_fom_outputs(run_id; make_movie = make_movie, force = force)
    @info "Postprocessed FOM outputs" run_id
end

function run_rom_from_fom(
    parent_fom_id::String,
    rom_path::AbstractString;
    make_plots::Bool = true,
    make_movie::Bool = true,
)
    romcfg = load_rom_config(rom_path)
    run_id = run_rom(parent_fom_id, romcfg; skip_if_exists = true)
    if make_plots
        plot_rom_outputs(run_id; make_movie)
    end
    println("ROM run id: ", run_id)
    println("Saved in: ", rom_dir(run_id))
end

function rerun_rom_by_id(rom_run_id::String)
    rom_run_id = String(strip(rom_run_id))
    dir = rom_dir(rom_run_id)
    if !isdir(dir)
        error("ROM run directory not found: $dir")
    end

    config_path = joinpath(dir, "config.toml")
    if !isfile(config_path)
        error("Missing ROM config.toml in: $dir")
    end
    romcfg = load_rom_config(config_path)

    parent_path = joinpath(dir, "parent_fom_id.txt")
    parent_fom_id = if isfile(parent_path)
        String(strip(read(parent_path, String)))
    else
        meta_path = joinpath(dir, "metadata.jld2")
        if isfile(meta_path)
            meta = load_object_iostream(meta_path)
            hasproperty(meta, :parent_fom_id) || error("ROM metadata missing parent_fom_id: $meta_path")
            string(meta.parent_fom_id)
        else
            error("Missing parent_fom_id.txt and metadata.jld2 in: $dir")
        end
    end
    isempty(parent_fom_id) && error("Empty parent_fom_id for ROM run: $rom_run_id")

    recomputed = rom_id(parent_fom_id, romcfg)
    if recomputed != rom_run_id
        @warn "ROM id mismatch; rerunning with computed id instead of requested id" (
            requested = rom_run_id,
            computed = recomputed,
            dir = dir,
        )
    end

    run_id = run_rom(parent_fom_id, romcfg; force = true, skip_if_exists = false)
    plot_rom_outputs(run_id; make_movie = true)
    println("ROM run id: ", run_id)
    println("Saved in: ", rom_dir(run_id))
end

function run_plot_10_modes(fom_run_id::String)
    plot_path = plot_fom_selected_modes(fom_run_id; modelist = 1:10, plot_name = "pod_modes_10.png")
    println("Saved in: ", plot_path)
end

function rerun_fom_with_tag(parent_fom_id::String, tag::AbstractString)
    cfg = load_fom_config_from_run(parent_fom_id)
    run_id = run_fom(cfg; skip_if_exists = true, run_tag = tag)
    println("FOM run id: ", run_id)
    println("Saved in: ", fom_dir(run_id))
end

if isempty(ARGS)
    usage()
    error("Missing command")
end

cmd = ARGS[1]
rest = ARGS[2:end]

if cmd == "run-all" || cmd == "run-pipeline" || cmd == "pod"
    if length(rest) < 2
        usage()
        error("Missing <fom_base.toml> <rom_base.toml>")
    end
    run_all_from_configs(rest[1], rest[2])
elseif cmd == "sweep-fom" || cmd == "sweep"
    if length(rest) < 2
        usage()
        error("Missing <grid_study.toml> <fom_base.toml>")
    end
    run_sweep(rest[1], rest[2])
elseif cmd == "run-fom"
    if isempty(rest)
        usage()
        error("Missing <fom_base.toml>")
    end
    run_fom_from_config(rest[1])
elseif cmd == "diagnostics-fom"
    if isempty(rest)
        run_all_fom_diagnostics()
    else
        run_compare_fom_diagnostics(rest)
    end
elseif cmd == "postprocess-fom" || cmd == "postprocess"
    force, remaining = parse_postprocess_flags(rest)
    if !isempty(remaining)
        error("Unexpected arguments: $(join(remaining, " "))")
    end
    postprocess_fom_outputs(; force = force)
elseif cmd == "postprocess-one"
    force, remaining = parse_postprocess_flags(rest)
    if isempty(remaining)
        usage()
        error("Missing <fom_run_id>")
    end
    if length(remaining) > 1
        error("Unexpected arguments: $(join(remaining[2:end], " "))")
    end
    postprocess_one_fom(remaining[1]; force = force)
elseif cmd == "run-rom"
    if length(rest) < 2
        usage()
        error("Missing <fom_run_id> <rom_base.toml>")
    end
    parent_fom_id, rom_path, make_plots, make_movie = parse_run_rom_flags(rest)
    run_rom_from_fom(parent_fom_id, rom_path; make_plots, make_movie)
elseif cmd == "rerun-rom"
    if length(rest) != 1
        usage()
        error("Usage: rerun-rom <rom_run_id>")
    end
    rerun_rom_by_id(rest[1])
elseif cmd == "plot-10-modes"
    if length(rest) < 1
        usage()
        error("Missing <fom_run_id>")
    end
    if length(rest) > 1
        error("Unexpected arguments: $(join(rest[2:end], " "))")
    end
    run_plot_10_modes(rest[1])
elseif cmd == "rerun-fom"
    if length(rest) < 2
        usage()
        error("Missing <fom_run_id> <tag>")
    end
    rerun_fom_with_tag(rest[1], rest[2])
else
    usage()
    error("Unknown command: $cmd")
end
