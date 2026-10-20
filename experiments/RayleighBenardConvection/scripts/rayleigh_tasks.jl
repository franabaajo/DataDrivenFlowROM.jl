include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Config.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Paths.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Registry.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Setup.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "FOM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "POD.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "ROM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Plotting.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Stats.jl"))

function usage()
    println("Usage:")
    println("  julia --project=. experiments/RayleighBenardConvection/scripts/rayleigh_tasks.jl run-fom <fom.toml>")
    println("  julia --project=. experiments/RayleighBenardConvection/scripts/rayleigh_tasks.jl sweep-fom [grid_study.toml] [fom.toml]")
    println("  julia --project=. experiments/RayleighBenardConvection/scripts/rayleigh_tasks.jl diagnostics-fom [fom_run_id ...]")
    println("  julia --project=. experiments/RayleighBenardConvection/scripts/rayleigh_tasks.jl postprocess-fom [flags]")
    println("  julia --project=. experiments/RayleighBenardConvection/scripts/rayleigh_tasks.jl postprocess-one <fom_run_id> [flags]")
    println("  julia --project=. experiments/RayleighBenardConvection/scripts/rayleigh_tasks.jl rerun-fom <fom_run_id> <tag>")
    println("  julia --project=. experiments/RayleighBenardConvection/scripts/rayleigh_tasks.jl run-rom <fom_run_id> <rom.toml>")
    println("  julia --project=. experiments/RayleighBenardConvection/scripts/rayleigh_tasks.jl run-switched-rom <fom_run_id> <switched_rom.toml> [--force]")
    println("  julia --project=. experiments/RayleighBenardConvection/scripts/rayleigh_tasks.jl run-auto-switched-rom <fom_run_id> <auto_switched_rom.toml> [--force]")
    println("  julia --project=. experiments/RayleighBenardConvection/scripts/rayleigh_tasks.jl benchmark-advance <fom_run_id> <switched_rom.toml>")
    println("  julia --project=. experiments/RayleighBenardConvection/scripts/rayleigh_tasks.jl benchmark-rom-advance <fom_run_id> <switched_rom.toml>")
    println("")
    println("Postprocess flags:")
    println("  --force-all")
    println("  --force-temperature")
    println("  --force-vorticity")
    println("  --force-singular-values")
    println("  --force-pod-modes")
    println("  --force-nusselt")
    println("  --force-movie")
    println("")
    println("Run flags:")
    println("  --force    Recompute an existing ROM run id instead of reusing saved output")
end

default_grid_study_path() = joinpath(@__DIR__, "..", "configs", "grid_study.toml")
default_fom_base_path() = joinpath(@__DIR__, "..", "configs", "fom_base.toml")

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

function limits_from_toml(d::AbstractDict, key::AbstractString, default)
    if !haskey(d, key)
        return default
    end
    vals = d[key]
    length(vals) == 2 || error("$key must contain exactly two values")
    return (Float64(vals[1]), Float64(vals[2]))
end

function grid_pairs_from_toml(sweep::AbstractDict, basecfg::FOMConfig)
    if haskey(sweep, "n")
        if haskey(sweep, "nx") || haskey(sweep, "ny")
            error("Use either square-grid key `n` or rectangular-grid keys `nx`/`ny`, not both.")
        end
        n_list = Int.(list_from_toml(sweep, "n", [basecfg.nx]))
        return [(n, n) for n in n_list]
    end

    nx_list = Int.(list_from_toml(sweep, "nx", [basecfg.nx]))
    ny_list = Int.(list_from_toml(sweep, "ny", [basecfg.ny]))
    return [(nx, ny) for nx in nx_list for ny in ny_list]
end

function stretch_pair_lists_from_toml(sweep::AbstractDict, basecfg::FOMConfig)
    if haskey(sweep, "stretch")
        if haskey(sweep, "stretch_x") || haskey(sweep, "stretch_y")
            error("Use either shared key `stretch` or directional keys `stretch_x`/`stretch_y`, not both.")
        end
        stretch_list = stretch_from_toml_list(list_from_toml(sweep, "stretch", [basecfg.stretch_x]))
        return stretch_list, stretch_list
    end

    stretch_x_list = stretch_from_toml_list(list_from_toml(sweep, "stretch_x", [basecfg.stretch_x]))
    stretch_y_list = stretch_from_toml_list(list_from_toml(sweep, "stretch_y", [basecfg.stretch_y]))
    return stretch_x_list, stretch_y_list
end

function parse_postprocess_flags(args)
    force = (
        temperature = false,
        vorticity = false,
        singular_values = false,
        pod_modes = false,
        nusselt = false,
        movie = false,
    )
    remaining = String[]
    for arg in args
        if arg == "--force-all"
            force = (
                temperature = true,
                vorticity = true,
                singular_values = true,
                pod_modes = true,
                nusselt = true,
                movie = true,
            )
        elseif arg == "--force-temperature"
            force = (; force..., temperature = true)
        elseif arg == "--force-vorticity"
            force = (; force..., vorticity = true)
        elseif arg == "--force-singular-values"
            force = (; force..., singular_values = true)
        elseif arg == "--force-pod-modes"
            force = (; force..., pod_modes = true)
        elseif arg == "--force-nusselt"
            force = (; force..., nusselt = true)
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
    postprocess_one_fom(run_id)
    println("FOM run id: ", run_id)
    println("Saved in: ", fom_dir(run_id))
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

    grid_pairs = grid_pairs_from_toml(sweep, basecfg)
    xlims = limits_from_toml(sweep, "xlims", basecfg.xlims)
    ylims = limits_from_toml(sweep, "ylims", basecfg.ylims)
    cfl_list = Float64.(list_from_toml(sweep, "cfl", [basecfg.cfl]))
    stretch_x_list, stretch_y_list = stretch_pair_lists_from_toml(sweep, basecfg)
    nsnapshot_list = Int.(list_from_toml(sweep, "nsnapshot", [basecfg.nsnapshot]))
    viscosity_list = Float64.(list_from_toml(sweep, "viscosity", [basecfg.viscosity]))
    conductivity_list = Float64.(list_from_toml(sweep, "conductivity", [basecfg.conductivity]))

    for (nx, ny) in grid_pairs, cfl in cfl_list, stx in stretch_x_list, sty in stretch_y_list,
        ns in nsnapshot_list, nu in viscosity_list, kappa in conductivity_list

        cfg = FOMConfig(
            viscosity = nu,
            conductivity = kappa,
            gravity = basecfg.gravity,
            gdir = basecfg.gdir,
            dodissipation = basecfg.dodissipation,
            nx = nx,
            ny = ny,
            xlims = xlims,
            ylims = ylims,
            stretch_x = stx,
            stretch_y = sty,
            temp_hot = basecfg.temp_hot,
            temp_cold = basecfg.temp_cold,
            initial_temperature_profile = basecfg.initial_temperature_profile,
            perturbation_amplitude = basecfg.perturbation_amplitude,
            perturbation_wavenumber = basecfg.perturbation_wavenumber,
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
        joinpath(dir, "plots", "final_temperature.png"),
        joinpath(dir, "plots", "final_vorticity.png"),
        joinpath(dir, "plots", "singular_values.pdf"),
        joinpath(dir, "plots", "pod_modes.png"),
        joinpath(dir, "plots", "nusselt_number.png"),
        joinpath(dir, "movies", "fom_evolution.mp4"),
    ]
    return all(isfile, needed)
end

function postprocess_fom_outputs(;
    force = (
        temperature = false,
        vorticity = false,
        singular_values = false,
        pod_modes = false,
        nusselt = false,
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
        temperature = false,
        vorticity = false,
        singular_values = false,
        pod_modes = false,
        nusselt = false,
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

function run_rom_from_fom(parent_fom_id::String, rom_path::AbstractString)
    romcfg = load_rom_config(rom_path)
    run_id = run_rom(parent_fom_id, romcfg; skip_if_exists = true)
    plot_rom_outputs(run_id; make_movie = true)
    println("ROM run id: ", run_id)
    println("Saved in: ", rom_dir(run_id))
end

function parse_run_flags(args)
    force = false
    remaining = String[]
    for arg in args
        if arg == "--force"
            force = true
        else
            push!(remaining, arg)
        end
    end
    return force, remaining
end

function run_switched_rom_from_fom(parent_fom_id::String, switch_path::AbstractString; force::Bool = false)
    switchcfg = load_switched_rom_config(switch_path)
    run_id = run_switched_rom(parent_fom_id, switchcfg; force = force, skip_if_exists = !force)
    plot_rom_outputs(run_id; make_movie = true)
    println("Switched ROM run id: ", run_id)
    println("Saved in: ", rom_dir(run_id))
end

function _toml_float(d::AbstractDict, key::AbstractString, default)
    return haskey(d, key) ? Float64(d[key]) : default
end

function _toml_int(d::AbstractDict, key::AbstractString, default)
    return haskey(d, key) ? Int(d[key]) : default
end

function _toml_bool(d::AbstractDict, key::AbstractString, default)
    return haskey(d, key) ? Bool(d[key]) : default
end

function _toml_string(d::AbstractDict, key::AbstractString, default)
    return haskey(d, key) ? String(d[key]) : default
end

function _toml_optional_float(d::AbstractDict, key::AbstractString, default)
    return haskey(d, key) ? Float64(d[key]) : default
end

function _first_downcrossing_time(
    times::AbstractVector,
    values::AbstractVector,
    threshold::Real;
    start_time::Real,
)
    istart = searchsortedfirst(times, Float64(start_time))
    istart = clamp(istart, 1, length(times))

    if values[istart] <= threshold
        return Float64(times[istart])
    end

    for i in (istart + 1):length(times)
        if values[i] <= threshold
            t0 = Float64(times[i - 1])
            t1 = Float64(times[i])
            y0 = Float64(values[i - 1])
            y1 = Float64(values[i])
            if y1 == y0
                return t1
            end
            α = clamp((Float64(threshold) - y0) / (y1 - y0), 0.0, 1.0)
            return t0 + α * (t1 - t0)
        end
    end

    error("No downward crossing found for threshold $threshold after t = $start_time")
end

function _auto_switched_local_windows(parent_fom_id::String, cfgdict::AbstractDict)
    fomcfg = load_fom_config_from_run(parent_fom_id)
    start_time = _toml_float(cfgdict, "start_time", 0.0)
    end_time = _toml_optional_float(cfgdict, "end_time", fomcfg.tsim)
    min_switch_time = _toml_float(cfgdict, "min_switch_time", 20.0)
    pre_to_transition_threshold = _toml_float(cfgdict, "pre_to_transition_threshold", 0.15)
    transition_to_post_threshold = _toml_float(cfgdict, "transition_to_post_threshold", -0.50)
    padding = _toml_float(cfgdict, "window_padding", 3.0)

    data = compute_fom_diagnostics(
        parent_fom_id;
        diagnostics = (:angular_momentum_z_normalized,),
    )
    times = collect(data.times)
    lz = data.results[:angular_momentum_z_normalized]

    t_pre_to_transition = _first_downcrossing_time(
        times,
        lz,
        pre_to_transition_threshold;
        start_time = max(start_time, min_switch_time),
    )
    t_transition_to_post = _first_downcrossing_time(
        times,
        lz,
        transition_to_post_threshold;
        start_time = t_pre_to_transition,
    )
    t_transition_to_post > t_pre_to_transition ||
        error("Detected transition_to_post crossing before pre_to_transition crossing")

    pre_start = start_time
    pre_end = min(end_time, t_pre_to_transition + padding)
    transition_start = max(start_time, t_pre_to_transition - padding)
    transition_end = min(end_time, t_transition_to_post + padding)
    post_start = max(start_time, t_transition_to_post - padding)
    post_end = end_time

    return (;
        start_time,
        end_time,
        min_switch_time,
        pre_to_transition_threshold,
        transition_to_post_threshold,
        t_pre_to_transition,
        t_transition_to_post,
        padding,
        windows = (
            pre = (pre_start, pre_end),
            transition = (transition_start, transition_end),
            post = (post_start, post_end),
        ),
    )
end

function _auto_local_rom_config(cfgdict::AbstractDict, tag::String, window)
    return ROMConfig(
        run_tag = tag,
        npod = _toml_int(cfgdict, "npod", 120),
        rank_selection = _toml_string(cfgdict, "rank_selection", "energy_fraction"),
        basis_centered = _toml_bool(cfgdict, "basis_centered", true),
        lambda = _toml_float(cfgdict, "lambda", 0.01),
        safety = _toml_float(cfgdict, "safety", 0.03),
        growthfac = _toml_float(cfgdict, "growthfac", 1.05),
        energy_fraction = _toml_float(cfgdict, "energy_fraction", 0.99),
        training_start_time = window[1],
        training_end_time = window[2],
        use_intrusive_rom = _toml_bool(cfgdict, "use_intrusive_rom", false),
        use_exact_opinf = _toml_bool(cfgdict, "use_exact_opinf", true),
        exact_opinf_force_explicit_euler = _toml_bool(cfgdict, "exact_opinf_force_explicit_euler", true),
        fallback_to_standard_opinf = _toml_bool(cfgdict, "fallback_to_standard_opinf", true),
        include_constant_term = _toml_bool(cfgdict, "include_constant_term", true),
        basis_mode = _toml_string(cfgdict, "basis_mode", "global"),
        basis_local_strategy = _toml_string(cfgdict, "basis_local_strategy", "S3"),
        basis_partition_method = _toml_string(cfgdict, "basis_partition_method", "uniform_x"),
        basis_partition_nx = _toml_int(cfgdict, "basis_partition_nx", 1),
        basis_partition_ny = _toml_int(cfgdict, "basis_partition_ny", 1),
        basis_local_energy_fraction = _toml_float(cfgdict, "basis_local_energy_fraction", 0.99),
        basis_space_local_sparse = _toml_bool(cfgdict, "basis_space_local_sparse", true),
    )
end

function run_auto_switched_rom_from_fom(parent_fom_id::String, auto_path::AbstractString; force::Bool = false)
    cfgdict = TOML.parsefile(auto_path)
    windows = _auto_switched_local_windows(parent_fom_id, cfgdict)
    tag_prefix = _toml_string(cfgdict, "run_tag", "auto_lz_state_switch")

    local_specs = [
        ("pre", windows.windows.pre),
        ("transition", windows.windows.transition),
        ("post", windows.windows.post),
    ]
    local_rom_ids = String[]
    for (name, window) in local_specs
        tag = "$(tag_prefix)_local_$(name)_$(fmt_token(window[1]))_$(fmt_token(window[2]))"
        romcfg = _auto_local_rom_config(cfgdict, tag, window)
        @info "Training auto-window local ROM" name window tag
        run_id = run_rom(parent_fom_id, romcfg; force = force, skip_if_exists = !force)
        push!(local_rom_ids, run_id)
        plot_rom_outputs(run_id; make_movie = false)
    end

    switching_policy = _toml_string(cfgdict, "switching_policy", "scheduled")
    switching_policy == "scheduled" ||
        error("switching_policy=\"$switching_policy\" is disabled; auto switched ROMs now use scheduled switching only")
    pre_to_transition_time = haskey(cfgdict, "pre_to_transition_time") ?
        Float64(cfgdict["pre_to_transition_time"]) :
        windows.t_pre_to_transition
    transition_to_post_time = haskey(cfgdict, "transition_to_post_time") ?
        Float64(cfgdict["transition_to_post_time"]) :
        windows.t_transition_to_post

    switchcfg = SwitchedROMConfig(
        run_tag = tag_prefix,
        local_rom_ids = local_rom_ids,
        regime_names = ["pre", "transition", "post"],
        start_time = windows.start_time,
        end_time = windows.end_time,
        initial_regime = _toml_string(cfgdict, "initial_regime", "pre"),
        switching_policy = switching_policy,
        min_switch_time = windows.min_switch_time,
        pre_to_transition_threshold = windows.pre_to_transition_threshold,
        transition_to_post_threshold = windows.transition_to_post_threshold,
        transition_to_pre_threshold = _toml_float(cfgdict, "transition_to_pre_threshold", 0.20),
        post_to_transition_threshold = _toml_float(cfgdict, "post_to_transition_threshold", -0.35),
        allow_backward_switches = false,
        pre_to_transition_time = pre_to_transition_time,
        transition_to_post_time = transition_to_post_time,
        safety = _toml_float(cfgdict, "safety", 0.03),
        growthfac = _toml_float(cfgdict, "growthfac", 1.05),
    )

    @info "Running auto-window switched ROM" (
        t_pre_to_transition = windows.t_pre_to_transition,
        t_transition_to_post = windows.t_transition_to_post,
        padding = windows.padding,
        local_rom_ids = local_rom_ids,
    )
    run_id = run_switched_rom(parent_fom_id, switchcfg; force = force, skip_if_exists = !force)
    plot_rom_outputs(run_id; make_movie = true)

    println("Detected FOM threshold crossings:")
    println("  pre -> transition: ", windows.t_pre_to_transition)
    println("  transition -> post: ", windows.t_transition_to_post)
    println("Local training windows:")
    println("  pre: ", windows.windows.pre)
    println("  transition: ", windows.windows.transition)
    println("  post: ", windows.windows.post)
    println("Local ROM run ids:")
    for run_id_local in local_rom_ids
        println("  ", run_id_local)
    end
    println("Switched ROM run id: ", run_id)
    println("Saved in: ", rom_dir(run_id))
end

function benchmark_fom_advance_cpu(parent_fom_id::String, switchcfg::SwitchedROMConfig)
    fomcfg = load_fom_config_from_run(parent_fom_id)
    snapshots = load_fom_snapshots(parent_fom_id)
    Kfull = size(snapshots.state, 2)
    prediction_window = rom_training_window(fomcfg, Kfull, switchcfg.start_time, switchcfg.end_time)

    setup = getsetup(fomcfg; backend = KernelAbstractions.CPU())
    psolver = rayleigh_psolver(fomcfg, setup)
    state = initial_state(fomcfg, setup, psolver)
    copy_vector_to_state!(state, @view(snapshots.state[:, first(prediction_window.snapshot_indices)]))
    NS.apply_bc_u!(state.u, prediction_window.actual_start_time, setup)
    NS.apply_bc_temp!(state.temp, prediction_window.actual_start_time, setup)

    method = NS.LMWray3()
    ode_cache = NS.get_cache(method, state, setup)
    force_cache = NS.get_cache(NS.boussinesq!, setup)
    stepper = NS.create_stepper(method; setup, psolver, state, t = prediction_window.actual_start_time)
    params = boussinesq_params(fomcfg)

    nsteps = Ref(0)
    elapsed = @elapsed begin
        while stepper.t < prevfloat(prediction_window.end_time)
            dt = fomcfg.cfl * NS.propose_timestep(NS.boussinesq!, stepper.state, setup, params)
            dt = min(dt, prediction_window.end_time - stepper.t)
            stepper = NS.timestep!(
                method,
                NS.boussinesq!,
                stepper,
                dt;
                params,
                ode_cache,
                force_cache,
            )
            nsteps[] += 1
        end
    end

    return (;
        runtime_seconds = elapsed,
        nsteps = nsteps[],
        start_time = prediction_window.actual_start_time,
        end_time = prediction_window.end_time,
        duration = prediction_window.duration,
        backend = "cpu",
    )
end

function benchmark_switched_rom_advance_cpu(parent_fom_id::String, switchcfg::SwitchedROMConfig)
    validate_switched_rom_config(switchcfg)
    fomcfg = load_fom_config_from_run(parent_fom_id)
    snapshots = load_fom_snapshots(parent_fom_id)
    Kfull = size(snapshots.state, 2)
    prediction_window = rom_training_window(fomcfg, Kfull, switchcfg.start_time, switchcfg.end_time)
    prediction_indices = prediction_window.snapshot_indices
    prediction_times = collect(prediction_window.times)

    models = [
        _load_local_rom_model(run_id, switchcfg.regime_names[i])
        for (i, run_id) in enumerate(switchcfg.local_rom_ids)
    ]
    for model in models
        local_meta = load_rom_metadata(model.run_id)
        _check_local_rom_parent_compatible!(model, local_meta, parent_fom_id, fomcfg, switchcfg)
    end

    initial_active = findfirst(==(switchcfg.initial_regime), switchcfg.regime_names)
    isnothing(initial_active) && error("Could not find initial_regime=$(switchcfg.initial_regime)")
    active = Int(initial_active)
    a = _project_state_to_local(models[active], @view(snapshots.state[:, first(prediction_indices)]))

    setup = getsetup(fomcfg; backend = KernelAbstractions.CPU())
    psolver = NS.default_psolver(setup)
    state = initial_state(fomcfg, setup, psolver)
    angular_cache = angular_momentum_z_cache(setup, state.u, fomcfg)
    state_weights = state_mass_vector(setup, state.u, state.temp)
    xbuf = Vector{Float64}(undef, size(snapshots.state, 1))
    xproj = similar(xbuf)

    switch_events = NamedTuple[]
    nintervals = 0
    advance_runtime_seconds = 0.0
    elapsed = @elapsed begin
        for k in 1:(length(prediction_times) - 1)
            t0 = prediction_times[k]
            t1 = prediction_times[k + 1]
            model = models[active]
            advance_runtime_seconds += @elapsed _advance_local_rom_interval!(
                a,
                model,
                t1 - t0;
                safety = switchcfg.safety,
                growthfac = switchcfg.growthfac,
                absolute_time = t0,
            )
            nintervals += 1

            _reconstruct_local_state!(xbuf, model, a)
            lz_before = _normalized_lz_from_state_vector!(xbuf, state, setup, angular_cache, t1)
            target = _switched_target_index(active, lz_before, t1, switchcfg)
            if target != active
                old_active = active
                new_model = models[target]
                anew = _project_state_to_local(new_model, xbuf)
                _reconstruct_local_state!(xproj, new_model, anew)
                lz_after = _normalized_lz_from_state_vector!(xproj, state, setup, angular_cache, t1)
                jump_rel = _weighted_relative_state_difference(xproj, xbuf, state_weights)
                push!(
                    switch_events,
                    (
                        time = t1,
                        from_index = old_active,
                        to_index = target,
                        from_regime = switchcfg.regime_names[old_active],
                        to_regime = switchcfg.regime_names[target],
                        lz_before_projection = lz_before,
                        lz_after_projection = lz_after,
                        projection_jump_rel = jump_rel,
                    ),
                )
                active = target
                a = anew
            end
        end
    end

    return (;
        runtime_seconds = elapsed,
        advance_runtime_seconds,
        overhead_seconds = elapsed - advance_runtime_seconds,
        nintervals,
        start_time = prediction_window.actual_start_time,
        end_time = prediction_window.end_time,
        duration = prediction_window.duration,
        local_rfits = [model.rfit for model in models],
        rfit = maximum(model.rfit for model in models),
        nswitches = length(switch_events),
        switch_events,
        backend = "cpu",
    )
end

function benchmark_advance(parent_fom_id::String, switch_path::AbstractString)
    switchcfg = load_switched_rom_config(switch_path)
    @info "Benchmarking CPU FOM advance without snapshot output"
    fom = benchmark_fom_advance_cpu(parent_fom_id, switchcfg)
    @info "Benchmarking CPU switched ROM advance without trajectory diagnostics"
    rom = benchmark_switched_rom_advance_cpu(parent_fom_id, switchcfg)

    result = (;
        created_at = string(Dates.now()),
        parent_fom_id,
        switch_path = String(switch_path),
        fom,
        rom,
        speedup_rom_runtime = fom.runtime_seconds / rom.runtime_seconds,
        speedup_rom_advance_only = fom.runtime_seconds / rom.advance_runtime_seconds,
    )

    outdir = joinpath(output_root(), "diagnostics", "advance_benchmarks")
    mkpath(outdir)
    outpath = joinpath(outdir, "benchmark_$(short_hash(parent_fom_id * switch_path * string(Dates.now()))).jld2")
    save_object(outpath, result)

    println("CPU advance benchmark:")
    println("  interval: ", fom.start_time, "..", fom.end_time, " (duration ", fom.duration, ")")
    println("  FOM runtime_seconds: ", fom.runtime_seconds)
    println("  FOM time steps: ", fom.nsteps)
    println("  ROM runtime_seconds: ", rom.runtime_seconds)
    println("  ROM advance_runtime_seconds: ", rom.advance_runtime_seconds)
    println("  ROM overhead_seconds: ", rom.overhead_seconds)
    println("  ROM local_rfits: ", rom.local_rfits)
    println("  ROM switches: ", rom.nswitches)
    println("  speedup FOM/ROM runtime: ", result.speedup_rom_runtime)
    println("  speedup FOM/ROM advance-only: ", result.speedup_rom_advance_only)
    println("  saved: ", outpath)
end

function benchmark_rom_advance(parent_fom_id::String, switch_path::AbstractString)
    switchcfg = load_switched_rom_config(switch_path)
    @info "Benchmarking CPU switched ROM advance without trajectory diagnostics"
    rom = benchmark_switched_rom_advance_cpu(parent_fom_id, switchcfg)

    result = (;
        created_at = string(Dates.now()),
        parent_fom_id,
        switch_path = String(switch_path),
        rom,
    )

    outdir = joinpath(output_root(), "diagnostics", "advance_benchmarks")
    mkpath(outdir)
    outpath = joinpath(outdir, "rom_benchmark_$(short_hash(parent_fom_id * switch_path * string(Dates.now()))).jld2")
    save_object(outpath, result)

    println("CPU ROM advance benchmark:")
    println("  interval: ", rom.start_time, "..", rom.end_time, " (duration ", rom.duration, ")")
    println("  ROM runtime_seconds: ", rom.runtime_seconds)
    println("  ROM advance_runtime_seconds: ", rom.advance_runtime_seconds)
    println("  ROM overhead_seconds: ", rom.overhead_seconds)
    println("  ROM local_rfits: ", rom.local_rfits)
    println("  ROM switches: ", rom.nswitches)
    println("  saved: ", outpath)
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

if cmd == "sweep-fom" || cmd == "sweep"
    if length(rest) > 2
        usage()
        error("Unexpected arguments: $(join(rest[3:end], " "))")
    end
    grid_path = length(rest) >= 1 ? rest[1] : default_grid_study_path()
    fom_path = length(rest) >= 2 ? rest[2] : default_fom_base_path()
    run_sweep(grid_path, fom_path)
elseif cmd == "run-fom"
    if isempty(rest)
        usage()
        error("Missing <fom_base.toml>")
    end
    run_fom_from_config(rest[1])
elseif cmd == "diagnostics-fom"
    if isempty(rest)
        run_all_fom_diagnostics()
    elseif length(rest) == 1
        run_fom_diagnostics(rest[1])
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
    parent_fom_id = rest[1]
    run_rom_from_fom(parent_fom_id, rest[2])
elseif cmd == "run-switched-rom"
    force, remaining = parse_run_flags(rest)
    if length(remaining) < 2
        usage()
        error("Missing <fom_run_id> <switched_rom.toml>")
    end
    if length(remaining) > 2
        error("Unexpected arguments: $(join(remaining[3:end], " "))")
    end
    parent_fom_id = remaining[1]
    run_switched_rom_from_fom(parent_fom_id, remaining[2]; force = force)
elseif cmd == "run-auto-switched-rom"
    force, remaining = parse_run_flags(rest)
    if length(remaining) < 2
        usage()
        error("Missing <fom_run_id> <auto_switched_rom.toml>")
    end
    if length(remaining) > 2
        error("Unexpected arguments: $(join(remaining[3:end], " "))")
    end
    parent_fom_id = remaining[1]
    run_auto_switched_rom_from_fom(parent_fom_id, remaining[2]; force = force)
elseif cmd == "benchmark-advance"
    if length(rest) < 2
        usage()
        error("Missing <fom_run_id> <switched_rom.toml>")
    end
    if length(rest) > 2
        error("Unexpected arguments: $(join(rest[3:end], " "))")
    end
    benchmark_advance(rest[1], rest[2])
elseif cmd == "benchmark-rom-advance"
    if length(rest) < 2
        usage()
        error("Missing <fom_run_id> <switched_rom.toml>")
    end
    if length(rest) > 2
        error("Unexpected arguments: $(join(rest[3:end], " "))")
    end
    benchmark_rom_advance(rest[1], rest[2])
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
