using Adapt
using CairoMakie
using KernelAbstractions
using LinearAlgebra
import IncompressibleNavierStokes as NS

include(joinpath(@__DIR__, "ThermalDiagnostics.jl"))

function copy_state_vector_to_fields!(x::AbstractVector, u, temp)
    nu = length(u)
    nt = length(temp)
    length(x) == nu + nt || error("State vector length mismatch")

    copyto!(u, reshape(view(x, 1:nu), size(u)))
    copyto!(temp, reshape(view(x, (nu + 1):(nu + nt)), size(temp)))
    return u, temp
end

function temperature_coords_and_values(temp, setup)
    coords = getindex.(Array.(setup.xp), setup.Ip.indices)
    values = Array(temp)[setup.Ip]
    return coords, values
end

function save_fom_nusselt_plot(
    plot_path::AbstractString,
    snapshots,
    cfg::FOMConfig,
    setup,
    u,
    temp,
)
    if cfg.dodissipation
        @warn(
            "This FOM includes viscous heating, so statistically averaged top and bottom plate Nusselt numbers need not be equal.",
            conductivity = cfg.conductivity,
            viscosity = cfg.viscosity,
        )
    end

    times = collect(range(0.0, cfg.tsim; length = cfg.nsnapshot + 1))
    nu_bottom = zeros(Float64, length(times))
    nu_top = zeros(Float64, length(times))
    nu_volume = zeros(Float64, length(times))

    for i in eachindex(times)
        copy_state_vector_to_fields!(@view(snapshots.state[:, i]), u, temp)
        NS.apply_bc_u!(u, times[i], setup)
        NS.apply_bc_temp!(temp, times[i], setup)
        nu = plate_nusselt_numbers(temp, setup, cfg)
        nu_bottom[i] = nu.bottom
        nu_top[i] = nu.top
        nu_volume[i] = volume_nusselt_number(u, temp, setup, cfg)
    end

    fig = Figure(size = (1000, 600))
    ax = Axis(
        fig[1, 1];
        title = "Nusselt number",
        xlabel = "t",
        ylabel = "Nu",
    )
    lines!(ax, times, nu_bottom; label = "Bottom plate")
    lines!(ax, times, nu_top; label = "Top plate")
    lines!(ax, times, nu_volume; label = "Volume")
    axislegend(ax; position = :rb)
    save(plot_path, fig)
    return plot_path
end

function cpu_setup_and_fields(cfg::FOMConfig)
    setup = getsetup(cfg; backend = KernelAbstractions.CPU())
    psolver = NS.default_psolver(setup)
    state = initial_state(cfg, setup, psolver)
    return setup, psolver, state.u, state.temp
end

function plot_fom_outputs(
    fom_run_id::String;
    make_movie::Bool = true,
    modelist = [1, 2, 3, 5, 10, 20],
    skip_existing::Bool = true,
    force = (
        temperature = false,
        vorticity = false,
        singular_values = false,
        pod_modes = false,
        nusselt = false,
        movie = false,
    ),
)
    cfg = load_fom_config_from_run(fom_run_id)
    snapshots = load_fom_snapshots(fom_run_id)
    outdir = fom_dir(fom_run_id)

    setup, psolver, u, temp = cpu_setup_and_fields(cfg)
    singular_values_path = joinpath(outdir, "plots", "singular_values.pdf")
    pod_modes_path = joinpath(outdir, "plots", "pod_modes.png")
    need_singular_values = !(skip_existing && isfile(singular_values_path) && !force.singular_values)
    need_pod_modes = !(skip_existing && isfile(pod_modes_path) && !force.pod_modes)
    spectrum = if need_singular_values || need_pod_modes
        w = state_mass_vector(setup, u, temp)
        length(w) == size(snapshots.state, 1) || error("State mass vector length mismatch with snapshot dimension")
        pod_spectrum_from_weighted_gram(snapshots.state, w)
    else
        nothing
    end

    let
        plot_path = joinpath(outdir, "plots", "final_temperature.png")
        if !(skip_existing && isfile(plot_path) && !force.temperature)
            copy_state_vector_to_fields!(@view(snapshots.state[:, end]), u, temp)
            NS.apply_bc_u!(u, cfg.tsim, setup)
            NS.apply_bc_temp!(temp, cfg.tsim, setup)

            fig = Figure()
            ax = Axis(fig[1, 1]; title = "Temperature at final time", xlabel = "x", ylabel = "y", aspect = DataAspect())
            coords, temp_values = temperature_coords_and_values(temp, setup)
            colorrange = NS.get_lims(temp_values, 1.0)
            heatmap!(ax, coords..., temp_values; colormap = :seaborn_icefire_gradient, colorrange)
            save(plot_path, fig)
        end
    end

    let
        plot_path = joinpath(outdir, "plots", "final_vorticity.png")
        if !(skip_existing && isfile(plot_path) && !force.vorticity)
            copy_state_vector_to_fields!(@view(snapshots.state[:, end]), u, temp)
            NS.apply_bc_u!(u, cfg.tsim, setup)

            fig = Figure()
            ax = Axis(fig[1, 1]; title = "Vorticity at final time", xlabel = "x", ylabel = "y", aspect = DataAspect())
            vort = NS.vorticity(u, setup) |> Array
            coords = getindex.(Array.(setup.xp), setup.Ip.indices)
            vort = vort[setup.Ip]
            colorrange = NS.get_lims(vort, 3.0)
            heatmap!(ax, coords..., vort; colormap = :RdBu, colorrange)
            save(plot_path, fig)
        end
    end

    let
        if need_singular_values
            fig = Figure()
            ax = Axis(fig[1, 1]; title = "Singular values", xlabel = "Mode number", ylabel = "sigma_i / sigma_1", yscale = log10)
            if isempty(spectrum.S)
                @warn "No POD singular values available for plotting" fom_run_id
            else
                denom = max(spectrum.S[1], eps(eltype(spectrum.S)))
                scatter!(ax, spectrum.S / denom)
            end
            save(singular_values_path, fig)
        end
    end

    let
        if need_pod_modes
            fig = Figure(size = (1000, 600))
            avail = isempty(spectrum.S) ? 0 : length(spectrum.S)
            filtered = [m for m in modelist if m <= avail]
            if isempty(filtered)
                @warn "No POD modes available for plotting" fom_run_id avail
            else
                nshow = min(length(filtered), 6)
                ncol = 3
                nrow = ceil(Int, nshow / ncol)

                for idx in 1:nshow
                    irow = ceil(Int, idx / ncol)
                    icol = idx - (irow - 1) * ncol
                    imode = filtered[idx]

                    ax = Axis(
                        fig[irow, icol];
                        title = "Mode $imode (temperature)",
                        xlabel = "x",
                        ylabel = "y",
                        aspect = DataAspect(),
                    )

                    mode = pod_mode_from_gram(snapshots.state, spectrum, imode)
                    copy_state_vector_to_fields!(mode, u, temp)
                    coords, temp_values = temperature_coords_and_values(temp, setup)
                    colorrange = NS.get_lims(temp_values, 5.0)
                    heatmap!(ax, coords..., temp_values; colormap = :seaborn_icefire_gradient, colorrange)
                end

                save(pod_modes_path, fig)
            end
        end
    end

    let
        plot_path = joinpath(outdir, "plots", "nusselt_number.png")
        if !(skip_existing && isfile(plot_path) && !force.nusselt)
            save_fom_nusselt_plot(plot_path, snapshots, cfg, setup, u, temp)
        end
    end

    movie_path = joinpath(outdir, "movies", "fom_evolution.mp4")
    should_make_movie = make_movie || force.movie
    if should_make_movie && !(skip_existing && isfile(movie_path) && !force.movie)
        copy_state_vector_to_fields!(@view(snapshots.state[:, end]), u, temp)
        NS.apply_bc_u!(u, cfg.tsim, setup)
        NS.apply_bc_temp!(temp, cfg.tsim, setup)

        coords = getindex.(Array.(setup.xp), setup.Ip.indices)
        temp_view = view(temp, setup.Ip)
        temp_cpu = Array(temp_view)
        temp_obs = Observable(temp_cpu)

        vort = NS.scalarfield(setup)
        NS.vorticity!(vort, u, setup)
        vort_view = view(vort, setup.Ip)
        vort_cpu = Array(vort_view)
        vort_obs = Observable(vort_cpu)

        temp_cr = (min(cfg.temp_hot, cfg.temp_cold), max(cfg.temp_hot, cfg.temp_cold))
        vort_cr = NS.get_lims(vort_cpu, 5.0)

        time_label = Observable("")
        fig = Figure(size = (1000, 450))
        Label(fig[0, 1:2], time_label)

        axT = Axis(fig[1, 1]; title = "Temperature", xlabel = "x", ylabel = "y", aspect = DataAspect())
        heatmap!(axT, coords..., temp_obs; colormap = :seaborn_icefire_gradient, colorrange = temp_cr)

        axW = Axis(fig[1, 2]; title = "Vorticity", xlabel = "x", ylabel = "y", aspect = DataAspect())
        heatmap!(axW, coords..., vort_obs; colormap = :RdBu, colorrange = vort_cr)

        Makie.record(
            fig,
            movie_path,
            0:cfg.nsnapshot;
            framerate = max(1, round(Int, cfg.nsnapshot / 5)),
        ) do i
            copy_state_vector_to_fields!(@view(snapshots.state[:, i + 1]), u, temp)
            t = snapshot_time(cfg, i)
            NS.apply_bc_u!(u, t, setup)
            NS.apply_bc_temp!(temp, t, setup)
            NS.vorticity!(vort, u, setup)

            copyto!(temp_cpu, temp_view)
            copyto!(vort_cpu, vort_view)
            notify(temp_obs)
            notify(vort_obs)
            time_label[] = "t = $(round(t; digits = 3))"
        end
    end

    return nothing
end

function plot_fom_temperature_only(fom_run_id::String)
    cfg = load_fom_config_from_run(fom_run_id)
    snapshots = load_fom_snapshots(fom_run_id)
    outdir = fom_dir(fom_run_id)

    setup, _, u, temp = cpu_setup_and_fields(cfg)

    plot_path = joinpath(outdir, "plots", "final_temperature.png")
    ensure_run_dirs!(outdir)

    copy_state_vector_to_fields!(@view(snapshots.state[:, end]), u, temp)
    NS.apply_bc_temp!(temp, cfg.tsim, setup)

    fig = Figure()
    ax = Axis(fig[1, 1]; title = "Temperature at final time", xlabel = "x", ylabel = "y", aspect = DataAspect())
    coords, temp_values = temperature_coords_and_values(temp, setup)
    colorrange = NS.get_lims(temp_values, 1.0)
    heatmap!(ax, coords..., temp_values; colormap = :seaborn_icefire_gradient, colorrange)
    save(plot_path, fig)

    return plot_path
end

function plot_fom_movie_only(fom_run_id::String)
    cfg = load_fom_config_from_run(fom_run_id)
    snapshots = load_fom_snapshots(fom_run_id)
    outdir = fom_dir(fom_run_id)

    setup, _, u, temp = cpu_setup_and_fields(cfg)
    ensure_run_dirs!(outdir)

    movie_path = joinpath(outdir, "movies", "fom_evolution.mp4")

    copy_state_vector_to_fields!(@view(snapshots.state[:, end]), u, temp)
    coords, temp_values = temperature_coords_and_values(temp, setup)

    temp_obs = Observable(copy(temp_values))
    colorrange = Observable(NS.get_lims(temp_values, 1.0))
    ttl = Observable("")

    fig = Figure()
    ax = Axis(fig[1, 1]; title = ttl, xlabel = "x", ylabel = "y", aspect = DataAspect())
    heatmap!(ax, coords..., temp_obs; colormap = :seaborn_icefire_gradient, colorrange)

    Makie.record(
        fig,
        movie_path,
        0:cfg.nsnapshot;
        framerate = max(1, round(Int, cfg.nsnapshot / 5)),
    ) do i
        copy_state_vector_to_fields!(@view(snapshots.state[:, i + 1]), u, temp)
        NS.apply_bc_temp!(temp, snapshot_time(cfg, i), setup)
        _, temp_values = temperature_coords_and_values(temp, setup)
        temp_obs[] = copy(temp_values)
        colorrange[] = NS.get_lims(temp_values, 1.0)
        ttl[] = "Temperature at t = $(round(snapshot_time(cfg, i); digits = 3))"
    end

    return movie_path
end

function rom_trajectory_alignment(traj, fomcfg::FOMConfig, snapshots)
    default_count = min(size(traj.a_rom, 2), size(traj.a_exact, 2), size(snapshots.state, 2))
    fom_snapshot_indices = hasproperty(traj, :fom_snapshot_indices) ?
        Int.(traj.fom_snapshot_indices) : collect(1:default_count)
    times = hasproperty(traj, :times) ?
        Float64.(traj.times) : collect(range(0.0, fomcfg.tsim; length = default_count))

    K = min(default_count, length(fom_snapshot_indices), length(times))
    K >= 1 || error("ROM trajectory contains no aligned snapshots")
    indices = fom_snapshot_indices[1:K]
    all(i -> 1 <= i <= size(snapshots.state, 2), indices) ||
        error("ROM trajectory contains an out-of-range parent FOM snapshot index")
    return (; K, fom_snapshot_indices = indices, times = times[1:K])
end

function plot_rom_outputs(rom_run_id::String; make_movie::Bool = true)
    rommeta = load_rom_metadata(rom_run_id)
    parent_fom_id = rommeta.parent_fom_id

    fomcfg = load_fom_config_from_run(parent_fom_id)
    snapshots = load_fom_snapshots(parent_fom_id)
    traj = load_rom_trajectory(rom_run_id)
    fit = load_rom_fit(rom_run_id)
    is_switched = hasproperty(fit, :switched) && fit.switched
    local_models = if is_switched
        names = hasproperty(fit, :local_model_names) ? fit.local_model_names : fit.local_rom_ids
        [_load_local_rom_model(run_id, names[i]) for (i, run_id) in enumerate(fit.local_rom_ids)]
    else
        nothing
    end
    rom_basis = is_switched ? local_models[1].rom_basis : fit.rom_basis
    basis_state_mean = !is_switched && hasproperty(fit, :basis_state_mean) ? fit.basis_state_mean : nothing

    outdir = rom_dir(rom_run_id)
    ensure_run_dirs!(outdir)

    # The trajectory file stores reduced coordinates (`a_rom`) but older
    # postprocessing assumed full-state reconstructions (`state_rom`). Rebuild
    # ROM states on the fly to avoid storing a potentially huge dense matrix on
    # disk for large FOMs.
    xrom = Vector{eltype(rom_basis)}(undef, size(rom_basis, 1))

    function reconstruct_state_at!(k::Integer)
        if is_switched
            active_idx = Int(traj.active_model_index[k])
            model = local_models[active_idx]
            r = model.rfit
            _reconstruct_local_state!(xrom, model, @view(traj.a_rom[1:r, k]))
            return xrom
        end
        mul!(xrom, rom_basis, @view(traj.a_rom[:, k]))
        isnothing(basis_state_mean) || (xrom .+= basis_state_mean)
        return xrom
    end

    model_label(k::Integer) = if is_switched
        ", model = $(local_models[Int(traj.active_model_index[k])].name)"
    else
        ""
    end

    function rel_state_error(x_rom::AbstractVector, x_fom::AbstractVector)
        sdiff = 0.0
        sfom = 0.0
        @inbounds for i in eachindex(x_fom, x_rom)
            xi = Float64(x_fom[i])
            di = Float64(x_rom[i]) - xi
            sdiff += di * di
            sfom += xi * xi
        end
        return sqrt(sdiff) / max(sqrt(sfom), eps(Float64))
    end

    alignment = rom_trajectory_alignment(traj, fomcfg, snapshots)
    K = alignment.K
    fom_snapshot_indices = alignment.fom_snapshot_indices
    trajectory_times = alignment.times

    training_horizon = if is_switched &&
            hasproperty(rommeta, :local_training_end_time) &&
            !isnothing(rommeta.local_training_end_time)
        h = Float64(rommeta.local_training_end_time)
        first(trajectory_times) < h < last(trajectory_times) ? h : nothing
    else
        nothing
    end

    function mark_training_horizon!(ax)
        isnothing(training_horizon) && return nothing
        vlines!(
            ax,
            [training_horizon];
            color = :red,
            linestyle = :dashdot,
            linewidth = 2,
            label = "training horizon",
        )
        return nothing
    end

    let
        errs = Vector{Float64}(undef, K)
        for k in 1:K
            xk = reconstruct_state_at!(k)
            errs[k] = rel_state_error(xk, @view(snapshots.state[:, fom_snapshot_indices[k]]))
        end

        fig = Figure()
        ax = Axis(fig[1, 1]; title = "Relative state error", xlabel = "t", ylabel = "norm(x_ROM - x_FOM) / norm(x_FOM)", yscale = log10)
        lines!(ax, trajectory_times, errs)
        mark_training_horizon!(ax)
        isnothing(training_horizon) || axislegend(ax; position = :rb)
        save(joinpath(outdir, "plots", "relative_state_error.png"), fig)
    end

    let
        nmodes = min(4, size(traj.a_exact, 1))
        fig = Figure(size = (1000, 700))

        for i in 1:nmodes
            ax = Axis(fig[i, 1]; title = "Mode coefficient a_$i", xlabel = "t", ylabel = "amplitude")
            lines!(ax, trajectory_times, @view(traj.a_exact[i, 1:K]), label = "Projected FOM")
            lines!(ax, trajectory_times, @view(traj.a_rom[i, 1:K]), label = "ROM")
            if is_switched && hasproperty(traj, :switch_events)
                for ev in traj.switch_events
                    vlines!(ax, [ev.time]; color = :gray, linestyle = :dash)
                end
            end
            mark_training_horizon!(ax)
            axislegend(ax)
        end

        save(joinpath(outdir, "plots", "modal_coefficients.png"), fig)
    end

    let
        setup = getsetup(fomcfg; backend = KernelAbstractions.CPU())
        psolver = NS.default_psolver(setup)
        state = initial_state(fomcfg, setup, psolver)
        u_fom = state.u
        temp_fom = state.temp
        u_rom = similar(u_fom)
        temp_rom = similar(temp_fom)

        ke_fom = zeros(Float64, K)
        ke_rom = zeros(Float64, K)
        nu_bottom_fom = zeros(Float64, K)
        nu_bottom_rom = zeros(Float64, K)
        nu_top_fom = zeros(Float64, K)
        nu_top_rom = zeros(Float64, K)
        lz_fom = zeros(Float64, K)
        lz_rom = zeros(Float64, K)
        angular_cache = angular_momentum_z_cache(setup, u_fom, fomcfg)

        for k in 1:K
            t = trajectory_times[k]
            copy_state_vector_to_fields!(@view(snapshots.state[:, fom_snapshot_indices[k]]), u_fom, temp_fom)
            xk = reconstruct_state_at!(k)
            copy_state_vector_to_fields!(xk, u_rom, temp_rom)

            NS.apply_bc_u!(u_fom, t, setup)
            NS.apply_bc_u!(u_rom, t, setup)
            NS.apply_bc_temp!(temp_fom, t, setup)
            NS.apply_bc_temp!(temp_rom, t, setup)

            ke_fom[k] = NS.total_kinetic_energy(u_fom, setup)
            ke_rom[k] = NS.total_kinetic_energy(u_rom, setup)

            nu_fom = plate_nusselt_numbers(temp_fom, setup, fomcfg)
            nu_rom = plate_nusselt_numbers(temp_rom, setup, fomcfg)
            nu_bottom_fom[k] = nu_fom.bottom
            nu_bottom_rom[k] = nu_rom.bottom
            nu_top_fom[k] = nu_fom.top
            nu_top_rom[k] = nu_rom.top
            lz_fom[k] = normalized_angular_momentum_z(u_fom, setup, angular_cache)
            lz_rom[k] = normalized_angular_momentum_z(u_rom, setup, angular_cache)
        end

        let
            fig = Figure(size = (1000, 600))
            ax = Axis(fig[1, 1]; title = "Total kinetic energy", xlabel = "t", ylabel = "kinetic energy")
            lines!(ax, trajectory_times, ke_fom; label = "FOM")
            lines!(ax, trajectory_times, ke_rom; label = "ROM")
            mark_training_horizon!(ax)
            axislegend(ax; position = :rb)
            save(joinpath(outdir, "plots", "kinetic_energy.png"), fig)
        end

        let
            fig = Figure(size = (1000, 600))
            ax = Axis(fig[1, 1]; title = "Nusselt number (bottom plate)", xlabel = "t", ylabel = "Nu")
            lines!(ax, trajectory_times, nu_bottom_fom; label = "FOM")
            lines!(ax, trajectory_times, nu_bottom_rom; label = "ROM")
            mark_training_horizon!(ax)
            axislegend(ax; position = :rb)
            save(joinpath(outdir, "plots", "nusselt_bottom.png"), fig)
        end

        let
            fig = Figure(size = (1000, 600))
            ax = Axis(fig[1, 1]; title = "Nusselt number (top plate)", xlabel = "t", ylabel = "Nu")
            lines!(ax, trajectory_times, nu_top_fom; label = "FOM")
            lines!(ax, trajectory_times, nu_top_rom; label = "ROM")
            mark_training_horizon!(ax)
            axislegend(ax; position = :rb)
            save(joinpath(outdir, "plots", "nusselt_top.png"), fig)
        end

        let
            fig = Figure(size = (1000, 600))
            ax = Axis(fig[1, 1]; title = "Normalized angular momentum", xlabel = "t", ylabel = "Lz normalized")
            lines!(ax, trajectory_times, lz_fom; label = "FOM")
            lines!(ax, trajectory_times, lz_rom; label = "ROM")
            hlines!(ax, [0.0]; color = :black, linestyle = :dash)
            if is_switched && hasproperty(traj, :switch_events)
                for ev in traj.switch_events
                    vlines!(ax, [ev.time]; color = :gray, linestyle = :dash)
                end
            end
            mark_training_horizon!(ax)
            axislegend(ax; position = :rb)
            save(joinpath(outdir, "plots", "angular_momentum_z_normalized.png"), fig)
        end

        if is_switched
            fig = Figure(size = (1000, 600))
            ax = Axis(fig[1, 1]; title = "Active local ROM", xlabel = "t", ylabel = "model index")
            lines!(ax, trajectory_times, Float64.(traj.active_model_index[1:K]))
            mark_training_horizon!(ax)
            isnothing(training_horizon) || axislegend(ax; position = :rb)
            save(joinpath(outdir, "plots", "active_model.png"), fig)
        end
    end

    if make_movie
        setup = getsetup(fomcfg; backend = KernelAbstractions.CPU())
        psolver = NS.default_psolver(setup)
        state = initial_state(fomcfg, setup, psolver)
        u_fom = state.u
        temp_fom = state.temp
        u_rom = similar(u_fom)
        temp_rom = similar(temp_fom)

        nframe = K - 1

        copy_state_vector_to_fields!(@view(snapshots.state[:, fom_snapshot_indices[K]]), u_fom, temp_fom)
        xK = reconstruct_state_at!(K)
        copy_state_vector_to_fields!(xK, u_rom, temp_rom)

        coords, temp_fom_vals = temperature_coords_and_values(temp_fom, setup)
        _, temp_rom_vals = temperature_coords_and_values(temp_rom, setup)

        temp_fom_obs = Observable(copy(temp_fom_vals))
        temp_rom_obs = Observable(copy(temp_rom_vals))
        suptitle = Observable("")
        colorrange = (min(fomcfg.temp_hot, fomcfg.temp_cold), max(fomcfg.temp_hot, fomcfg.temp_cold))

        movie_temp_path = joinpath(outdir, "movies", "fom_vs_rom.mp4")
        isfile(movie_temp_path) && rm(movie_temp_path; force = true)

        fig = Figure(size = (1200, 520))
        ax1 = Axis(fig[1, 1]; title = "FOM temperature", xlabel = "x", ylabel = "y", aspect = DataAspect())
        ax2 = Axis(fig[1, 2]; title = "ROM temperature", xlabel = "x", ylabel = "y", aspect = DataAspect())
        Label(fig[0, :], suptitle)

        heatmap!(ax1, coords..., temp_fom_obs; colormap = :seaborn_icefire_gradient, colorrange)
        heatmap!(ax2, coords..., temp_rom_obs; colormap = :seaborn_icefire_gradient, colorrange)

        Makie.record(
            fig,
            movie_temp_path,
            0:nframe;
            framerate = max(1, round(Int, nframe / 5)),
        ) do i
            k = i + 1
            t = trajectory_times[k]
            copy_state_vector_to_fields!(@view(snapshots.state[:, fom_snapshot_indices[k]]), u_fom, temp_fom)
            xk = reconstruct_state_at!(k)
            copy_state_vector_to_fields!(xk, u_rom, temp_rom)
            NS.apply_bc_temp!(temp_fom, t, setup)
            NS.apply_bc_temp!(temp_rom, t, setup)

            _, temp_fom_vals = temperature_coords_and_values(temp_fom, setup)
            _, temp_rom_vals = temperature_coords_and_values(temp_rom, setup)

            temp_fom_obs[] = copy(temp_fom_vals)
            temp_rom_obs[] = copy(temp_rom_vals)
            suptitle[] = "FOM vs ROM, t = $(round(t; digits = 3)), r = $(rommeta.rfit)$(model_label(k))"
        end

        vort = NS.scalarfield(setup)
        vort_view = view(vort, setup.Ip)

        copy_state_vector_to_fields!(@view(snapshots.state[:, fom_snapshot_indices[K]]), u_fom, temp_fom)
        xK = reconstruct_state_at!(K)
        copy_state_vector_to_fields!(xK, u_rom, temp_rom)
        NS.apply_bc_u!(u_fom, trajectory_times[K], setup)
        NS.apply_bc_u!(u_rom, trajectory_times[K], setup)

        NS.vorticity!(vort, u_fom, setup)
        vort_fom_cpu = Array(vort_view)
        NS.vorticity!(vort, u_rom, setup)
        vort_rom_cpu = Array(vort_view)

        vort_fom_obs = Observable(vort_fom_cpu)
        vort_rom_obs = Observable(vort_rom_cpu)
        suptitle_vort = Observable("")

        lims_fom = NS.get_lims(vort_fom_cpu, 5.0)
        lims_rom = NS.get_lims(vort_rom_cpu, 5.0)
        vort_cr = (min(lims_fom[1], lims_rom[1]), max(lims_fom[2], lims_rom[2]))

        fig = Figure(size = (1200, 520))
        ax1 = Axis(fig[1, 1]; title = "FOM vorticity", xlabel = "x", ylabel = "y", aspect = DataAspect())
        ax2 = Axis(fig[1, 2]; title = "ROM vorticity", xlabel = "x", ylabel = "y", aspect = DataAspect())
        Label(fig[0, :], suptitle_vort)

        heatmap!(ax1, coords..., vort_fom_obs; colormap = :RdBu, colorrange = vort_cr)
        heatmap!(ax2, coords..., vort_rom_obs; colormap = :RdBu, colorrange = vort_cr)

        movie_vort_path = joinpath(outdir, "movies", "fom_vs_rom_vorticity.mp4")
        isfile(movie_vort_path) && rm(movie_vort_path; force = true)

        Makie.record(
            fig,
            movie_vort_path,
            0:nframe;
            framerate = max(1, round(Int, nframe / 5)),
        ) do i
            k = i + 1
            t = trajectory_times[k]
            copy_state_vector_to_fields!(@view(snapshots.state[:, fom_snapshot_indices[k]]), u_fom, temp_fom)
            xk = reconstruct_state_at!(k)
            copy_state_vector_to_fields!(xk, u_rom, temp_rom)

            NS.apply_bc_u!(u_fom, t, setup)
            NS.apply_bc_u!(u_rom, t, setup)
            NS.vorticity!(vort, u_fom, setup)
            copyto!(vort_fom_cpu, vort_view)
            NS.vorticity!(vort, u_rom, setup)
            copyto!(vort_rom_cpu, vort_view)

            notify(vort_fom_obs)
            notify(vort_rom_obs)
            suptitle_vort[] = "FOM vs ROM, t = $(round(t; digits = 3)), r = $(rommeta.rfit)$(model_label(k))"
        end
    end

    return nothing
end
