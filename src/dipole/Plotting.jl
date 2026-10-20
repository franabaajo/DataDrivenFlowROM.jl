using Adapt
using CairoMakie
using KernelAbstractions
using LinearAlgebra
import IncompressibleNavierStokes as NS

const INTERPOLATE_OMEGA_P = getproperty(NS, Symbol("interpolate_ω_p"))

function total_enstrophy(u, setup)
    omega = NS.vorticity(u, setup)
    omega_p = INTERPOLATE_OMEGA_P(omega, setup)
    enstrophy_density = @. omega_p^2 / 2
    enstrophy_density = NS.scalewithvolume(enstrophy_density, setup)
    return sum(view(enstrophy_density, setup.Ip))
end

function cpu_setup_and_fields(cfg::FOMConfig)
    setup = getsetup(cfg; backend = KernelAbstractions.CPU())
    psolver = NS.default_psolver(setup)
    u = initial_velocity(cfg, setup, psolver)
    return setup, psolver, u
end

function save_fom_pod_modes_plot(
    snapshots_u::AbstractMatrix,
    setup,
    utemplate,
    plot_path::AbstractString;
    modelist = [1, 2, 3, 5, 10, 20],
    pod_spectrum = nothing,
)
    spectrum = if isnothing(pod_spectrum)
        w = velocity_mass_vector(setup, utemplate)
        pod_spectrum_from_weighted_gram(snapshots_u, w)
    else
        pod_spectrum
    end
    requested_modes = Int.(collect(modelist))
    available_modes = [m for m in requested_modes if m <= length(spectrum.S)]

    nshow = length(available_modes)
    ncol = max(1, min(5, nshow))
    nrow = max(1, ceil(Int, nshow / ncol))
    fig = Figure(size = (360 * ncol, 300 * nrow))

    if nshow == 0
        @warn "No POD modes available for plotting" requested_modes
        save(plot_path, fig)
        return plot_path
    end

    coords = adapt(Array, setup.xp)
    for (idx, imode) in enumerate(available_modes)
        irow = ceil(Int, idx / ncol)
        icol = idx - (irow - 1) * ncol

        ax = Axis(
            fig[irow, icol];
            title = "Mode $imode",
            xlabel = "x",
            ylabel = "y",
            aspect = DataAspect(),
        )

        mode = reshape(pod_mode_from_gram(snapshots_u, spectrum, imode), size(utemplate))
        vort = NS.vorticity(mode, setup) |> Array
        colorrange = NS.get_lims(vort, 5.0)
        heatmap!(ax, coords..., vort; colormap = :RdBu, colorrange)
    end

    save(plot_path, fig)
    return plot_path
end

function plot_fom_selected_modes(
    fom_run_id::String;
    modelist = 1:10,
    plot_name::AbstractString = "pod_modes_10.png",
)
    cfg = load_fom_config_from_run(fom_run_id)
    snapshots = load_fom_snapshots(fom_run_id)
    outdir = fom_dir(fom_run_id)
    ensure_run_dirs!(outdir)

    setup, _, u = cpu_setup_and_fields(cfg)
    plot_path = joinpath(outdir, "plots", plot_name)
    return save_fom_pod_modes_plot(snapshots.u, setup, u, plot_path; modelist = modelist)
end

function plot_fom_outputs(
    fom_run_id::String;
    make_movie::Bool = true,
    modelist = [1, 2, 3, 5, 10, 20],
    skip_existing::Bool = true,
    force = (
        vorticity = false,
        vorticity_zoom = false,
        pressure_zoom = false,
        singular_values = false,
        pod_modes = false,
        movie = false,
    ),
)
    cfg = load_fom_config_from_run(fom_run_id)
    snapshots = load_fom_snapshots(fom_run_id)
    outdir = fom_dir(fom_run_id)

    setup, psolver, u = cpu_setup_and_fields(cfg)

    singular_values_path = joinpath(outdir, "plots", "singular_values.pdf")
    pod_modes_path = joinpath(outdir, "plots", "pod_modes.png")
    needs_singular_values = !(skip_existing && isfile(singular_values_path) && !force.singular_values)
    needs_pod_modes = !(skip_existing && isfile(pod_modes_path) && !force.pod_modes)
    needs_pod = needs_singular_values || needs_pod_modes
    pod_spectrum = if needs_pod
        w = velocity_mass_vector(setup, u)
        pod_spectrum_from_weighted_gram(snapshots.u, w)
    else
        nothing
    end

    let
        plot_path = joinpath(outdir, "plots", "final_vorticity.png")
        if !(skip_existing && isfile(plot_path) && !force.vorticity)
            copyto!(u, reshape(@view(snapshots.u[:, end]), size(u)))
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
        plot_path = joinpath(outdir, "plots", "final_vorticity_zoom.png")
        if !(skip_existing && isfile(plot_path) && !force.vorticity_zoom)
            copyto!(u, reshape(@view(snapshots.u[:, end]), size(u)))
            NS.apply_bc_u!(u, cfg.tsim, setup)
            fig = Figure()
            ax = Axis(fig[1, 1]; title = "Final vorticity zoom", xlabel = "x", ylabel = "y", aspect = DataAspect())
            vort = NS.vorticity(u, setup) |> Array
            coords = getindex.(Array.(setup.xp), setup.Ip.indices)
            vort = vort[setup.Ip]
            contour!(ax, coords..., vort; levels = 80, color = :grey)
            xlims!(ax, 0.4, 1.0)
            ylims!(ax, 0.0, 0.6)
            save(plot_path, fig)
        end
    end

    let
        plot_path = joinpath(outdir, "plots", "final_pressure_zoom.png")
        if !(skip_existing && isfile(plot_path) && !force.pressure_zoom)
            copyto!(u, reshape(@view(snapshots.u[:, end]), size(u)))
            NS.apply_bc_u!(u, cfg.tsim, setup)

            fig = Figure()
            ax = Axis(
                fig[1, 1];
                title = "Final pressure zoom",
                xlabel = "x",
                ylabel = "y",
                aspect = DataAspect(),
            )

            # Compute pressure from final velocity snapshot (not stored on disk)
            # Use ∇²p = ∇·(convection + diffusion) for a meaningful field.
            f = NS.navierstokes((; u), cfg.tsim; setup, viscosity = cfg.viscosity)
            div = NS.divergence(f.u, setup)
            NS.scalewithvolume!(div, setup)
            p = NS.poisson(psolver, div)
            NS.apply_bc_p!(p, cfg.tsim, setup)

            coords = getindex.(Array.(setup.xp), setup.Ip.indices)
            p = Array(p)[setup.Ip]

            # pressure is only defined up to a constant -> shift it
            ix0 = argmin(abs.(coords[1]))
            iy0 = argmin(abs.(coords[2]))
            p .-= p[ix0, iy0]

            pmax = maximum(abs, p)
            println("pmax = ", pmax)
            levels = collect(floor(-pmax):1:ceil(pmax))

            contour!(ax, coords..., p; levels = levels, color = :grey)

            xlims!(ax, 0.4, 1.0)
            ylims!(ax, 0.0, 0.6)

            save(plot_path, fig)
        end
    end

    let
        plot_path = singular_values_path
        if needs_singular_values
            fig = Figure()
            ax = Axis(fig[1, 1]; title = "Singular values", xlabel = "Mode number", ylabel = "σᵢ / σ₁", yscale = log10)
            if isnothing(pod_spectrum) || isempty(pod_spectrum.S)
                @warn "No POD singular values available for plotting" fom_run_id
            else
                denom = max(pod_spectrum.S[1], eps(eltype(pod_spectrum.S)))
                scatter!(ax, pod_spectrum.S / denom)
            end
            save(plot_path, fig)
        end
    end

    let
        plot_path = pod_modes_path
        if needs_pod_modes
            save_fom_pod_modes_plot(
                snapshots.u,
                setup,
                u,
                plot_path;
                modelist = modelist,
                pod_spectrum = pod_spectrum,
            )
        end
    end

    movie_path = joinpath(outdir, "movies", "fom_evolution.mp4")
    should_make_movie = make_movie || force.movie
    if should_make_movie && !(skip_existing && isfile(movie_path) && !force.movie)
        vort = NS.scalarfield(setup)
        vort_cpu = Array(vort)
        vort_obs = Observable(vort_cpu)
        coords = adapt(Array, setup.xp)
        copyto!(u, reshape(@view(snapshots.u[:, end]), size(u)))
        NS.vorticity!(vort, u, setup)
        copyto!(vort_cpu, vort)
        colorrange = Observable(NS.get_lims(vort, 5.0))
        ttl = Observable("")

        fig = Figure()
        ax = Axis(fig[1, 1]; title = ttl, xlabel = "x", ylabel = "y", aspect = DataAspect())
        heatmap!(ax, coords..., vort_obs; colormap = :RdBu, colorrange)

        Makie.record(
            fig,
            movie_path,
            0:cfg.nsnapshot;
            framerate = max(1, round(Int, cfg.nsnapshot / 5)),
        ) do i
            copyto!(u, reshape(@view(snapshots.u[:, i + 1]), size(u)))
            NS.vorticity!(vort, u, setup)
            copyto!(vort_cpu, vort)
            notify(vort_obs)
            ttl[] = "Vorticity at t = $(round(cfg.tsim * i / cfg.nsnapshot; digits = 3))"
        end
    end

    return nothing
end

function plot_fom_vorticity_zoom_only(fom_run_id::String)
    cfg = load_fom_config_from_run(fom_run_id)
    snapshots = load_fom_snapshots(fom_run_id)
    outdir = fom_dir(fom_run_id)

    setup, psolver, u = cpu_setup_and_fields(cfg)

    plot_path = joinpath(outdir, "plots", "final_vorticity_zoom.png")
    ensure_run_dirs!(outdir)

    copyto!(u, reshape(@view(snapshots.u[:, end]), size(u)))
    NS.apply_bc_u!(u, cfg.tsim, setup)
    fig = Figure()
    ax = Axis(fig[1, 1]; title = "Final vorticity zoom", xlabel = "x", ylabel = "y", aspect = DataAspect())
    vort = NS.vorticity(u, setup) |> Array
    coords = getindex.(Array.(setup.xp), setup.Ip.indices)
    vort = vort[setup.Ip]
    contour!(ax, coords..., vort; levels = 80, color = :grey)
    xlims!(ax, 0.4, 1.0)
    ylims!(ax, 0.0, 0.6)
    save(plot_path, fig)

    return plot_path
end

function plot_fom_movie_only(fom_run_id::String)
    cfg = load_fom_config_from_run(fom_run_id)
    snapshots = load_fom_snapshots(fom_run_id)
    outdir = fom_dir(fom_run_id)

    setup, _, u = cpu_setup_and_fields(cfg)
    ensure_run_dirs!(outdir)

    movie_path = joinpath(outdir, "movies", "fom_evolution.mp4")

    vort = NS.scalarfield(setup)
    vort_cpu = Array(vort)
    vort_obs = Observable(vort_cpu)
    coords = adapt(Array, setup.xp)
    copyto!(u, reshape(@view(snapshots.u[:, end]), size(u)))
    NS.vorticity!(vort, u, setup)
    copyto!(vort_cpu, vort)
    colorrange = Observable(NS.get_lims(vort, 5.0))
    ttl = Observable("")

    fig = Figure()
    ax = Axis(fig[1, 1]; title = ttl, xlabel = "x", ylabel = "y", aspect = DataAspect())
    heatmap!(ax, coords..., vort_obs; colormap = :RdBu, colorrange)

    Makie.record(
        fig,
        movie_path,
        0:cfg.nsnapshot;
        framerate = max(1, round(Int, cfg.nsnapshot / 5)),
    ) do i
        copyto!(u, reshape(@view(snapshots.u[:, i + 1]), size(u)))
        NS.vorticity!(vort, u, setup)
        copyto!(vort_cpu, vort)
        notify(vort_obs)
        ttl[] = "Vorticity at t = $(round(cfg.tsim * i / cfg.nsnapshot; digits = 3))"
    end

    return movie_path
end

function plot_rom_outputs(rom_run_id::String; make_movie::Bool = true)
    romcfg = load_rom_config_from_run(rom_run_id)
    rommeta = load_rom_metadata(rom_run_id)
    parent_fom_id = rommeta.parent_fom_id

    fomcfg = load_fom_config_from_run(parent_fom_id)
    snapshots = load_fom_snapshots(parent_fom_id)
    fit = load_rom_fit(rom_run_id)
    traj = load_rom_trajectory(rom_run_id)

    outdir = rom_dir(rom_run_id)

    let
        K = min(size(traj.a_rom, 2), size(snapshots.u, 2))
        errs = zeros(Float64, K)

        if hasproperty(traj, :u_rom)
            for k in 1:K
                errs[k] =
                    norm(traj.u_rom[:, k] - snapshots.u[:, k]) / max(norm(snapshots.u[:, k]), eps(Float64))
            end
        else
            rom_basis = fit.rom_basis
            u_rom_vec = zeros(eltype(rom_basis), size(rom_basis, 1))

            for k in 1:K
                mul!(u_rom_vec, rom_basis, @view(traj.a_rom[:, k]))
                u_fom_vec = @view(snapshots.u[:, k])

                num = 0.0
                den = 0.0
                @inbounds for i in eachindex(u_rom_vec, u_fom_vec)
                    d = Float64(u_rom_vec[i]) - Float64(u_fom_vec[i])
                    f = Float64(u_fom_vec[i])
                    num += d * d
                    den += f * f
                end

                errs[k] = sqrt(num) / max(sqrt(den), eps(Float64))
            end
        end
        times = range(0.0, fomcfg.tsim; length = length(errs))

        fig = Figure()
        ax = Axis(fig[1, 1]; title = "Relative state error", xlabel = "t", ylabel = "‖u_ROM - u_FOM‖ / ‖u_FOM‖", yscale = log10)
        lines!(ax, times, errs)
        save(joinpath(outdir, "plots", "relative_state_error.png"), fig)
    end

    let
        nmodes = min(4, size(traj.a_exact, 1))
        fig = Figure(size = (1000, 700))
        tt_exact = range(0.0, fomcfg.tsim; length = size(traj.a_exact, 2))
        tt_rom = range(0.0, fomcfg.tsim; length = size(traj.a_rom, 2))

        for i in 1:nmodes
            ax = Axis(fig[i, 1]; title = "Mode coefficient a_$i", xlabel = "t", ylabel = "amplitude")
            lines!(ax, tt_exact, traj.a_exact[i, :], label = "Projected FOM")
            lines!(ax, tt_rom, traj.a_rom[i, :], label = "ROM")
            axislegend(ax)
        end

        save(joinpath(outdir, "plots", "modal_coefficients.png"), fig)
    end

    # ROM vs FOM diagnostics (kinetic energy and enstrophy)
    let
        setup = getsetup(fomcfg; backend = KernelAbstractions.CPU())
        psolver = NS.default_psolver(setup)
        u_fom = initial_velocity(fomcfg, setup, psolver)
        u_rom = initial_velocity(fomcfg, setup, psolver)
        u_fom_flat = vec(u_fom)
        u_rom_flat = vec(u_rom)

        rom_basis = fit.rom_basis
        u_rom_vec = zeros(eltype(rom_basis), size(rom_basis, 1))

        K = min(size(traj.a_rom, 2), size(snapshots.u, 2))
        times = range(0.0, fomcfg.tsim; length = K)
        ke_fom = zeros(Float64, K)
        ke_rom = zeros(Float64, K)
        en_fom = zeros(Float64, K)
        en_rom = zeros(Float64, K)

        for k in 1:K
            copyto!(u_fom_flat, @view(snapshots.u[:, k]))
            if hasproperty(traj, :u_rom)
                copyto!(u_rom_flat, @view(traj.u_rom[:, k]))
            else
                mul!(u_rom_vec, rom_basis, @view(traj.a_rom[:, k]))
                copyto!(u_rom_flat, u_rom_vec)
            end
            NS.apply_bc_u!(u_fom, times[k], setup)
            NS.apply_bc_u!(u_rom, times[k], setup)
            ke_fom[k] = NS.total_kinetic_energy(u_fom, setup)
            ke_rom[k] = NS.total_kinetic_energy(u_rom, setup)
            en_fom[k] = total_enstrophy(u_fom, setup)
            en_rom[k] = total_enstrophy(u_rom, setup)
        end

        let
            fig = Figure(size = (1000, 600))
            ax = Axis(fig[1, 1]; title = "Total kinetic energy", xlabel = "t", ylabel = "kinetic energy")
            lines!(ax, times, ke_fom; label = "FOM")
            lines!(ax, times, ke_rom; label = "ROM")
            axislegend(ax; position = :rb)
            save(joinpath(outdir, "plots", "kinetic_energy.png"), fig)
        end

        let
            fig = Figure(size = (1000, 600))
            ax = Axis(fig[1, 1]; title = "Total enstrophy", xlabel = "t", ylabel = "enstrophy")
            lines!(ax, times, en_fom; label = "FOM")
            lines!(ax, times, en_rom; label = "ROM")
            axislegend(ax; position = :rb)
            save(joinpath(outdir, "plots", "enstrophy.png"), fig)
        end
    end

    if make_movie
        setup = getsetup(fomcfg; backend = KernelAbstractions.CPU())
        psolver = NS.default_psolver(setup)
        u_fom = initial_velocity(fomcfg, setup, psolver)
        u_rom = initial_velocity(fomcfg, setup, psolver)
        u_fom_flat = vec(u_fom)
        u_rom_flat = vec(u_rom)

        vort_fom = NS.scalarfield(setup)
        vort_rom = NS.scalarfield(setup)
        vort_fom_cpu = Array(vort_fom)
        vort_rom_cpu = Array(vort_rom)

        vort_fom_obs = Observable(vort_fom_cpu)
        vort_rom_obs = Observable(vort_rom_cpu)
        suptitle = Observable("")

        coords = adapt(Array, setup.xp)
        K = min(size(traj.a_rom, 2), size(snapshots.u, 2))
        rom_basis = fit.rom_basis
        u_rom_vec = zeros(eltype(rom_basis), size(rom_basis, 1))

        copyto!(u_fom_flat, @view(snapshots.u[:, K]))
        if hasproperty(traj, :u_rom)
            copyto!(u_rom_flat, @view(traj.u_rom[:, K]))
        else
            mul!(u_rom_vec, rom_basis, @view(traj.a_rom[:, K]))
            copyto!(u_rom_flat, u_rom_vec)
        end
        NS.vorticity!(vort_fom, u_fom, setup)
        NS.vorticity!(vort_rom, u_rom, setup)
        copyto!(vort_fom_cpu, vort_fom)
        copyto!(vort_rom_cpu, vort_rom)
        lims_fom = NS.get_lims(vort_fom, 5.0)
        lims_rom = NS.get_lims(vort_rom, 5.0)
        colorrange = Observable((min(lims_fom[1], lims_rom[1]), max(lims_fom[2], lims_rom[2])))

        fig = Figure(size = (1200, 520))
        ax1 = Axis(fig[1, 1]; title = "FOM", xlabel = "x", ylabel = "y", aspect = DataAspect())
        ax2 = Axis(fig[1, 2]; title = "ROM", xlabel = "x", ylabel = "y", aspect = DataAspect())
        Label(fig[0, :], suptitle)

        heatmap!(ax1, coords..., vort_fom_obs; colormap = :RdBu, colorrange)
        heatmap!(ax2, coords..., vort_rom_obs; colormap = :RdBu, colorrange)

        nframe = K - 1
        Makie.record(
            fig,
            joinpath(outdir, "movies", "fom_vs_rom.mp4"),
            0:nframe;
            framerate = max(1, round(Int, nframe / 5)),
        ) do i
            copyto!(u_fom_flat, @view(snapshots.u[:, i + 1]))
            if hasproperty(traj, :u_rom)
                copyto!(u_rom_flat, @view(traj.u_rom[:, i + 1]))
            else
                mul!(u_rom_vec, rom_basis, @view(traj.a_rom[:, i + 1]))
                copyto!(u_rom_flat, u_rom_vec)
            end

            NS.vorticity!(vort_fom, u_fom, setup)
            NS.vorticity!(vort_rom, u_rom, setup)

            copyto!(vort_fom_cpu, vort_fom)
            copyto!(vort_rom_cpu, vort_rom)
            notify(vort_fom_obs)
            notify(vort_rom_obs)

            suptitle[] = "FOM vs ROM, t = $(round(fomcfg.tsim * i / nframe; digits = 3)), r = $(rommeta.rfit)"
        end
    end

    return nothing
end
