using CairoMakie
using KernelAbstractions
using Printf
import IncompressibleNavierStokes as NS

set_theme!(
    Theme(
        fontsize = 24,
        Axis = (
            xlabelsize = 26,
            ylabelsize = 26,
            xticklabelsize = 22,
            yticklabelsize = 22,
            titlesize = 24,
        ),
        Colorbar = (
            labelsize = 24,
            ticklabelsize = 20,
        ),
    ),
)

include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Config.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Paths.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Registry.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Setup.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "FOM.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "POD.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Stats.jl"))

const THESIS_FIGURE_DIR = joinpath(@__DIR__, "..", "results", "figures")

const REVERSAL_FOM_ID =
    "fom_nx256_ny256_cfl0p67_stxuniform_styuniform_nu0p00025_kappa0p00025_c5e2b71f"

const SNAPSHOT_TARGETS = [
    (time = 15.0, phase = "transient"),
    (time = 60.0, phase = "pre-reversal"),
    (time = 92.1, phase = "reversal"),
    (time = 120.0, phase = "post-reversal"),
]

const DOMAIN_TICKS = [0.0, 0.5, 1.0]

figure_dir() = mkpath(THESIS_FIGURE_DIR)

function nearest_snapshot_index(times, target)
    return argmin(abs.(times .- target))
end

function downsample_indices(n; target = 15)
    step = max(1, round(Int, n / target))
    return 1:step:n
end

function snapshot_fields!(state, snapshots, setup, cfg, isnap, t)
    copy_state_vector_to_fields!(@view(snapshots.state[:, isnap]), state.u, state.temp)
    NS.apply_bc_u!(state.u, t, setup)
    NS.apply_bc_temp!(state.temp, t, setup)

    x, y = getindex.(Array.(setup.xp), setup.Ip.indices)
    temperature = Matrix(Array(state.temp)[setup.Ip])

    up = Array(NS.interpolate_u_p(state.u, setup))
    ux = Matrix(up[setup.Ip, 1])
    uy = Matrix(up[setup.Ip, 2])
    speed = sqrt.(ux .^ 2 .+ uy .^ 2)

    return (; x, y, temperature, ux, uy, maxspeed = maximum(speed))
end

function arrow_vectors(field, ix, iy, maxspeed; arrow_scale = 0.078)
    xs = Float64[]
    ys = Float64[]
    us = Float64[]
    vs = Float64[]
    scale = arrow_scale / max(maxspeed, eps(Float64))

    for j in iy, i in ix
        push!(xs, field.x[i])
        push!(ys, field.y[j])
        push!(us, scale * field.ux[i, j])
        push!(vs, scale * field.uy[i, j])
    end

    return xs, ys, us, vs
end

function plot_reversal_snapshots(run_id, outpath; targets = SNAPSHOT_TARGETS)
    cfg = load_fom_config_from_run(run_id)
    snapshots = load_fom_snapshots(run_id)
    setup = getsetup(cfg; backend = KernelAbstractions.CPU())
    psolver = NS.default_psolver(setup)
    state = initial_state(cfg, setup, psolver)
    angular_cache = angular_momentum_z_cache(setup, state.u, cfg)

    stored_times = collect(range(0.0, cfg.tsim; length = cfg.nsnapshot + 1))
    selected = []
    for target in targets
        isnap = nearest_snapshot_index(stored_times, target.time)
        t = stored_times[isnap]
        field = snapshot_fields!(state, snapshots, setup, cfg, isnap, t)
        lz = normalized_angular_momentum_z(state.u, setup, angular_cache)
        push!(selected, (; target.phase, time = t, lz, field))
    end

    maxspeed = maximum(item.field.maxspeed for item in selected)
    ix = downsample_indices(length(selected[1].field.x))
    iy = downsample_indices(length(selected[1].field.y))

    fig = Figure(size = (2300, 560))
    hm = nothing
    for (i, item) in enumerate(selected)
        ax = Axis(
            fig[1, i];
            title = @sprintf("t = %.2f", item.time),
            xlabel = "x",
            ylabel = "y",
            xticks = DOMAIN_TICKS,
            yticks = DOMAIN_TICKS,
            aspect = DataAspect(),
        )
        hm = heatmap!(
            ax,
            item.field.x,
            item.field.y,
            item.field.temperature;
            colormap = :seaborn_icefire_gradient,
            colorrange = (min(cfg.temp_hot, cfg.temp_cold), max(cfg.temp_hot, cfg.temp_cold)),
        )
        xs, ys, us, vs = arrow_vectors(item.field, ix, iy, maxspeed)
        arrows!(
            ax,
            xs,
            ys,
            us,
            vs;
            color = (:white, 0.9),
            linewidth = 1.8,
            arrowsize = 13,
            lengthscale = 1.0,
            align = :center,
        )
        xlims!(ax, cfg.xlims...)
        ylims!(ax, cfg.ylims...)
    end
    Colorbar(fig[1, 5], hm; label = "temperature")

    save(outpath, fig)

    for item in selected
        @info "Selected Rayleigh snapshot" phase = item.phase time = item.time lz = item.lz
    end
    return outpath
end

function main()
    figs = figure_dir()
    outpath = joinpath(figs, "rayleigh_fom_reversal_snapshots.png")
    @info "Writing Rayleigh reversal snapshots" REVERSAL_FOM_ID outpath
    plot_reversal_snapshots(REVERSAL_FOM_ID, outpath)
    println("Wrote Rayleigh reversal snapshots to: ", outpath)
end

main()
