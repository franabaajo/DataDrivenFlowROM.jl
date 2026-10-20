using CairoMakie
using KernelAbstractions
using Printf
using Statistics
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

include(joinpath(@__DIR__, "..", "src", "dipole", "Config.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "Paths.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "Setup.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "FOM.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "Stats.jl"))

const THESIS_FIGURE_DIR = joinpath(@__DIR__, "..", "results", "figures")

const SUPERVISOR_ASSET_DIR = THESIS_FIGURE_DIR

const VALIDATION_RUNS = [
    (n = 768, run_id = "fom_n768_cfl0p4_st1p0_ns100_nu0p0004_19d3d5a0"),
    (n = 1024, run_id = "fom_n1024_cfl0p67_st1p2_ns100_nu0p0004_b702e129"),
    (n = 1536, run_id = "fom_n1536_cfl0p67_st1p2_ns50_nu0p0004_e592f3a8"),
]

const REFERENCE_RUN_ID = VALIDATION_RUNS[end].run_id
const ZOOM_XLIMS = (0.4, 1.0)
const ZOOM_YLIMS = (0.0, 0.6)
const ZOOM_XTICKS = [0.4, 0.6, 0.8, 1.0]
const ZOOM_YTICKS = [0.0, 0.2, 0.4, 0.6]

function figure_dir()
    mkpath(THESIS_FIGURE_DIR)
    return THESIS_FIGURE_DIR
end

function copy_required_asset(src, dst)
    if !isfile(src)
        error("Required validation asset not found: $src")
    end
    cp(src, dst; force = true)
    return dst
end

function setup_and_velocity(run_id)
    cfg = load_fom_config_from_run(run_id)
    setup = getsetup(cfg; backend = KernelAbstractions.CPU())
    psolver = NS.default_psolver(setup)
    u = initial_velocity(cfg, setup, psolver)
    return cfg, setup, psolver, u
end

function vorticity_matrix!(u, xvec, setup, t)
    copyto!(vec(u), xvec)
    NS.apply_bc_u!(u, t, setup)
    return Matrix((NS.vorticity(u, setup) |> Array)[setup.Ip])
end

function pressure_matrix!(u, xvec, cfg, setup, psolver, t)
    copyto!(vec(u), xvec)
    NS.apply_bc_u!(u, t, setup)

    # Pressure is reconstructed from the final velocity because it is not stored in snapshots.
    f = NS.navierstokes((; u), t; setup, viscosity = cfg.viscosity)
    div = NS.divergence(f.u, setup)
    NS.scalewithvolume!(div, setup)
    p = NS.poisson(psolver, div)
    NS.apply_bc_p!(p, t, setup)

    coords = getindex.(Array.(setup.xp), setup.Ip.indices)
    mat = Matrix((p |> Array)[setup.Ip])

    ix0 = argmin(abs.(coords[1]))
    iy0 = argmin(abs.(coords[2]))
    mat .-= mat[ix0, iy0]
    return mat
end

function downsample_indices(n; target = 520)
    step = max(1, ceil(Int, n / target))
    return 1:step:n
end

const DOMAIN_TICKS = [-1.0, -0.5, 0.0, 0.5, 1.0]

function plot_final_vorticity_zoom(run_id, outpath)
    cfg, setup, _psolver, u = setup_and_velocity(run_id)
    snapshots = load_fom_snapshots(run_id)
    x, y = getindex.(Array.(setup.xp), setup.Ip.indices)
    mat = vorticity_matrix!(u, @view(snapshots.u[:, end]), setup, cfg.tsim)

    fig = Figure(size = (560, 560), figure_padding = 25)
    ax = Axis(
        fig[1, 1];
        aspect = DataAspect(),
        xticks = ZOOM_XTICKS,
        yticks = ZOOM_YTICKS,
        xticklabelsize = 16,
        yticklabelsize = 16,
        xgridvisible = false,
        ygridvisible = false,
    )
    contour!(ax, x, y, mat; levels = 80, color = :grey35, linewidth = 1.25)
    xlims!(ax, ZOOM_XLIMS...)
    ylims!(ax, ZOOM_YLIMS...)
    save(outpath, fig)
    return outpath
end

function plot_final_pressure_zoom(run_id, outpath)
    cfg, setup, psolver, u = setup_and_velocity(run_id)
    snapshots = load_fom_snapshots(run_id)
    x, y = getindex.(Array.(setup.xp), setup.Ip.indices)
    mat = pressure_matrix!(u, @view(snapshots.u[:, end]), cfg, setup, psolver, cfg.tsim)

    pmax = maximum(abs, mat)
    levels = collect(floor(-pmax):1:ceil(pmax))

    fig = Figure(size = (560, 560), figure_padding = 25)
    ax = Axis(
        fig[1, 1];
        aspect = DataAspect(),
        xticks = ZOOM_XTICKS,
        yticks = ZOOM_YTICKS,
        xticklabelsize = 16,
        yticklabelsize = 16,
        xgridvisible = false,
        ygridvisible = false,
    )
    contour!(ax, x, y, mat; levels = levels, color = :grey35, linewidth = 1.25)
    xlims!(ax, ZOOM_XLIMS...)
    ylims!(ax, ZOOM_YLIMS...)
    save(outpath, fig)
    return outpath
end

function plot_reference_snapshots(run_id, outpath; target_times = [0.0, 0.33, 0.61, 1.0])
    cfg, setup, _psolver, u = setup_and_velocity(run_id)
    snapshots = load_fom_snapshots(run_id)
    x, y = getindex.(Array.(setup.xp), setup.Ip.indices)
    ix = downsample_indices(length(x))
    iy = downsample_indices(length(y))

    stored_times = collect(range(0.0, cfg.tsim; length = size(snapshots.u, 2)))
    indices = [argmin(abs.(stored_times .- t)) for t in target_times]

    fields = Matrix{Float64}[]
    actual_times = Float64[]
    for k in indices
        mat = vorticity_matrix!(u, @view(snapshots.u[:, k]), setup, stored_times[k])
        push!(fields, mat[ix, iy])
        push!(actual_times, stored_times[k])
    end

    absvals = reduce(vcat, (vec(abs.(field)) for field in fields))
    clim = quantile(absvals, 0.995)

    fig = Figure(size = (2300, 560))
    hm = nothing
    for i in eachindex(fields)
        ax = Axis(
            fig[1, i];
            title = @sprintf("t = %.2f", actual_times[i]),
            xlabel = "x",
            ylabel = "y",
            xticks = DOMAIN_TICKS,
            yticks = DOMAIN_TICKS,
            aspect = DataAspect(),
        )
        hm = heatmap!(
            ax,
            x[ix],
            y[iy],
            fields[i];
            colormap = :RdBu,
            colorrange = (-clim, clim),
        )
        xlims!(ax, cfg.lims...)
        ylims!(ax, cfg.lims...)
    end
    Colorbar(fig[1, 5], hm; label = "vorticity")
    save(outpath, fig)
    return outpath
end

function main()
    figs = figure_dir()

    copy_required_asset(
        joinpath(SUPERVISOR_ASSET_DIR, "vorticity benchmark.png"),
        joinpath(figs, "benchmark_vorticity_zoom.png"),
    )
    copy_required_asset(
        joinpath(SUPERVISOR_ASSET_DIR, "pressure benchmark.png"),
        joinpath(figs, "benchmark_pressure_zoom.png"),
    )

    for item in VALIDATION_RUNS
        outpath = joinpath(figs, "final_vorticity_zoom_$(item.n).png")
        @info "Writing validation vorticity zoom" item.run_id outpath
        plot_final_vorticity_zoom(item.run_id, outpath)

        pressure_outpath = joinpath(figs, "final_pressure_zoom_$(item.n).png")
        @info "Writing validation pressure zoom" item.run_id pressure_outpath
        plot_final_pressure_zoom(item.run_id, pressure_outpath)
    end

    snapshot_path = joinpath(figs, "dipole_fom_vorticity_snapshots.png")
    @info "Writing validation vorticity snapshots" REFERENCE_RUN_ID snapshot_path
    plot_reference_snapshots(REFERENCE_RUN_ID, snapshot_path)

    println("Wrote validation figures to: ", figs)
end

main()
