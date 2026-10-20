# # Dipole FOM enstrophy comparison
#
# Compare enstrophy histories for multiple FOM runs, and compute absolute enstrophy
# errors relative to a reference run.
#
# Usage:
#   julia --project=. experiments/DipoleWallCollision/studies/DipoleFOMEnstrophyCompare.jl \
#       --ref <reference_fom_run_id> <fom_run_id_1> [fom_run_id_2 ...]
#
# Outputs (in `output/DipoleWallCollision/dipolepod/diagnostics` by default):
#   - enstrophy.png
#   - enstrophy_abs_error.png

include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Config.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Paths.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Setup.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "FOM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Stats.jl"))

using CairoMakie

set_theme!(
    Theme(
        fontsize = 32,
        Axis = (
            xlabelsize = 34,
            ylabelsize = 34,
            xticklabelsize = 30,
            yticklabelsize = 30,
        ),
        Legend = (
            labelsize = 27,
            patchsize = (48, 22),
        ),
    ),
)

const THESIS_LINEWIDTH = 5.0
const THESIS_REFERENCE_LINEWIDTH = 5.6

function usage_and_exit()
    error(
        "Usage:\n" *
        "  julia --project=. experiments/DipoleWallCollision/studies/DipoleFOMEnstrophyCompare.jl " *
        "--ref <reference_fom_run_id> <fom_run_id_1> [fom_run_id_2 ...]\n",
    )
end

function parse_args(args::Vector{String})
    ref_run_id = nothing
    run_ids = String[]

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help", "help")
            usage_and_exit()
        elseif arg == "--ref"
            if i == length(args)
                usage_and_exit()
            end
            ref_run_id = args[i + 1]
            i += 2
            continue
        elseif startswith(arg, "--")
            error("Unknown flag: $arg")
        else
            push!(run_ids, arg)
        end
        i += 1
    end

    if isnothing(ref_run_id) || isempty(run_ids)
        usage_and_exit()
    end

    return ref_run_id::String, run_ids
end

function unique_preserve_order(ids::AbstractVector{String})
    seen = Set{String}()
    out = String[]
    for id in ids
        if id in seen
            continue
        end
        push!(out, id)
        push!(seen, id)
    end
    return out
end

function linear_interpolate(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, xq::AbstractVector{<:Real})
    @assert length(x) == length(y)
    @assert issorted(x)
    yq = Vector{Float64}(undef, length(xq))

    x1 = x[1]
    xend = x[end]
    for (i, t) in pairs(xq)
        if t <= x1
            yq[i] = y[1]
            continue
        elseif t >= xend
            yq[i] = y[end]
            continue
        end

        j = searchsortedlast(x, t)
        @inbounds begin
            x0 = x[j]
            x2 = x[j + 1]
            y0 = y[j]
            y2 = y[j + 1]
        end
        w = (t - x0) / (x2 - x0)
        yq[i] = (1 - w) * y0 + w * y2
    end

    return yq
end

function diagnostics_root_dir()
    outdir = joinpath(output_root(), "diagnostics")
    mkpath(outdir)
    return outdir
end

function compute_enstrophy_dataset(run_id::String)
    data = compute_fom_diagnostics(run_id; diagnostics = (:enstrophy,))
    return (; run_id = data.run_id, cfg = data.cfg, times = data.times, enstrophy = data.results[:enstrophy])
end

function grid_label(cfg::FOMConfig; reference = false)
    label = "$(cfg.n) x $(cfg.n)"
    return reference ? "$label reference" : label
end

function save_enstrophy_compare_plot(datasets; outpath::AbstractString)
    fig = Figure(size = (1250, 780), figure_padding = 18)
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "total enstrophy")

    colors = Makie.wong_colors()
    for (i, data) in enumerate(datasets)
        color = colors[mod1(i, length(colors))]
        label = grid_label(data.cfg; reference = i == 1)
        linewidth = i == 1 ? THESIS_REFERENCE_LINEWIDTH : THESIS_LINEWIDTH
        lines!(ax, data.times, data.enstrophy; label, color, linewidth)
    end

    xlims!(ax, 0.0, maximum(d -> d.times[end], datasets))
    axislegend(ax; position = :rt)
    save(outpath, fig)
    return fig
end

function save_enstrophy_abs_error_plot(datasets; ref_run_id::String, outpath::AbstractString)
    ref = only(filter(d -> d.run_id == ref_run_id, datasets))
    tmax_plot = maximum(d -> d.times[end], datasets)

    fig = Figure(size = (1250, 780), figure_padding = 18)
    ax = Axis(
        fig[1, 1];
        xlabel = "time",
        ylabel = "absolute enstrophy error",
    )

    colors = Makie.wong_colors()
    plotted = 0
    for (i, data) in enumerate(datasets)
        if data.run_id == ref_run_id
            continue
        end

        t = data.times
        en = data.enstrophy
        ref_en = linear_interpolate(ref.times, ref.enstrophy, t)
        err = abs.(en .- ref_en)

        color = colors[mod1(i, length(colors))]
        label = grid_label(data.cfg)
        lines!(ax, t, err; label, color, linewidth = THESIS_LINEWIDTH)
        plotted += 1
    end

    if plotted == 0
        error("No non-reference runs provided; nothing to plot for absolute error.")
    end

    xlims!(ax, 0.0, tmax_plot)
    axislegend(ax; position = :rt)
    save(outpath, fig)
    return fig
end

function main()
    ref_run_id, run_ids_in = parse_args(copy(ARGS))
    run_ids = unique_preserve_order(vcat([ref_run_id], run_ids_in))

    @info "Computing enstrophy diagnostics for FOM runs" ref_run_id run_ids
    datasets = [compute_enstrophy_dataset(run_id) for run_id in run_ids]

    outdir = diagnostics_root_dir()
    compare_path = joinpath(outdir, "enstrophy.png")
    err_path = joinpath(outdir, "enstrophy_abs_error.png")

    save_enstrophy_compare_plot(datasets; outpath = compare_path)
    save_enstrophy_abs_error_plot(datasets; ref_run_id, outpath = err_path)

    println("Saved: ", compare_path)
    println("Saved: ", err_path)
end

main()
