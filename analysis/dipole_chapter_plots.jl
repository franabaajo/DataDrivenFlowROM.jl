using CairoMakie
using JLD2
using Statistics

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
        Legend = (
            labelsize = 20,
            patchsize = (34, 18),
        ),
    ),
)

const THESIS_LINEWIDTH = 3.4
const THESIS_REFERENCE_LINEWIDTH = 3.8

analysis_dir() = mkpath(joinpath(@__DIR__, "..", "results", "tables", "dipole"))
figure_dir() = mkpath(joinpath(@__DIR__, "..", "results", "figures"))

function load_plot_data()
    return load_object(joinpath(analysis_dir(), "plot_data.jld2"))
end

function plot_method_comparison(data)
    figs = figure_dir()
    labels = data.series.labels
    colors = Makie.wong_colors()
    styles = Dict(
        "Exact OpInf, eta=0.95" => (color = colors[4], linestyle = :dash, linewidth = THESIS_LINEWIDTH),
        "Exact OpInf, eta=0.99" => (color = colors[1], linestyle = :solid, linewidth = THESIS_LINEWIDTH),
        "Intrusive Galerkin, eta=0.99" => (color = colors[2], linestyle = :dash, linewidth = THESIS_LINEWIDTH),
        "Trajectory OpInf, eta=0.99 request" => (color = colors[3], linestyle = :solid, linewidth = THESIS_LINEWIDTH),
    )
    pretty = Dict(
        "Exact OpInf, eta=0.95" => "Exact OpInf, eta=0.95",
        "Exact OpInf, eta=0.99" => "Exact OpInf, eta=0.99",
        "Intrusive Galerkin, eta=0.99" => "Intrusive, eta=0.99",
        "Trajectory OpInf, eta=0.99 request" => "Trajectory OpInf",
    )

    fig = Figure(size = (1150, 720))
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "relative state error", yscale = log10)
    for label in labels
        lines!(ax, data.series.times, data.series.rel_error[label]; label = pretty[label], styles[label]...)
    end
    axislegend(ax; position = :lt)
    save(joinpath(figs, "dipole_rom_relative_error.png"), fig)

    fig = Figure(size = (1150, 900))
    ax1 = Axis(fig[1, 1]; xlabel = "time", ylabel = "kinetic energy")
    ax2 = Axis(fig[2, 1]; xlabel = "time", ylabel = "enstrophy")
    lines!(ax1, data.series.times, data.series.ke_fom; label = "FOM", color = :black, linewidth = THESIS_REFERENCE_LINEWIDTH)
    lines!(ax2, data.series.times, data.series.en_fom; label = "FOM", color = :black, linewidth = THESIS_REFERENCE_LINEWIDTH)
    for label in labels
        lines!(ax1, data.series.times, data.series.ke_rom[label]; label = pretty[label], styles[label]...)
        lines!(ax2, data.series.times, data.series.en_rom[label]; label = pretty[label], styles[label]...)
    end
    axislegend(ax1; position = :rt)
    axislegend(ax2; position = :rt)
    save(joinpath(figs, "dipole_rom_diagnostics.png"), fig)
    save(joinpath(figs, "dipole_rom_energy_enstrophy.png"), fig)
end

function plot_final_vorticity(data)
    figs = figure_dir()
    labels = data.vort.labels
    fields = [data.vort.fields[label] for label in labels]
    clim = maximum(quantile(abs.(vec(field)), 0.99) for field in fields)
    clim = max(clim, eps(Float64))

    fig = Figure(size = (1450, 900), figure_padding = 18)
    hm = nothing
    for (i, label) in enumerate(labels)
        row = i <= 2 ? 1 : 2
        col = i <= 2 ? i : i - 2
        ax = Axis(fig[row, col]; title = label, xlabel = "x", ylabel = "y", aspect = DataAspect())
        plot = heatmap!(
            ax,
            data.vort.x,
            data.vort.y,
            data.vort.fields[label];
            colormap = :RdBu,
            colorrange = (-clim, clim),
        )
        hm = isnothing(hm) ? plot : hm
    end
    Colorbar(fig[:, 3], hm; label = "vorticity")
    save(joinpath(figs, "dipole_rom_final_vorticity_zoom.png"), fig)
end

function main()
    data = load_plot_data()
    plot_method_comparison(data)
    plot_final_vorticity(data)
    println("Wrote thesis figures to: ", figure_dir())
end

main()
