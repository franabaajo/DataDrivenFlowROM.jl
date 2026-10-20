using CairoMakie
using JLD2
using KernelAbstractions
import IncompressibleNavierStokes as NS

const INTERPOLATE_OMEGA_P = getproperty(NS, Symbol("interpolate_ω_p"))

function total_enstrophy(u, setup)
    omega = NS.vorticity(u, setup)
    omega_p = INTERPOLATE_OMEGA_P(omega, setup)
    enstrophy_density = @. omega_p^2 / 2
    enstrophy_density = NS.scalewithvolume(enstrophy_density, setup)
    return sum(view(enstrophy_density, setup.Ip))
end

function fom_times(cfg::FOMConfig)
    return range(0.0, cfg.tsim; length = cfg.nsnapshot + 1)
end

function compute_fom_diagnostics(run_id::String; diagnostics = (:kinetic_energy, :enstrophy))
    cfg = load_fom_config_from_run(run_id)
    snapshots = load_fom_snapshots(run_id)
    setup = getsetup(cfg; backend = KernelAbstractions.CPU())
    u = NS.vectorfield(setup)

    times = collect(fom_times(cfg))
    K = length(times)
    results = Dict{Symbol,Vector{Float64}}()

    if :kinetic_energy in diagnostics
        results[:kinetic_energy] = zeros(Float64, K)
    end
    if :enstrophy in diagnostics
        results[:enstrophy] = zeros(Float64, K)
    end

    for i in 1:K
        copyto!(u, reshape(@view(snapshots.u[:, i]), size(u)))
        NS.apply_bc_u!(u, times[i], setup)

        if haskey(results, :kinetic_energy)
            results[:kinetic_energy][i] = NS.total_kinetic_energy(u, setup)
        end
        if haskey(results, :enstrophy)
            results[:enstrophy][i] = total_enstrophy(u, setup)
        end
    end

    return (; run_id, cfg, times, results)
end

function diagnostics_output_dir(run_id::String)
    outdir = joinpath(fom_dir(run_id), "plots", "diagnostics")
    mkpath(outdir)
    return outdir
end

function save_fom_diagnostic_plots(run_id::String; diagnostics = (:kinetic_energy, :enstrophy))
    data = compute_fom_diagnostics(run_id; diagnostics)
    outdir = diagnostics_output_dir(run_id)

    for (name, values) in data.results
        fig = Figure(size = (900, 500))
        ax = Axis(fig[1, 1]; title = string(name), xlabel = "t", ylabel = string(name))
        lines!(ax, data.times, values)
        save(joinpath(outdir, "$(name).png"), fig)
    end

    save_object(joinpath(outdir, "diagnostics.jld2"), data)
    return outdir
end

function compare_output_dir()
    outdir = joinpath(output_root(), "diagnostics", "fom_compare")
    mkpath(outdir)
    return outdir
end

function run_label(cfg::FOMConfig, run_id::String)
    return "n=$(cfg.n), cfl=$(cfg.cfl), st=$(cfg.stretch), nu=$(cfg.viscosity)"
end

function save_fom_diagnostics_compare(run_ids::AbstractVector{String};
    diagnostics = (:kinetic_energy, :enstrophy),
)
    datasets = [compute_fom_diagnostics(run_id; diagnostics) for run_id in run_ids]
    outdir = compare_output_dir()
    colors = Makie.wong_colors()

    for diag in diagnostics
        fig = Figure(size = (1000, 600))
        ax = Axis(fig[1, 1]; title = string(diag), xlabel = "t", ylabel = string(diag))
        for (i, data) in enumerate(datasets)
            vals = data.results[diag]
            color = colors[mod1(i, length(colors))]
            label = run_label(data.cfg, data.run_id)
            lines!(ax, data.times, vals; label, color)
        end
        axislegend(ax; position = :rb)
        save(joinpath(outdir, "$(diag).png"), fig)
    end

    save_object(joinpath(outdir, "diagnostics_compare.jld2"), (; datasets))
    return outdir
end
