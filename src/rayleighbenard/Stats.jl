using CairoMakie
using JLD2
using KernelAbstractions
import IncompressibleNavierStokes as NS

include(joinpath(@__DIR__, "ThermalDiagnostics.jl"))

mean_temperature(temp, setup, cfg::FOMConfig) = scalar_volume_average(temp, setup, cfg)

function fom_times(cfg::FOMConfig)
    return range(0.0, cfg.tsim; length = cfg.nsnapshot + 1)
end

function angular_momentum_z_cache(setup, utemplate, cfg::FOMConfig)
    setup.dimension() == 2 || error("angular_momentum_z currently supports 2D setups")

    w = velocity_mass_vector(setup, utemplate)
    xcoords, ycoords = velocity_xycoords(setup; T = Float64)
    dims = size(utemplate)
    D = setup.dimension()
    lastdim = ndims(utemplate)
    comp_field = fill(0, dims)
    for α in 1:D
        selectdim(comp_field, lastdim, α) .= α
    end
    comp = vec(comp_field)

    xc = 0.5 * (cfg.xlims[1] + cfg.xlims[2])
    yc = 0.5 * (cfg.ylims[1] + cfg.ylims[2])
    coeff = zeros(Float64, length(w))
    inertia = 0.0

    @inbounds for k in eachindex(w)
        wk = Float64(w[k])
        wk == 0.0 && continue
        dx = xcoords[k] - xc
        dy = ycoords[k] - yc
        inertia += wk * (dx^2 + dy^2)
        if comp[k] == 1
            coeff[k] = -dy * wk
        elseif comp[k] == 2
            coeff[k] = dx * wk
        end
    end

    return (; coeff, inertia)
end

angular_momentum_z(u, cache) = dot(cache.coeff, vec(Array(u)))

function normalized_angular_momentum_z(u, setup, cache)
    lz = angular_momentum_z(u, cache)
    speed_norm2 = 2.0 * NS.total_kinetic_energy(u, setup)
    denom = sqrt(max(cache.inertia * speed_norm2, eps(Float64)))
    return lz / denom
end

function compute_fom_diagnostics(
    run_id::String;
    diagnostics = (
        :kinetic_energy,
        :angular_momentum_z,
        :angular_momentum_z_normalized,
        :nusselt_bottom,
        :nusselt_top,
        :nusselt_volume,
        :mean_temperature,
        :mean_temperature_rate,
        :temperature_source,
        :heat_budget_residual,
    ),
)
    cfg = load_fom_config_from_run(run_id)
    snapshots = load_fom_snapshots(run_id)
    setup = getsetup(cfg; backend = KernelAbstractions.CPU())
    psolver = NS.default_psolver(setup)
    state = initial_state(cfg, setup, psolver)
    u = state.u
    temp = state.temp

    times = collect(fom_times(cfg))
    K = length(times)
    results = Dict{Symbol,Vector{Float64}}()

    if :kinetic_energy in diagnostics
        results[:kinetic_energy] = zeros(Float64, K)
    end
    if :angular_momentum_z in diagnostics
        results[:angular_momentum_z] = zeros(Float64, K)
    end
    if :angular_momentum_z_normalized in diagnostics
        results[:angular_momentum_z_normalized] = zeros(Float64, K)
    end
    if :nusselt_bottom in diagnostics
        results[:nusselt_bottom] = zeros(Float64, K)
    end
    if :nusselt_top in diagnostics
        results[:nusselt_top] = zeros(Float64, K)
    end
    if :nusselt_volume in diagnostics
        results[:nusselt_volume] = zeros(Float64, K)
    end
    if :mean_temperature in diagnostics
        results[:mean_temperature] = zeros(Float64, K)
    end
    if :mean_temperature_rate in diagnostics
        results[:mean_temperature_rate] = zeros(Float64, K)
    end
    if :temperature_source in diagnostics
        results[:temperature_source] = zeros(Float64, K)
    end
    if :heat_budget_residual in diagnostics
        results[:heat_budget_residual] = zeros(Float64, K)
    end

    dtemp = similar(temp)
    needs_angular_momentum = :angular_momentum_z in diagnostics ||
        :angular_momentum_z_normalized in diagnostics
    angular_cache = needs_angular_momentum ? angular_momentum_z_cache(setup, u, cfg) : nothing

    for i in 1:K
        copy_state_vector_to_fields!(@view(snapshots.state[:, i]), u, temp)
        copyto!(dtemp, reshape(@view(snapshots.dtemp[:, i]), size(dtemp)))
        NS.apply_bc_u!(u, times[i], setup)
        NS.apply_bc_temp!(temp, times[i], setup)

        if haskey(results, :kinetic_energy)
            results[:kinetic_energy][i] = NS.total_kinetic_energy(u, setup)
        end
        if haskey(results, :angular_momentum_z)
            results[:angular_momentum_z][i] = angular_momentum_z(u, angular_cache)
        end
        if haskey(results, :angular_momentum_z_normalized)
            results[:angular_momentum_z_normalized][i] = normalized_angular_momentum_z(u, setup, angular_cache)
        end
        if haskey(results, :nusselt_bottom) || haskey(results, :nusselt_top)
            nu = plate_nusselt_numbers(temp, setup, cfg)
            haskey(results, :nusselt_bottom) && (results[:nusselt_bottom][i] = nu.bottom)
            haskey(results, :nusselt_top) && (results[:nusselt_top][i] = nu.top)
        end
        if haskey(results, :nusselt_volume)
            results[:nusselt_volume][i] = volume_nusselt_number(u, temp, setup, cfg)
        end
        if haskey(results, :mean_temperature)
            results[:mean_temperature][i] = mean_temperature(temp, setup, cfg)
        end
        if haskey(results, :mean_temperature_rate) ||
           haskey(results, :temperature_source) ||
           haskey(results, :heat_budget_residual)
            budget = temperature_heat_budget(u, temp, dtemp, setup, cfg)
            haskey(results, :mean_temperature_rate) &&
                (results[:mean_temperature_rate][i] = budget.rate)
            haskey(results, :temperature_source) &&
                (results[:temperature_source][i] = budget.source)
            haskey(results, :heat_budget_residual) &&
                (results[:heat_budget_residual][i] = budget.residual)
        end
    end

    return (; run_id, cfg, times, results)
end

function diagnostics_output_dir(run_id::String)
    outdir = joinpath(fom_dir(run_id), "plots", "diagnostics")
    mkpath(outdir)
    return outdir
end

function save_fom_diagnostic_plots(
    run_id::String;
    diagnostics = (
        :kinetic_energy,
        :angular_momentum_z,
        :angular_momentum_z_normalized,
        :nusselt_bottom,
        :nusselt_top,
        :nusselt_volume,
        :mean_temperature,
        :mean_temperature_rate,
        :temperature_source,
        :heat_budget_residual,
    ),
)
    data = compute_fom_diagnostics(run_id; diagnostics)
    outdir = diagnostics_output_dir(run_id)

    for (name, values) in data.results
        fig = Figure(size = (900, 500))
        ax = Axis(fig[1, 1]; title = string(name), xlabel = "t", ylabel = string(name))
        lines!(ax, data.times, values)
        save(joinpath(outdir, "$(name).png"), fig)
    end

    if haskey(data.results, :angular_momentum_z_normalized)
        fig = Figure(size = (1000, 900))
        ax1 = Axis(fig[1, 1]; title = "Normalized angular momentum", xlabel = "t", ylabel = "Lz normalized")
        lines!(ax1, data.times, data.results[:angular_momentum_z_normalized])
        hlines!(ax1, [0.0]; color = :black, linestyle = :dash)

        ax2 = Axis(fig[2, 1]; title = "Nusselt number", xlabel = "t", ylabel = "Nu")
        if haskey(data.results, :nusselt_volume)
            lines!(ax2, data.times, data.results[:nusselt_volume]; label = "volume")
        end
        if haskey(data.results, :nusselt_bottom)
            lines!(ax2, data.times, data.results[:nusselt_bottom]; label = "bottom")
        end
        if haskey(data.results, :nusselt_top)
            lines!(ax2, data.times, data.results[:nusselt_top]; label = "top")
        end
        axislegend(ax2; position = :rt)

        ax3 = Axis(fig[3, 1]; title = "Kinetic energy", xlabel = "t", ylabel = "KE")
        if haskey(data.results, :kinetic_energy)
            lines!(ax3, data.times, data.results[:kinetic_energy])
        end
        save(joinpath(outdir, "rotation_diagnostics.png"), fig)
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
    return "nx=$(cfg.nx), ny=$(cfg.ny), cfl=$(cfg.cfl), nu=$(cfg.viscosity), kappa=$(cfg.conductivity)"
end

function save_fom_diagnostics_compare(
    run_ids::AbstractVector{String};
    diagnostics = (
        :kinetic_energy,
        :angular_momentum_z,
        :angular_momentum_z_normalized,
        :nusselt_bottom,
        :nusselt_top,
        :nusselt_volume,
        :mean_temperature,
    ),
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
