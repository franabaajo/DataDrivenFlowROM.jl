using CairoMakie
using KernelAbstractions
using Statistics

include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Config.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Paths.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Setup.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "FOM.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Stats.jl"))

import IncompressibleNavierStokes as NS

const RUN_IDS = [
    "fom_nx256_ny256_cfl0p67_stxuniform_styuniform_nu0p00025_kappa0p00025_0fca1778",
    "fom_nx256_ny256_cfl0p67_stxuniform_styuniform_nu0p00025_kappa0p00025_10f4df3e",
    "fom_nx512_ny512_cfl0p67_stxuniform_styuniform_nu0p00025_kappa0p00025_24aad466",
]

const PAPER_BENCHMARK = [
    (Ra = 1e3, NuH = 1.0004, NuC = 1.0004, umax = 3.6464e-6, vmax = 3.9684e-6, y_umax = NaN, x_vmax = NaN, Numax_hot = NaN, x_Numax_hot = NaN),
    (Ra = 1e4, NuH = 2.1581, NuC = 2.1580, umax = 0.25228, vmax = 0.26369, y_umax = 0.8023, x_vmax = 0.8263, Numax_hot = 3.023, x_Numax_hot = 0.7183),
    (Ra = 1e5, NuH = 3.9103, NuC = 3.9103, umax = 0.34434, vmax = 0.37569, y_umax = 0.8636, x_vmax = 0.8973, Numax_hot = 6.065, x_Numax_hot = 0.6993),
    (Ra = 1e6, NuH = 6.3092, NuC = 6.3092, umax = 0.37088, vmax = 0.40600, y_umax = 0.9036, x_vmax = 0.9359, Numax_hot = 11.69, x_Numax_hot = 0.6448),
]

function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            vals = map(header) do key
                v = getproperty(row, Symbol(key))
                if v isa AbstractString
                    "\"" * replace(v, "\"" => "\"\"") * "\""
                elseif v isa Real
                    isfinite(v) ? string(v) : ""
                else
                    string(v)
                end
            end
            println(io, join(vals, ","))
        end
    end
    return path
end

function nearest_index(coords, inds, x)
    vals = collect(inds)
    _, k = findmin(abs.(coords[vals] .- x))
    return vals[k]
end

function max_abs_on_indices(field, inds, xcoords, ycoords)
    best = -Inf
    bestx = NaN
    besty = NaN
    for I in inds
        val = abs(field[I])
        if val > best
            best = val
            bestx = xcoords[I[1]]
            besty = ycoords[I[2]]
        end
    end
    return best, bestx, besty
end

function centerline_umax(u_cpu, setup)
    ux = @view u_cpu[:, :, 1]
    irange, jrange = setup.Iu[1].indices
    i = nearest_index(Array(setup.x[1]), irange, 0.5)
    best = -Inf
    besty = NaN
    for j in jrange
        val = abs(ux[i, j])
        if val > best
            best = val
            besty = Array(setup.xp[2])[j]
        end
    end
    return best, besty
end

function centerline_vmax(u_cpu, setup)
    uy = @view u_cpu[:, :, 2]
    irange, jrange = setup.Iu[2].indices
    j = nearest_index(Array(setup.x[2]), jrange, 0.5)
    best = -Inf
    bestx = NaN
    for i in irange
        val = abs(uy[i, j])
        if val > best
            best = val
            bestx = Array(setup.xp[1])[i]
        end
    end
    return best, bestx
end

function hot_wall_local_nusselt(temp, setup)
    temp_cpu = Array(temp)
    dy_bottom = sum(Array(setup.Δu[2][1:1]))
    local_nu = @. -(temp_cpu[:, 2] - temp_cpu[:, 1]) / dy_bottom
    inds = 2:(length(local_nu) - 1)
    vals = local_nu[inds]
    k = argmax(vals)
    idx = collect(inds)[k]
    return vals[k], Array(setup.xp[1])[idx]
end

function rel_drift(times, vals)
    length(times) < 2 && return NaN
    t = collect(times)
    v = collect(vals)
    denom = sum((t .- mean(t)).^2)
    denom == 0 && return NaN
    slope = sum((t .- mean(t)) .* (v .- mean(v))) / denom
    return slope * (maximum(t) - minimum(t)) / max(abs(mean(v)), eps(Float64))
end

function summarize_tail(times, vals; frac_start = 0.8)
    tail_start = frac_start * maximum(times)
    mask = times .>= tail_start
    ttail = times[mask]
    vtail = vals[mask]
    return (
        tail_start = tail_start,
        mean = mean(vtail),
        std = std(vtail),
        min = minimum(vtail),
        max = maximum(vtail),
        final = vals[end],
        rel_drift = rel_drift(ttail, vtail),
    )
end

function compute_run(run_id)
    cfg = load_fom_config_from_run(run_id)
    snapshots = load_fom_snapshots(run_id)
    setup = getsetup(cfg; backend = KernelAbstractions.CPU())
    psolver = NS.default_psolver(setup)
    state = initial_state(cfg, setup, psolver)
    u = state.u
    temp = state.temp
    times = collect(fom_times(cfg))
    K = length(times)

    rows = Vector{NamedTuple}(undef, K)
    for i in 1:K
        t = times[i]
        copy_state_vector_to_fields!(@view(snapshots.state[:, i]), u, temp)
        NS.apply_bc_u!(u, t, setup)
        NS.apply_bc_temp!(temp, t, setup)

        nu = plate_nusselt_numbers(temp, setup, cfg)
        u_cpu = Array(u)
        ux = @view u_cpu[:, :, 1]
        uy = @view u_cpu[:, :, 2]

        gux, gux_x, gux_y = max_abs_on_indices(ux, setup.Iu[1], Array(setup.x[1]), Array(setup.xp[2]))
        guy, guy_x, guy_y = max_abs_on_indices(uy, setup.Iu[2], Array(setup.xp[1]), Array(setup.x[2]))
        cux, cux_y = centerline_umax(u_cpu, setup)
        cvy, cvy_x = centerline_vmax(u_cpu, setup)
        numax, x_numax = hot_wall_local_nusselt(temp, setup)

        rows[i] = (
            run_id = run_id,
            t = t,
            NuH = nu.bottom,
            NuC = nu.top,
            Nu_abs_mismatch = abs(nu.bottom - nu.top),
            Numax_hot = numax,
            x_Numax_hot = x_numax,
            kinetic_energy = NS.total_kinetic_energy(u, setup),
            centerline_umax = cux,
            y_centerline_umax = cux_y,
            centerline_vmax = cvy,
            x_centerline_vmax = cvy_x,
            global_umax = gux,
            x_global_umax = gux_x,
            y_global_umax = gux_y,
            global_vmax = guy,
            x_global_vmax = guy_x,
            y_global_vmax = guy_y,
        )
    end

    return cfg, rows
end

function main()
    outdir = mkpath(joinpath(output_root(), "paper_comparison_ouertatani"))

    paper_header = [
        "Ra", "NuH", "NuC", "umax", "vmax", "y_umax", "x_vmax", "Numax_hot", "x_Numax_hot",
    ]
    write_csv(joinpath(outdir, "ouertatani_2008_table.csv"), paper_header, PAPER_BENCHMARK)

    ts_header = [
        "run_id", "t", "NuH", "NuC", "Nu_abs_mismatch", "Numax_hot", "x_Numax_hot",
        "kinetic_energy", "centerline_umax", "y_centerline_umax", "centerline_vmax",
        "x_centerline_vmax", "global_umax", "x_global_umax", "y_global_umax",
        "global_vmax", "x_global_vmax", "y_global_vmax",
    ]
    summary_rows = NamedTuple[]
    all_series = Dict{String,Vector{NamedTuple}}()

    for run_id in RUN_IDS
        if !isdir(fom_dir(run_id)) || !has_streamed_snapshots(fom_dir(run_id))
            @warn "Skipping missing/incomplete run" run_id
            continue
        end
        @info "Computing paper-style statistics" run_id
        cfg, rows = compute_run(run_id)
        all_series[run_id] = rows
        write_csv(joinpath(outdir, "time_series_$(run_id).csv"), ts_header, rows)

        times = [r.t for r in rows]
        nuh = [r.NuH for r in rows]
        nuc = [r.NuC for r in rows]
        numax = [r.Numax_hot for r in rows]
        ke = [r.kinetic_energy for r in rows]
        cux = [r.centerline_umax for r in rows]
        cvy = [r.centerline_vmax for r in rows]
        gux = [r.global_umax for r in rows]
        guy = [r.global_vmax for r in rows]
        mismatch = [r.Nu_abs_mismatch for r in rows]

        s_nuh = summarize_tail(times, nuh)
        s_nuc = summarize_tail(times, nuc)
        s_numax = summarize_tail(times, numax)
        s_ke = summarize_tail(times, ke)
        s_cux = summarize_tail(times, cux)
        s_cvy = summarize_tail(times, cvy)
        s_gux = summarize_tail(times, gux)
        s_guy = summarize_tail(times, guy)
        s_mismatch = summarize_tail(times, mismatch)

        push!(summary_rows, (
            run_id = run_id,
            nx = cfg.nx,
            ny = cfg.ny,
            x_length = cfg.xlims[2] - cfg.xlims[1],
            y_height = cfg.ylims[2] - cfg.ylims[1],
            tsim = cfg.tsim,
            nsnapshot = cfg.nsnapshot,
            dt_snapshot = cfg.tsim / cfg.nsnapshot,
            Ra = cfg.gravity * (cfg.temp_hot - cfg.temp_cold) * (cfg.ylims[2] - cfg.ylims[1])^3 / (cfg.viscosity * cfg.conductivity),
            Pr = cfg.viscosity / cfg.conductivity,
            tail_start = s_nuh.tail_start,
            NuH_tail_mean = s_nuh.mean,
            NuH_tail_std = s_nuh.std,
            NuH_tail_rel_drift = s_nuh.rel_drift,
            NuH_final = s_nuh.final,
            NuC_tail_mean = s_nuc.mean,
            NuC_tail_std = s_nuc.std,
            NuC_tail_rel_drift = s_nuc.rel_drift,
            NuC_final = s_nuc.final,
            Nu_abs_mismatch_tail_mean = s_mismatch.mean,
            Numax_hot_tail_mean = s_numax.mean,
            Numax_hot_tail_std = s_numax.std,
            Numax_hot_tail_rel_drift = s_numax.rel_drift,
            Numax_hot_final = s_numax.final,
            kinetic_energy_tail_mean = s_ke.mean,
            kinetic_energy_tail_std = s_ke.std,
            kinetic_energy_tail_rel_drift = s_ke.rel_drift,
            kinetic_energy_final = s_ke.final,
            centerline_umax_tail_mean = s_cux.mean,
            centerline_umax_tail_std = s_cux.std,
            centerline_umax_final = s_cux.final,
            centerline_vmax_tail_mean = s_cvy.mean,
            centerline_vmax_tail_std = s_cvy.std,
            centerline_vmax_final = s_cvy.final,
            global_umax_tail_mean = s_gux.mean,
            global_umax_tail_std = s_gux.std,
            global_umax_final = s_gux.final,
            global_vmax_tail_mean = s_guy.mean,
            global_vmax_tail_std = s_guy.std,
            global_vmax_final = s_guy.final,
        ))
    end

    summary_header = [
        "run_id", "nx", "ny", "x_length", "y_height", "tsim", "nsnapshot", "dt_snapshot",
        "Ra", "Pr", "tail_start", "NuH_tail_mean", "NuH_tail_std", "NuH_tail_rel_drift",
        "NuH_final", "NuC_tail_mean", "NuC_tail_std", "NuC_tail_rel_drift", "NuC_final",
        "Nu_abs_mismatch_tail_mean", "Numax_hot_tail_mean", "Numax_hot_tail_std",
        "Numax_hot_tail_rel_drift", "Numax_hot_final", "kinetic_energy_tail_mean",
        "kinetic_energy_tail_std", "kinetic_energy_tail_rel_drift", "kinetic_energy_final",
        "centerline_umax_tail_mean", "centerline_umax_tail_std", "centerline_umax_final",
        "centerline_vmax_tail_mean", "centerline_vmax_tail_std", "centerline_vmax_final",
        "global_umax_tail_mean", "global_umax_tail_std", "global_umax_final",
        "global_vmax_tail_mean", "global_vmax_tail_std", "global_vmax_final",
    ]
    write_csv(joinpath(outdir, "simulation_summary.csv"), summary_header, summary_rows)

    fig = Figure(size = (1100, 700))
    ax1 = Axis(fig[1, 1]; xlabel = "t", ylabel = "Nu", title = "Wall-averaged Nusselt numbers")
    ax2 = Axis(fig[2, 1]; xlabel = "t", ylabel = "Kinetic energy", title = "Kinetic energy")
    colors = Makie.wong_colors()
    for (k, run_id) in enumerate(RUN_IDS)
        haskey(all_series, run_id) || continue
        rows = all_series[run_id]
        t = [r.t for r in rows]
        color = colors[mod1(k, length(colors))]
        label = replace(run_id, "fom_" => "")
        lines!(ax1, t, [r.NuH for r in rows]; color, label = "NuH " * label)
        lines!(ax1, t, [r.NuC for r in rows]; color, linestyle = :dash, label = "NuC " * label)
        lines!(ax2, t, [r.kinetic_energy for r in rows]; color, label)
    end
    axislegend(ax1; position = :rt, nbanks = 1)
    axislegend(ax2; position = :lt, nbanks = 1)
    save(joinpath(outdir, "nusselt_and_ke_timeseries.png"), fig)

    fig2 = Figure(size = (1100, 700))
    ax3 = Axis(fig2[1, 1]; xlabel = "t", ylabel = "max abs velocity", title = "Paper-style centerline velocity maxima")
    ax4 = Axis(fig2[2, 1]; xlabel = "t", ylabel = "Nu max", title = "Hot-wall maximum local Nusselt number")
    for (k, run_id) in enumerate(RUN_IDS)
        haskey(all_series, run_id) || continue
        rows = all_series[run_id]
        t = [r.t for r in rows]
        color = colors[mod1(k, length(colors))]
        label = replace(run_id, "fom_" => "")
        lines!(ax3, t, [r.centerline_umax for r in rows]; color, label = "umax " * label)
        lines!(ax3, t, [r.centerline_vmax for r in rows]; color, linestyle = :dash, label = "vmax " * label)
        lines!(ax4, t, [r.Numax_hot for r in rows]; color, label)
    end
    axislegend(ax3; position = :rt, nbanks = 1)
    axislegend(ax4; position = :rt, nbanks = 1)
    save(joinpath(outdir, "velocity_and_local_nusselt_timeseries.png"), fig2)

    println("Saved comparison outputs in: ", outdir)
end

main()
