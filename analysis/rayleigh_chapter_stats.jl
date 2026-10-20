using CairoMakie
using JLD2
using KernelAbstractions
using LinearAlgebra
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
        ),
        Legend = (
            labelsize = 20,
            patchsize = (34, 18),
        ),
    ),
)

const THESIS_LINEWIDTH = 3.4
const THESIS_REFERENCE_LINEWIDTH = 3.8
const THESIS_TIME_TICK_STEP = 50.0

include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Config.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Paths.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Registry.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Setup.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "FOM.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "POD.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "ROM.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Plotting.jl"))
include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Stats.jl"))

const FOM_ID = "fom_nx256_ny256_cfl0p67_stxuniform_styuniform_nu0p00025_kappa0p00025_c5e2b71f"
const PREDICTION_FOM_ID = "fom_nx256_ny256_cfl0p67_stxuniform_styuniform_nu0p00025_kappa0p00025_84e447fc"
const PREDICTION_TRAINING_END = 150.0

const FOM_GRID_RUNS = [
    (label = "256x256", nx = 256, id = "fom_nx256_ny256_cfl0p67_stxuniform_styuniform_nu0p00025_kappa0p00025_c5e2b71f"),
    (label = "384x384", nx = 384, id = "fom_nx384_ny384_cfl0p67_stxuniform_styuniform_nu0p00025_kappa0p00025_cpu_0cfb09fc"),
    (label = "512x512", nx = 512, id = "fom_nx512_ny512_cfl0p67_stxuniform_styuniform_nu0p00025_kappa0p00025_cpu_0173147d"),
    (label = "768x768", nx = 768, id = "fom_nx768_ny768_cfl0p67_stxuniform_styuniform_nu0p00025_kappa0p00025_7002799e"),
]

const METHOD_RUNS = [
    (label = "Exact OpInf, scheduled", id = "rom_from_c5e2b71f_switched_lz_ts10p0_te150p0_tagauto_lz_aligned_transition_10_150_6ff90faa"),
    (label = "Intrusive Galerkin, scheduled", id = "rom_from_c5e2b71f_switched_lz_ts10p0_te150p0_tagauto_lz_aligned_transition_10_150_intrusive_2f1d39ed"),
    (label = "Trajectory OpInf, scheduled", id = "rom_from_c5e2b71f_switched_lz_ts10p0_te150p0_tagauto_lz_aligned_transition_10_150_trajopinf_fd2e68f4"),
    (label = "Exact OpInf, threshold", id = "rom_from_c5e2b71f_switched_lz_ts10p0_te150p0_tagauto_lz_aligned_transition_10_150_06ee3e6f"),
]

const MATCHED_RANK_RUNS = [
    (label = "Exact OpInf, matched ranks", id = "rom_from_c5e2b71f_switched_lz_ts10p0_te150p0_tagmatched_trajrank_exact_10_150_ce8c850d"),
    (label = "Intrusive Galerkin, matched ranks", id = "rom_from_c5e2b71f_switched_lz_ts10p0_te150p0_tagmatched_trajrank_intrusive_10_150_12d29989"),
    (label = "Trajectory OpInf, matched ranks", id = "rom_from_c5e2b71f_switched_lz_ts10p0_te150p0_tagauto_lz_aligned_transition_10_150_trajopinf_fd2e68f4"),
]

const RETENTION_RUNS = [
    (eta = 0.90, label = "eta = 0.90", id = "rom_from_c5e2b71f_switched_lz_ts10p0_te150p0_tagauto_lz_aligned_transition_10_150_e0p90_212f65d9"),
    (eta = 0.95, label = "eta = 0.95", id = "rom_from_c5e2b71f_switched_lz_ts10p0_te150p0_tagauto_lz_aligned_transition_10_150_e0p95_f9a9a6a3"),
    (eta = 0.97, label = "eta = 0.97", id = "rom_from_c5e2b71f_switched_lz_ts10p0_te150p0_tagauto_lz_aligned_transition_10_150_e0p97_efc8da04"),
    (eta = 0.99, label = "eta = 0.99", id = "rom_from_c5e2b71f_switched_lz_ts10p0_te150p0_tagauto_lz_aligned_transition_10_150_6ff90faa"),
]

const PREDICTION_RUN = (
    label = "Exact OpInf, scheduled prediction",
    id = "rom_from_84e447fc_switched_lz_ts10p0_te300p0_tagauto_lz_aligned_transition_10_150_predict_300_scheduled_091aca87",
)

const EXACT_LOCAL = [
    "rom_from_c5e2b71f_exactopinf_bglob_e0p99_ts10p0_te85p66362054663915_tagauto_lz_aligned_transition_10_150_local_pre_10p0_85p66362054663915_e2b595b5",
    "rom_from_c5e2b71f_exactopinf_bglob_e0p99_ts77p66362054663915_te100p38686024952808_tagauto_lz_aligned_transition_10_150_local_transition_77p66362054663915_100p38686024952808_91e6094a",
    "rom_from_c5e2b71f_exactopinf_bglob_e0p99_ts92p38686024952808_te150p0_tagauto_lz_aligned_transition_10_150_local_post_92p38686024952808_150p0_d52178b1",
]

const INTRUSIVE_LOCAL = [
    "rom_from_c5e2b71f_intrusive_bglob_e0p99_ts10p0_te85p66362054663915_tagauto_lz_aligned_transition_10_150_intrusive_local_pre_10p0_85p66362054663915_a50472ea",
    "rom_from_c5e2b71f_intrusive_bglob_e0p99_ts77p66362054663915_te100p38686024952808_tagauto_lz_aligned_transition_10_150_intrusive_local_transition_77p66362054663915_100p38686024952808_bab4d343",
    "rom_from_c5e2b71f_intrusive_bglob_e0p99_ts92p38686024952808_te150p0_tagauto_lz_aligned_transition_10_150_intrusive_local_post_92p38686024952808_150p0_c8d0e178",
]

analysis_dir() = mkpath(joinpath(@__DIR__, "..", "results", "tables", "rayleigh_benard"))
figure_dir() = mkpath(joinpath(@__DIR__, "..", "results", "figures"))

function set_time_axis_from_zero!(ax, times)
    tmax = THESIS_TIME_TICK_STEP * ceil(Int, maximum(times) / THESIS_TIME_TICK_STEP)
    ticks = collect(0:THESIS_TIME_TICK_STEP:tmax)
    xlims!(ax, 0.0, tmax)
    ax.xticks = ticks
    return ax
end

function set_angular_momentum_ticks!(ax)
    ax.yticks = [-0.5, 0.0, 0.5]
    return ax
end

function csv_escape(x)
    s = string(x)
    if occursin(",", s) || occursin("\"", s) || occursin("\n", s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end

function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join(csv_escape.(row), ","))
        end
    end
end

function rel_state_error(x_rom, x_fom)
    sdiff = 0.0
    sfom = 0.0
    @inbounds for i in eachindex(x_rom, x_fom)
        xi = Float64(x_fom[i])
        di = Float64(x_rom[i]) - xi
        sdiff += di * di
        sfom += xi * xi
    end
    return sqrt(sdiff) / max(sqrt(sfom), eps(Float64))
end

function rel_series_error(y, yref)
    den = sqrt(mean(abs2, yref))
    den <= eps(Float64) && return sqrt(mean(abs2, y .- yref))
    return sqrt(mean(abs2, y .- yref)) / den
end

function load_local_models(fit)
    names = hasproperty(fit, :local_model_names) ? fit.local_model_names : fit.local_rom_ids
    return [_load_local_rom_model(run_id, names[i]) for (i, run_id) in enumerate(fit.local_rom_ids)]
end

function reconstruct_state!(xrom, fit, traj, k, local_models)
    if hasproperty(fit, :switched) && fit.switched
        model = local_models[Int(traj.active_model_index[k])]
        r = model.rfit
        _reconstruct_local_state!(xrom, model, @view(traj.a_rom[1:r, k]))
        return xrom
    end
    mul!(xrom, fit.rom_basis, @view(traj.a_rom[:, k]))
    if hasproperty(fit, :basis_state_mean) && !isnothing(fit.basis_state_mean)
        xrom .+= fit.basis_state_mean
    end
    return xrom
end

function run_diagnostics(run_id, snapshots, setup, fomcfg, state, angular_cache)
    meta = load_rom_metadata(run_id)
    fit = load_rom_fit(run_id)
    traj = load_rom_trajectory(run_id)
    local_models = hasproperty(fit, :switched) && fit.switched ? load_local_models(fit) : nothing
    rom_basis = isnothing(local_models) ? fit.rom_basis : local_models[1].rom_basis
    xrom = Vector{Float64}(undef, size(rom_basis, 1))

    alignment = rom_trajectory_alignment(traj, fomcfg, snapshots)
    K = alignment.K
    times = alignment.times
    indices = alignment.fom_snapshot_indices

    u_fom = state.u
    temp_fom = state.temp
    u_rom = similar(u_fom)
    temp_rom = similar(temp_fom)

    err = zeros(Float64, K)
    ke_fom = zeros(Float64, K)
    ke_rom = zeros(Float64, K)
    nu_bottom_fom = zeros(Float64, K)
    nu_bottom_rom = zeros(Float64, K)
    nu_top_fom = zeros(Float64, K)
    nu_top_rom = zeros(Float64, K)
    lz_fom = zeros(Float64, K)
    lz_rom = zeros(Float64, K)

    for k in 1:K
        t = times[k]
        x_fom = @view snapshots.state[:, indices[k]]
        xk = reconstruct_state!(xrom, fit, traj, k, local_models)
        err[k] = rel_state_error(xk, x_fom)

        copy_state_vector_to_fields!(x_fom, u_fom, temp_fom)
        copy_state_vector_to_fields!(xk, u_rom, temp_rom)
        NS.apply_bc_u!(u_fom, t, setup)
        NS.apply_bc_u!(u_rom, t, setup)
        NS.apply_bc_temp!(temp_fom, t, setup)
        NS.apply_bc_temp!(temp_rom, t, setup)

        ke_fom[k] = NS.total_kinetic_energy(u_fom, setup)
        ke_rom[k] = NS.total_kinetic_energy(u_rom, setup)
        nub_fom = plate_nusselt_numbers(temp_fom, setup, fomcfg)
        nub_rom = plate_nusselt_numbers(temp_rom, setup, fomcfg)
        nu_bottom_fom[k] = nub_fom.bottom
        nu_bottom_rom[k] = nub_rom.bottom
        nu_top_fom[k] = nub_fom.top
        nu_top_rom[k] = nub_rom.top
        lz_fom[k] = normalized_angular_momentum_z(u_fom, setup, angular_cache)
        lz_rom[k] = normalized_angular_momentum_z(u_rom, setup, angular_cache)
    end

    active_counts = if hasproperty(traj, :active_model_index)
        [count(==(i), traj.active_model_index[1:K]) for i in 1:maximum(traj.active_model_index[1:K])]
    else
        Int[]
    end

    return (;
        run_id,
        meta,
        fit,
        traj,
        times,
        rel_error = err,
        ke_fom,
        ke_rom,
        nu_bottom_fom,
        nu_bottom_rom,
        nu_top_fom,
        nu_top_rom,
        lz_fom,
        lz_rom,
        active_counts,
        lz_rmse = rel_series_error(lz_rom, lz_fom),
        ke_rmse = rel_series_error(ke_rom, ke_fom),
        nu_bottom_rmse = rel_series_error(nu_bottom_rom, nu_bottom_fom),
        nu_top_rmse = rel_series_error(nu_top_rom, nu_top_fom),
    )
end

function summarize_run(label, data)
    e = data.rel_error
    m = data.meta
    runtime_seconds = hasproperty(m, :runtime_seconds) ? m.runtime_seconds : NaN
    advance_runtime_seconds = hasproperty(m, :advance_runtime_seconds) ? m.advance_runtime_seconds : NaN
    setup_and_online_runtime_seconds = hasproperty(m, :setup_and_online_runtime_seconds) ?
        m.setup_and_online_runtime_seconds : NaN
    return [
        label,
        data.run_id,
        join(m.local_rfits, "/"),
        m.rfit,
        round(runtime_seconds; digits = 3),
        round(advance_runtime_seconds; digits = 3),
        round(setup_and_online_runtime_seconds; digits = 3),
        round(mean(e); digits = 5),
        round(median(e); digits = 5),
        round(maximum(e); digits = 5),
        round(last(e); digits = 5),
        round(data.ke_rmse; digits = 5),
        round(data.lz_rmse; digits = 5),
        round(data.nu_bottom_rmse; digits = 5),
        round(data.nu_top_rmse; digits = 5),
        join(data.active_counts, "/"),
        join(round.(getproperty.(data.meta.switch_events, :time); digits = 3), "/"),
    ]
end

function interval_summary_row(label, data, idx)
    if isempty(idx)
        return [
            label, "n/a", 0,
            "NaN", "NaN", "NaN", "NaN", "NaN", "NaN", "NaN",
        ]
    end
    e = data.rel_error[idx]
    return [
        label,
        @sprintf("%.3f--%.3f", first(data.times[idx]), last(data.times[idx])),
        length(idx),
        @sprintf("%.5f", mean(e)),
        @sprintf("%.5f", median(e)),
        @sprintf("%.5f", maximum(e)),
        @sprintf("%.5f", last(e)),
        @sprintf("%.5f", rel_series_error(data.ke_rom[idx], data.ke_fom[idx])),
        @sprintf("%.5f", rel_series_error(data.lz_rom[idx], data.lz_fom[idx])),
        @sprintf("%.5f", rel_series_error(data.nu_bottom_rom[idx], data.nu_bottom_fom[idx])),
    ]
end

function prediction_interval_rows(data)
    ttrain = PREDICTION_TRAINING_END
    all_idx = collect(eachindex(data.times))
    training_idx = findall(t -> t <= ttrain, data.times)
    prediction_idx = findall(t -> t > ttrain, data.times)
    return [
        interval_summary_row("full run", data, all_idx),
        interval_summary_row("training-window comparison", data, training_idx),
        interval_summary_row("post-training prediction", data, prediction_idx),
    ]
end

function operator_comparison_rows()
    regimes = ["pre", "transition", "post"]
    rows = Any[]
    for i in eachindex(regimes)
        ef = load_rom_fit(EXACT_LOCAL[i])
        inf = load_rom_fit(INTRUSIVE_LOCAL[i])
        r = size(ef.L, 1)
        c_abs = norm(ef.c - inf.c)
        c_rel = c_abs / max(norm(inf.c), eps(Float64))
        L_abs = norm(ef.L - inf.L)
        L_rel = L_abs / max(norm(inf.L), eps(Float64))
        Q_abs = norm(vec(ef.Q) - vec(inf.Q))
        Q_rel = Q_abs / max(norm(vec(inf.Q)), eps(Float64))
        push!(rows, [
            regimes[i],
            r,
            @sprintf("%.3e", c_abs),
            @sprintf("%.3e", c_rel),
            @sprintf("%.3e", L_abs),
            @sprintf("%.3e", L_rel),
            @sprintf("%.3e", Q_abs),
            @sprintf("%.3e", Q_rel),
        ])
    end
    return rows
end

function plot_method_comparison(datas)
    figs = figure_dir()
    colors = Makie.wong_colors()
    exact = datas["Exact OpInf, scheduled"]
    method_styles = Dict(
        "Exact OpInf, scheduled" => (color = colors[1], linestyle = :solid, linewidth = THESIS_LINEWIDTH),
        "Intrusive Galerkin, scheduled" => (color = colors[2], linestyle = :dash, linewidth = THESIS_LINEWIDTH),
        "Trajectory OpInf, scheduled" => (color = colors[3], linestyle = :solid, linewidth = THESIS_LINEWIDTH),
    )
    method_keys = ["Exact OpInf, scheduled", "Intrusive Galerkin, scheduled", "Trajectory OpInf, scheduled"]

    fig = Figure(size = (1150, 720))
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "relative state error", yscale = log10)
    for key in method_keys
        d = datas[key]
        lines!(ax, d.times, d.rel_error; label = key, method_styles[key]...)
    end
    for ev in exact.meta.switch_events
        vlines!(ax, [ev.time]; color = :gray35, linestyle = :dash, linewidth = 2.6)
    end
    set_time_axis_from_zero!(ax, exact.times)
    axislegend(ax; position = :lt)
    save(joinpath(figs, "rayleigh_method_relative_error.png"), fig)

    fig = Figure(size = (1150, 760))
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "kinetic energy")
    lines!(ax, exact.times, exact.ke_fom; label = "FOM", color = :black, linewidth = THESIS_REFERENCE_LINEWIDTH)
    for key in method_keys
        d = datas[key]
        lines!(ax, d.times, d.ke_rom; label = key, method_styles[key]...)
    end
    for ev in exact.meta.switch_events
        vlines!(ax, [ev.time]; color = :gray35, linestyle = :dash, linewidth = 2.6)
    end
    set_time_axis_from_zero!(ax, exact.times)
    axislegend(ax; position = :rb)
    save(joinpath(figs, "rayleigh_method_kinetic_energy.png"), fig)

    fig = Figure(size = (1150, 760))
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "normalized angular momentum")
    lines!(ax, exact.times, exact.lz_fom; label = "FOM", color = :black, linewidth = THESIS_REFERENCE_LINEWIDTH)
    for key in method_keys
        d = datas[key]
        lines!(ax, d.times, d.lz_rom; label = key, method_styles[key]...)
    end
    hlines!(ax, [0.0]; color = :gray50, linestyle = :dot, linewidth = 2.4)
    for ev in exact.meta.switch_events
        vlines!(ax, [ev.time]; color = :gray35, linestyle = :dash, linewidth = 2.6)
    end
    set_time_axis_from_zero!(ax, exact.times)
    set_angular_momentum_ticks!(ax)
    axislegend(ax; position = :rb)
    save(joinpath(figs, "rayleigh_method_angular_momentum.png"), fig)
end

function plot_matched_rank_comparison(datas)
    figs = figure_dir()
    colors = Makie.wong_colors()
    exact = datas["Exact OpInf, matched ranks"]
    method_styles = Dict(
        "Exact OpInf, matched ranks" => (color = colors[1], linestyle = :solid, linewidth = THESIS_LINEWIDTH),
        "Intrusive Galerkin, matched ranks" => (color = colors[2], linestyle = :dash, linewidth = THESIS_LINEWIDTH),
        "Trajectory OpInf, matched ranks" => (color = colors[3], linestyle = :solid, linewidth = THESIS_LINEWIDTH),
    )
    method_keys = [
        "Exact OpInf, matched ranks",
        "Intrusive Galerkin, matched ranks",
        "Trajectory OpInf, matched ranks",
    ]

    fig = Figure(size = (1150, 720))
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "relative state error", yscale = log10)
    for key in method_keys
        d = datas[key]
        lines!(ax, d.times, d.rel_error; label = key, method_styles[key]...)
    end
    for ev in exact.meta.switch_events
        vlines!(ax, [ev.time]; color = :gray35, linestyle = :dash, linewidth = 2.6)
    end
    set_time_axis_from_zero!(ax, exact.times)
    axislegend(ax; position = :lt)
    save(joinpath(figs, "rayleigh_matched_rank_relative_error.png"), fig)
end

function plot_retention_comparison(datas)
    figs = figure_dir()
    colors = Makie.wong_colors()
    e99 = datas[0.99]

    fig = Figure(size = (1150, 720))
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "relative state error", yscale = log10)
    for (i, item) in enumerate(RETENTION_RUNS)
        d = datas[item.eta]
        lines!(ax, d.times, d.rel_error; label = item.label, color = colors[i], linewidth = THESIS_LINEWIDTH)
    end
    for ev in e99.meta.switch_events
        vlines!(ax, [ev.time]; color = :gray35, linestyle = :dash, linewidth = 2.6)
    end
    set_time_axis_from_zero!(ax, e99.times)
    axislegend(ax; position = :lt)
    save(joinpath(figs, "rayleigh_retention_relative_error.png"), fig)

    fig = Figure(size = (1150, 760))
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "normalized angular momentum")
    lines!(ax, e99.times, e99.lz_fom; label = "FOM", color = :black, linewidth = THESIS_REFERENCE_LINEWIDTH)
    for (i, item) in enumerate(RETENTION_RUNS)
        d = datas[item.eta]
        lines!(ax, d.times, d.lz_rom; label = item.label, color = colors[i], linewidth = THESIS_LINEWIDTH)
    end
    hlines!(ax, [0.0]; color = :gray50, linestyle = :dot, linewidth = 2.4)
    for ev in e99.meta.switch_events
        vlines!(ax, [ev.time]; color = :gray35, linestyle = :dash, linewidth = 2.6)
    end
    set_time_axis_from_zero!(ax, e99.times)
    set_angular_momentum_ticks!(ax)
    axislegend(ax; position = :lb)
    save(joinpath(figs, "rayleigh_retention_angular_momentum.png"), fig)

    etas = [item.eta for item in RETENTION_RUNS]
    rmax = [datas[item.eta].meta.rfit for item in RETENTION_RUNS]
    advance = [datas[item.eta].meta.advance_runtime_seconds for item in RETENTION_RUNS]

    fig = Figure(size = (1150, 720))
    ax1 = Axis(fig[1, 1]; xlabel = "retained fluctuation energy", ylabel = "maximum local rank")
    scatterlines!(ax1, etas, rmax; color = colors[1], marker = :circle, markersize = 18, linewidth = THESIS_LINEWIDTH, label = "rank")
    ax2 = Axis(fig[1, 1]; yaxisposition = :right, ylabel = "online advance [s]")
    hidespines!(ax2, :l, :b, :t)
    hidexdecorations!(ax2)
    scatterlines!(ax2, etas, advance; color = colors[3], marker = :utriangle, markersize = 18, linewidth = THESIS_LINEWIDTH, label = "online advance")
    axislegend(ax1; position = :lt)
    axislegend(ax2; position = :rb)
    save(joinpath(figs, "rayleigh_retention_rank_runtime.png"), fig)
end

function plot_prediction(data)
    figs = figure_dir()
    switch_times = getproperty.(data.meta.switch_events, :time)

    fig = Figure(size = (1150, 720))
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "relative state error", yscale = log10)
    lines!(ax, data.times, data.rel_error; label = "ROM vs full-time FOM", color = Makie.wong_colors()[1], linewidth = THESIS_LINEWIDTH)
    for t in switch_times
        vlines!(ax, [t]; color = :gray35, linestyle = :dash, linewidth = 2.6)
    end
    vlines!(ax, [PREDICTION_TRAINING_END]; color = :red3, linestyle = :dot, linewidth = 3.0)
    axislegend(ax; position = :lt)
    save(joinpath(figs, "rayleigh_prediction_relative_error.png"), fig)

    fig = Figure(size = (1150, 980))
    ax1 = Axis(fig[1, 1]; xlabel = "time", ylabel = "kinetic energy")
    ax2 = Axis(fig[2, 1]; xlabel = "time", ylabel = "normalized angular momentum")
    ax3 = Axis(fig[3, 1]; xlabel = "time", ylabel = "bottom Nusselt number")

    lines!(ax1, data.times, data.ke_fom; label = "FOM", color = :black, linewidth = THESIS_REFERENCE_LINEWIDTH)
    lines!(ax1, data.times, data.ke_rom; label = "ROM", color = Makie.wong_colors()[1], linewidth = THESIS_LINEWIDTH)
    lines!(ax2, data.times, data.lz_fom; label = "FOM", color = :black, linewidth = THESIS_REFERENCE_LINEWIDTH)
    lines!(ax2, data.times, data.lz_rom; label = "ROM", color = Makie.wong_colors()[1], linewidth = THESIS_LINEWIDTH)
    lines!(ax3, data.times, data.nu_bottom_fom; label = "FOM", color = :black, linewidth = THESIS_REFERENCE_LINEWIDTH)
    lines!(ax3, data.times, data.nu_bottom_rom; label = "ROM", color = Makie.wong_colors()[1], linewidth = THESIS_LINEWIDTH)

    hlines!(ax2, [0.0]; color = :gray50, linestyle = :dot, linewidth = 2.4)
    for ax in (ax1, ax2, ax3)
        for t in switch_times
            vlines!(ax, [t]; color = :gray35, linestyle = :dash, linewidth = 2.6)
        end
        vlines!(ax, [PREDICTION_TRAINING_END]; color = :red3, linestyle = :dot, linewidth = 3.0)
    end
    axislegend(ax1; position = :rt)
    save(joinpath(figs, "rayleigh_prediction_diagnostics.png"), fig)
end

function common_time_indices(datasets)
    sets = [Set(round.(data.times; digits = 8)) for data in datasets]
    common = sort!(collect(intersect(sets...)))
    maps = [
        Dict(round(t; digits = 8) => i for (i, t) in enumerate(data.times))
        for data in datasets
    ]
    return common, [[maps[j][t] for t in common] for j in eachindex(datasets)]
end

function metric_values(data, key::Symbol)
    return data.results[key]
end

function relative_common_error(data, ref, idx, idxref, key::Symbol)
    y = metric_values(data, key)[idx]
    yref = metric_values(ref, key)[idxref]
    return rel_series_error(y, yref)
end

function fom_grid_rows(datasets)
    common_times, common_indices = common_time_indices(datasets)
    ref = datasets[end]
    idxref = common_indices[end]
    rows = Any[]
    for (i, item) in enumerate(FOM_GRID_RUNS)
        data = datasets[i]
        idx = common_indices[i]
        tail = findall(t -> t >= 100.0, data.times)
        push!(rows, [
            item.label,
            item.nx,
            data.cfg.nsnapshot,
            round(data.times[2] - data.times[1]; digits = 4),
            length(common_times),
            @sprintf("%.6f", last(metric_values(data, :kinetic_energy))),
            @sprintf("%.6f", last(metric_values(data, :angular_momentum_z_normalized))),
            @sprintf("%.6f", mean(metric_values(data, :nusselt_bottom)[tail])),
            @sprintf("%.6f", mean(metric_values(data, :nusselt_top)[tail])),
            @sprintf("%.6f", mean(metric_values(data, :nusselt_volume)[tail])),
            i == length(datasets) ? "0.00000" : @sprintf("%.5f", relative_common_error(data, ref, idx, idxref, :kinetic_energy)),
            i == length(datasets) ? "0.00000" : @sprintf("%.5f", relative_common_error(data, ref, idx, idxref, :angular_momentum_z_normalized)),
            i == length(datasets) ? "0.00000" : @sprintf("%.5f", relative_common_error(data, ref, idx, idxref, :nusselt_bottom)),
            i == length(datasets) ? "0.00000" : @sprintf("%.5f", relative_common_error(data, ref, idx, idxref, :nusselt_top)),
            i == length(datasets) ? "0.00000" : @sprintf("%.5f", relative_common_error(data, ref, idx, idxref, :nusselt_volume)),
        ])
    end
    return rows
end

function plot_fom_grid_convergence(datasets)
    figs = figure_dir()
    colors = Makie.wong_colors()
    fig = Figure(size = (1150, 980))
    ax1 = Axis(fig[1, 1]; xlabel = "time", ylabel = "kinetic energy")
    ax2 = Axis(fig[2, 1]; xlabel = "time", ylabel = "normalized angular momentum")
    ax3 = Axis(fig[3, 1]; xlabel = "time", ylabel = "volume Nusselt number")

    for (i, item) in enumerate(FOM_GRID_RUNS)
        data = datasets[i]
        color = colors[mod1(i, length(colors))]
        lines!(ax1, data.times, metric_values(data, :kinetic_energy); label = item.label, color, linewidth = THESIS_LINEWIDTH)
        lines!(ax2, data.times, metric_values(data, :angular_momentum_z_normalized); label = item.label, color, linewidth = THESIS_LINEWIDTH)
        lines!(ax3, data.times, metric_values(data, :nusselt_volume); label = item.label, color, linewidth = THESIS_LINEWIDTH)
    end
    axislegend(ax1; position = :lt)
    hlines!(ax2, [0.0]; color = :gray50, linestyle = :dot, linewidth = 2.4)
    save(joinpath(figs, "rayleigh_fom_grid_convergence.png"), fig)
end

function main()
    out = analysis_dir()
    fom_grid_data = Any[]
    for item in FOM_GRID_RUNS
        @info "Computing FOM grid diagnostics" item.label item.id
        push!(fom_grid_data, compute_fom_diagnostics(
            item.id;
            diagnostics = (
                :kinetic_energy,
                :angular_momentum_z_normalized,
                :nusselt_bottom,
                :nusselt_top,
                :nusselt_volume,
            ),
        ))
    end
    write_csv(
        joinpath(out, "fom_grid_convergence.csv"),
        [
            "grid", "nx", "nsnapshot", "saved_dt", "n_common_times",
            "final_ke", "final_lz_norm", "mean_t100_150_nu_bottom",
            "mean_t100_150_nu_top", "mean_t100_150_nu_volume",
            "rel_ke_error_vs_reference", "rel_lz_error_vs_reference",
            "rel_nu_bottom_error_vs_reference", "rel_nu_top_error_vs_reference",
            "rel_nu_volume_error_vs_reference",
        ],
        fom_grid_rows(fom_grid_data),
    )
    plot_fom_grid_convergence(fom_grid_data)

    fomcfg = load_fom_config_from_run(FOM_ID)
    snapshots = load_fom_snapshots(FOM_ID)
    setup = getsetup(fomcfg; backend = KernelAbstractions.CPU())
    psolver = NS.default_psolver(setup)
    state = initial_state(fomcfg, setup, psolver)
    angular_cache = angular_momentum_z_cache(setup, state.u, fomcfg)

    method_data = Dict{String,Any}()
    for item in METHOD_RUNS
        @info "Computing diagnostics" item.label item.id
        method_data[item.label] = run_diagnostics(item.id, snapshots, setup, fomcfg, state, angular_cache)
    end

    matched_rank_data = Dict{String,Any}()
    for item in MATCHED_RANK_RUNS
        @info "Computing matched-rank diagnostics" item.label item.id
        matched_rank_data[item.label] = run_diagnostics(item.id, snapshots, setup, fomcfg, state, angular_cache)
    end

    retention_data = Dict{Float64,Any}()
    for item in RETENTION_RUNS
        @info "Computing retention diagnostics" item.eta item.id
        retention_data[item.eta] = run_diagnostics(item.id, snapshots, setup, fomcfg, state, angular_cache)
    end

    method_rows = [summarize_run(item.label, method_data[item.label]) for item in METHOD_RUNS]
    write_csv(
        joinpath(out, "method_comparison.csv"),
        [
            "label", "run_id", "local_ranks", "rmax", "runtime_seconds",
            "advance_runtime_seconds", "setup_online_seconds", "mean_rel_state_error",
            "median_rel_state_error", "max_rel_state_error", "final_rel_state_error",
            "rel_ke_rmse", "rel_lz_rmse", "rel_nu_bottom_rmse", "rel_nu_top_rmse",
            "active_counts", "switch_times",
        ],
        method_rows,
    )

    matched_rank_rows = [
        summarize_run(item.label, matched_rank_data[item.label])
        for item in MATCHED_RANK_RUNS
    ]
    write_csv(
        joinpath(out, "matched_rank_comparison.csv"),
        [
            "label", "run_id", "local_ranks", "rmax", "runtime_seconds",
            "advance_runtime_seconds", "setup_online_seconds", "mean_rel_state_error",
            "median_rel_state_error", "max_rel_state_error", "final_rel_state_error",
            "rel_ke_rmse", "rel_lz_rmse", "rel_nu_bottom_rmse", "rel_nu_top_rmse",
            "active_counts", "switch_times",
        ],
        matched_rank_rows,
    )

    retention_rows = [summarize_run("eta=$(item.eta)", retention_data[item.eta]) for item in RETENTION_RUNS]
    write_csv(
        joinpath(out, "retention_comparison.csv"),
        [
            "label", "run_id", "local_ranks", "rmax", "runtime_seconds",
            "advance_runtime_seconds", "setup_online_seconds", "mean_rel_state_error",
            "median_rel_state_error", "max_rel_state_error", "final_rel_state_error",
            "rel_ke_rmse", "rel_lz_rmse", "rel_nu_bottom_rmse", "rel_nu_top_rmse",
            "active_counts", "switch_times",
        ],
        retention_rows,
    )

    write_csv(
        joinpath(out, "operator_exact_vs_intrusive.csv"),
        ["regime", "rank", "c_abs_fro", "c_rel_fro", "L_abs_fro", "L_rel_fro", "Q_abs_fro", "Q_rel_fro"],
        operator_comparison_rows(),
    )

    plot_method_comparison(method_data)
    plot_matched_rank_comparison(matched_rank_data)
    plot_retention_comparison(retention_data)

    prediction_fomcfg = load_fom_config_from_run(PREDICTION_FOM_ID)
    prediction_snapshots = load_fom_snapshots(PREDICTION_FOM_ID)
    prediction_setup = getsetup(prediction_fomcfg; backend = KernelAbstractions.CPU())
    prediction_psolver = NS.default_psolver(prediction_setup)
    prediction_state = initial_state(prediction_fomcfg, prediction_setup, prediction_psolver)
    prediction_angular_cache = angular_momentum_z_cache(
        prediction_setup,
        prediction_state.u,
        prediction_fomcfg,
    )
    @info "Computing prediction diagnostics" PREDICTION_RUN.label PREDICTION_RUN.id
    prediction_data = run_diagnostics(
        PREDICTION_RUN.id,
        prediction_snapshots,
        prediction_setup,
        prediction_fomcfg,
        prediction_state,
        prediction_angular_cache,
    )
    write_csv(
        joinpath(out, "prediction_comparison.csv"),
        [
            "interval", "time_range", "n_points", "mean_rel_state_error",
            "median_rel_state_error", "max_rel_state_error", "final_rel_state_error",
            "rel_ke_rmse", "rel_lz_rmse", "rel_nu_bottom_rmse",
        ],
        prediction_interval_rows(prediction_data),
    )
    plot_prediction(prediction_data)

    println("Wrote analysis tables to: ", out)
    println("Wrote thesis figures to: ", figure_dir())
end

main()
