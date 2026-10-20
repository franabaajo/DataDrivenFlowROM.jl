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

include(joinpath(@__DIR__, "..", "src", "dipole", "Config.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "Paths.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "Registry.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "Setup.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "FOM.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "POD.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "ROM.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "Stats.jl"))

const FOM_ID = "fom_n1536_cfl0p67_st1p2_ns50_nu0p0004_e592f3a8"

const METHOD_RUNS = [
    (label = "Exact OpInf, eta=0.95", id = "rom_from_e592f3a8_exactopinf_bglob_e0p95_d9547ec6"),
    (label = "Exact OpInf, eta=0.99", id = "rom_from_e592f3a8_exactopinf_bglob_e0p99_87d75a35"),
    (label = "Intrusive Galerkin, eta=0.99", id = "rom_from_e592f3a8_intrusive_bglob_e0p99_tagdirect_intrusive_e12a4370"),
    (label = "Trajectory OpInf, eta=0.99 request", id = "rom_from_e592f3a8_trajopinf_bglob_e0p99_lam0p01_tagtrajopinf_e0p99_2b364ac5"),
]

analysis_dir() = mkpath(joinpath(@__DIR__, "..", "results", "tables", "dipole"))
figure_dir() = mkpath(joinpath(@__DIR__, "..", "results", "figures"))

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
        f = Float64(x_fom[i])
        d = Float64(x_rom[i]) - f
        sdiff += d * d
        sfom += f * f
    end
    return sqrt(sdiff) / max(sqrt(sfom), eps(Float64))
end

function rel_series_error(y, yref)
    den = sqrt(mean(abs2, yref))
    den <= eps(Float64) && return sqrt(mean(abs2, y .- yref))
    return sqrt(mean(abs2, y .- yref)) / den
end

function peak_in_window(times, y, lo, hi)
    idx = findall(t -> lo <= t <= hi, times)
    isempty(idx) && error("No samples in peak window [$lo,$hi]")
    local_idx = argmax(@view(y[idx]))
    k = idx[local_idx]

    if firstindex(times) < k < lastindex(times)
        ts = Float64[times[k - 1], times[k], times[k + 1]]
        ys = Float64[y[k - 1], y[k], y[k + 1]]
        coeffs = [ts .^ 2 ts ones(3)] \ ys
        a, b, c = coeffs
        tpeak = -b / (2a)
        if isfinite(tpeak) && a < 0 && ts[1] <= tpeak <= ts[3] && lo <= tpeak <= hi
            zpeak = a * tpeak^2 + b * tpeak + c
            return (; index = k, time = tpeak, value = zpeak, interpolated = true)
        end
    end

    return (; index = k, time = times[k], value = y[k], interpolated = false)
end

function velocity_to_vorticity_vector!(u, xvec, setup, t)
    copyto!(vec(u), xvec)
    NS.apply_bc_u!(u, t, setup)
    vort = NS.vorticity(u, setup) |> Array
    return vec(vort[setup.Ip])
end

function near_wall_mask(setup)
    coords = getindex.(Array.(setup.xp), setup.Ip.indices)
    x = coords[1]
    y = coords[2]
    return vec([0.4 <= xi <= 1.0 && 0.0 <= yj <= 0.6 for xi in x, yj in y])
end

function spatial_correlation(a, b)
    aa = a .- mean(a)
    bb = b .- mean(b)
    return dot(aa, bb) / max(norm(aa) * norm(bb), eps(Float64))
end

function summarize_run(label, run_id, snapshots, setup, u_fom, u_rom, mask, energy_by_rank)
    meta = load_rom_metadata(run_id)
    fit = load_rom_fit(run_id)
    traj = load_rom_trajectory(run_id)
    cfg = load_fom_config_from_run(FOM_ID)
    times = collect(range(0.0, cfg.tsim; length = min(size(traj.a_rom, 2), size(snapshots.u, 2))))
    K = length(times)

    u_rom_vec = zeros(eltype(fit.rom_basis), size(fit.rom_basis, 1))
    x_fom = vec(u_fom)
    x_rom = vec(u_rom)

    err = zeros(Float64, K)
    ke_fom = zeros(Float64, K)
    ke_rom = zeros(Float64, K)
    en_fom = zeros(Float64, K)
    en_rom = zeros(Float64, K)

    final_vort_fom = Float64[]
    final_vort_rom = Float64[]

    for k in 1:K
        t = times[k]
        copyto!(x_fom, @view(snapshots.u[:, k]))
        mul!(u_rom_vec, fit.rom_basis, @view(traj.a_rom[:, k]))
        copyto!(x_rom, u_rom_vec)

        err[k] = rel_state_error(x_rom, x_fom)

        NS.apply_bc_u!(u_fom, t, setup)
        NS.apply_bc_u!(u_rom, t, setup)
        ke_fom[k] = NS.total_kinetic_energy(u_fom, setup)
        ke_rom[k] = NS.total_kinetic_energy(u_rom, setup)
        en_fom[k] = total_enstrophy(u_fom, setup)
        en_rom[k] = total_enstrophy(u_rom, setup)

        if k == K
            final_vort_fom = vec((NS.vorticity(u_fom, setup) |> Array)[setup.Ip])
            final_vort_rom = vec((NS.vorticity(u_rom, setup) |> Array)[setup.Ip])
        end
    end

    p1_fom = peak_in_window(times, en_fom, 0.20, 0.45)
    p2_fom = peak_in_window(times, en_fom, 0.45, 0.80)
    p1_rom = peak_in_window(times, en_rom, 0.20, 0.45)
    p2_rom = peak_in_window(times, en_rom, 0.45, 0.80)

    near_fom = final_vort_fom[mask]
    near_rom = final_vort_rom[mask]
    near_en_fom = 0.5 * mean(abs2, near_fom)
    near_en_rom = 0.5 * mean(abs2, near_rom)

    energy_retained = get(energy_by_rank, Int(meta.rfit), NaN)

    return (;
        label,
        run_id,
        meta,
        fit,
        traj,
        times,
        rel_error = err,
        ke_fom,
        ke_rom,
        en_fom,
        en_rom,
        final_vort_fom,
        final_vort_rom,
        final_vorticity_rel_error = norm(final_vort_rom - final_vort_fom) / max(norm(final_vort_fom), eps(Float64)),
        near_wall_vorticity_rel_error = norm(near_rom - near_fom) / max(norm(near_fom), eps(Float64)),
        near_wall_vorticity_corr = spatial_correlation(near_rom, near_fom),
        near_wall_enstrophy_ratio = near_en_rom / max(near_en_fom, eps(Float64)),
        ke_rmse = rel_series_error(ke_rom, ke_fom),
        en_rmse = rel_series_error(en_rom, en_fom),
        p1_fom,
        p2_fom,
        p1_rom,
        p2_rom,
        energy_retained,
    )
end

function method_rows(datas)
    rows = Any[]
    for item in METHOD_RUNS
        d = datas[item.label]
        m = d.meta
        push!(rows, [
            d.label,
            d.run_id,
            m.rfit,
            @sprintf("%.5f", d.energy_retained),
            hasproperty(m, :runtime_seconds) ? @sprintf("%.3f", m.runtime_seconds) : "NaN",
            @sprintf("%.5f", mean(d.rel_error)),
            @sprintf("%.5f", maximum(d.rel_error)),
            @sprintf("%.5f", last(d.rel_error)),
            @sprintf("%.5f", d.ke_rmse),
            @sprintf("%.5f", d.en_rmse),
            @sprintf("%.5f", abs(d.p1_rom.time - d.p1_fom.time)),
            @sprintf("%.5f", abs(d.p2_rom.time - d.p2_fom.time)),
            @sprintf("%.5f", abs(d.p1_rom.value - d.p1_fom.value) / max(abs(d.p1_fom.value), eps(Float64))),
            @sprintf("%.5f", abs(d.p2_rom.value - d.p2_fom.value) / max(abs(d.p2_fom.value), eps(Float64))),
            @sprintf("%.5f", d.final_vorticity_rel_error),
            @sprintf("%.5f", d.near_wall_vorticity_rel_error),
            @sprintf("%.5f", d.near_wall_vorticity_corr),
            @sprintf("%.5f", d.near_wall_enstrophy_ratio),
        ])
    end
    return rows
end

function operator_comparison_row(exact_id, intrusive_id)
    ef = load_rom_fit(exact_id)
    inf = load_rom_fit(intrusive_id)
    return [
        size(ef.L, 1),
        @sprintf("%.3e", norm(ef.L - inf.L)),
        @sprintf("%.3e", norm(ef.L - inf.L) / max(norm(inf.L), eps(Float64))),
        @sprintf("%.3e", norm(vec(ef.Q) - vec(inf.Q))),
        @sprintf("%.3e", norm(vec(ef.Q) - vec(inf.Q)) / max(norm(vec(inf.Q)), eps(Float64))),
    ]
end

function plot_method_comparison(datas)
    figs = figure_dir()
    colors = Makie.wong_colors()
    labels = [item.label for item in METHOD_RUNS]
    styles = Dict(
        labels[1] => (color = colors[4], linestyle = :dash, linewidth = THESIS_LINEWIDTH),
        labels[2] => (color = colors[1], linestyle = :solid, linewidth = THESIS_LINEWIDTH),
        labels[3] => (color = colors[2], linestyle = :dash, linewidth = THESIS_LINEWIDTH),
        labels[4] => (color = colors[3], linestyle = :solid, linewidth = THESIS_LINEWIDTH),
    )

    fig = Figure(size = (1150, 720))
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "relative state error", yscale = log10)
    for label in labels
        d = datas[label]
        lines!(ax, d.times, d.rel_error; label, styles[label]...)
    end
    axislegend(ax; position = :lt)
    save(joinpath(figs, "dipole_rom_relative_error.png"), fig)

    fig = Figure(size = (1150, 900))
    ax1 = Axis(fig[1, 1]; xlabel = "time", ylabel = "kinetic energy")
    ax2 = Axis(fig[2, 1]; xlabel = "time", ylabel = "enstrophy")
    ref = datas[labels[2]]
    lines!(ax1, ref.times, ref.ke_fom; label = "FOM", color = :black, linewidth = THESIS_REFERENCE_LINEWIDTH)
    lines!(ax2, ref.times, ref.en_fom; label = "FOM", color = :black, linewidth = THESIS_REFERENCE_LINEWIDTH)
    for label in labels
        d = datas[label]
        lines!(ax1, d.times, d.ke_rom; label, styles[label]...)
        lines!(ax2, d.times, d.en_rom; label, styles[label]...)
    end
    axislegend(ax1; position = :rt)
    axislegend(ax2; position = :rt)
    save(joinpath(figs, "dipole_rom_diagnostics.png"), fig)
    save(joinpath(figs, "dipole_rom_energy_enstrophy.png"), fig)
end

function plot_final_vorticity(datas, setup)
    figs = figure_dir()
    coords = getindex.(Array.(setup.xp), setup.Ip.indices)
    labels = [
        "FOM",
        "Exact OpInf, eta=0.99",
        "Exact OpInf, eta=0.95",
        "Trajectory OpInf",
    ]
    fields = [
        datas["Exact OpInf, eta=0.99"].final_vort_fom,
        datas["Exact OpInf, eta=0.99"].final_vort_rom,
        datas["Exact OpInf, eta=0.95"].final_vort_rom,
        datas["Trajectory OpInf, eta=0.99 request"].final_vort_rom,
    ]
    matrices = [reshape(f, length(coords[1]), length(coords[2])) for f in fields]
    clim = maximum(quantile(abs.(vec(m)), 0.99) for m in matrices)
    clim = max(clim, eps(Float64))

    fig = Figure(size = (1450, 900), figure_padding = 18)
    hm = nothing
    for i in 1:4
        row = i <= 2 ? 1 : 2
        col = i <= 2 ? i : i - 2
        ax = Axis(fig[row, col]; title = labels[i], xlabel = "x", ylabel = "y", aspect = DataAspect())
        plot = heatmap!(ax, coords..., matrices[i]; colormap = :RdBu, colorrange = (-clim, clim))
        hm = isnothing(hm) ? plot : hm
        xlims!(ax, 0.4, 1.0)
        ylims!(ax, 0.0, 0.6)
    end
    Colorbar(fig[:, 3], hm; label = "vorticity")
    save(joinpath(figs, "dipole_rom_final_vorticity_zoom.png"), fig)
end

function write_plot_data(datas, setup)
    out = analysis_dir()
    labels = [item.label for item in METHOD_RUNS]
    ref = datas["Exact OpInf, eta=0.99"]
    series = (;
        labels,
        times = ref.times,
        ke_fom = ref.ke_fom,
        en_fom = ref.en_fom,
        rel_error = Dict(label => datas[label].rel_error for label in labels),
        ke_rom = Dict(label => datas[label].ke_rom for label in labels),
        en_rom = Dict(label => datas[label].en_rom for label in labels),
    )

    coords_full = getindex.(Array.(setup.xp), setup.Ip.indices)
    xfull = coords_full[1]
    yfull = coords_full[2]
    ix = findall(x -> 0.4 <= x <= 1.0, xfull)
    iy = findall(y -> 0.0 <= y <= 0.6, yfull)
    ix = ix[1:4:end]
    iy = iy[1:4:end]

    fields = Dict(
        "FOM" => datas["Exact OpInf, eta=0.99"].final_vort_fom,
        "Exact OpInf, eta=0.99" => datas["Exact OpInf, eta=0.99"].final_vort_rom,
        "Exact OpInf, eta=0.95" => datas["Exact OpInf, eta=0.95"].final_vort_rom,
        "Trajectory OpInf" => datas["Trajectory OpInf, eta=0.99 request"].final_vort_rom,
    )
    nx = length(xfull)
    ny = length(yfull)
    vorticity_zoom = Dict{String,Matrix{Float64}}()
    for (label, vecfield) in fields
        mat = reshape(vecfield, nx, ny)
        vorticity_zoom[label] = Matrix{Float64}(mat[ix, iy])
    end
    vort = (;
        labels = ["FOM", "Exact OpInf, eta=0.99", "Exact OpInf, eta=0.95", "Trajectory OpInf"],
        x = xfull[ix],
        y = yfull[iy],
        fields = vorticity_zoom,
    )

    save_object(joinpath(out, "plot_data.jld2"), (; series, vort))
end

function main()
    out = analysis_dir()
    cfg = load_fom_config_from_run(FOM_ID)
    snapshots = load_fom_snapshots(FOM_ID)
    setup = getsetup(cfg; backend = KernelAbstractions.CPU())
    psolver = NS.default_psolver(setup)
    u_fom = initial_velocity(cfg, setup, psolver)
    u_rom = initial_velocity(cfg, setup, psolver)
    mask = near_wall_mask(setup)
    weights = velocity_mass_vector(setup, u_fom)
    spectrum = pod_spectrum_from_weighted_gram(snapshots.u, weights)
    energy_by_rank = Dict(
        r => _pod_energy_retained(spectrum.S, r)
        for r in 1:length(spectrum.S)
    )

    datas = Dict{String,Any}()
    for item in METHOD_RUNS
        @info "Computing dipole ROM diagnostics" item.label item.id
        datas[item.label] = summarize_run(item.label, item.id, snapshots, setup, u_fom, u_rom, mask, energy_by_rank)
    end

    write_csv(
        joinpath(out, "method_comparison.csv"),
        [
            "label", "run_id", "rank", "energy_retained", "runtime_seconds",
            "mean_rel_state_error", "max_rel_state_error", "final_rel_state_error",
            "rel_ke_rmse", "rel_enstrophy_rmse",
            "first_peak_time_abs_error", "second_peak_time_abs_error",
            "first_peak_value_rel_error", "second_peak_value_rel_error",
            "final_vorticity_rel_error", "near_wall_vorticity_rel_error",
            "near_wall_vorticity_correlation", "near_wall_enstrophy_ratio",
        ],
        method_rows(datas),
    )

    write_csv(
        joinpath(out, "operator_exact_vs_intrusive.csv"),
        ["rank", "L_abs_fro", "L_rel_fro", "Q_abs_fro", "Q_rel_fro"],
        [operator_comparison_row(
            "rom_from_e592f3a8_exactopinf_bglob_e0p99_87d75a35",
            "rom_from_e592f3a8_intrusive_bglob_e0p99_tagdirect_intrusive_e12a4370",
        )],
    )

    write_plot_data(datas, setup)

    println("Wrote analysis tables to: ", out)
    println("Wrote compact plot data to: ", joinpath(out, "plot_data.jld2"))
end

main()
