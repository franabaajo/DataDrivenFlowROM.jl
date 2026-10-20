using CairoMakie
using KernelAbstractions
using LinearAlgebra
using Printf
using Statistics
import IncompressibleNavierStokes as NS

include(joinpath(@__DIR__, "..", "src", "dipole", "Config.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "Paths.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "Registry.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "Setup.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "FOM.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "POD.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "ROM.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "Stats.jl"))

const FOM_ID = "fom_n1536_cfl0p67_st1p2_ns50_nu0p0004_e592f3a8"
const RUNS = [
    (label = "Exact OpInf, eta=0.95", short = "Exact OpInf, eta=0.95", id = "rom_from_e592f3a8_exactopinf_bglob_e0p95_d9547ec6"),
    (label = "Exact OpInf, eta=0.99", short = "Exact OpInf, eta=0.99", id = "rom_from_e592f3a8_exactopinf_bglob_e0p99_87d75a35"),
    (label = "Intrusive Galerkin, eta=0.99", short = "Intrusive, eta=0.99", id = "rom_from_e592f3a8_intrusive_bglob_e0p99_tagdirect_intrusive_e12a4370"),
    (label = "Trajectory OpInf, eta=0.99 request", short = "Trajectory OpInf", id = "rom_from_e592f3a8_trajopinf_bglob_e0p99_lam0p01_tagtrajopinf_e0p99_2b364ac5"),
]

figure_dir() = mkpath(joinpath(@__DIR__, "..", "results", "figures"))
analysis_dir() = mkpath(joinpath(@__DIR__, "..", "results", "tables", "dipole"))

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

function spatial_correlation(a, b)
    aa = a .- mean(a)
    bb = b .- mean(b)
    return dot(aa, bb) / max(norm(aa) * norm(bb), eps(Float64))
end

function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join(row, ","))
        end
    end
end

function main()
    cfg = load_fom_config_from_run(FOM_ID)
    snapshots = load_fom_snapshots(FOM_ID)
    setup = getsetup(cfg; backend = KernelAbstractions.CPU())
    psolver = NS.default_psolver(setup)
    u = initial_velocity(cfg, setup, psolver)
    u_fom = initial_velocity(cfg, setup, psolver)
    u_rom = initial_velocity(cfg, setup, psolver)
    weights = velocity_mass_vector(setup, u)
    spectrum = pod_spectrum_from_weighted_gram(snapshots.u, weights)
    rank_energy = Dict(r => _pod_energy_retained(spectrum.S, r) for r in 1:length(spectrum.S))

    times = collect(range(0.0, cfg.tsim; length = size(snapshots.u, 2)))
    errs = Dict{String,Vector{Float64}}()
    final_vorts = Dict{String,Matrix{Float64}}()
    diagnostics = Dict{String,Any}()
    rows = Any[]
    diag_rows = Any[]

    coords_full = getindex.(Array.(setup.xp), setup.Ip.indices)
    xfull = coords_full[1]
    yfull = coords_full[2]
    ix = findall(x -> 0.4 <= x <= 1.0, xfull)[1:4:end]
    iy = findall(y -> 0.0 <= y <= 0.6, yfull)[1:4:end]
    xzoom = xfull[ix]
    yzoom = yfull[iy]
    near_mask = vec([0.4 <= xi <= 1.0 && 0.0 <= yj <= 0.6 for xi in xfull, yj in yfull])

    copyto!(vec(u), @view(snapshots.u[:, end]))
    NS.apply_bc_u!(u, cfg.tsim, setup)
    vort_fom = (NS.vorticity(u, setup) |> Array)[setup.Ip]
    final_vorts["FOM"] = Matrix(vort_fom[ix, iy])

    for item in RUNS
        @info "Processing ROM" item.label item.id
        meta = load_rom_metadata(item.id)
        fit = load_rom_fit(item.id)
        traj = load_rom_trajectory(item.id)
        K = min(size(traj.a_rom, 2), size(snapshots.u, 2))
        e = zeros(Float64, K)
        ke_fom = zeros(Float64, K)
        ke_rom = zeros(Float64, K)
        en_fom = zeros(Float64, K)
        en_rom = zeros(Float64, K)
        urom_vec = zeros(eltype(fit.rom_basis), size(fit.rom_basis, 1))
        for k in 1:K
            t = times[k]
            copyto!(vec(u_fom), @view(snapshots.u[:, k]))
            mul!(urom_vec, fit.rom_basis, @view(traj.a_rom[:, k]))
            copyto!(vec(u_rom), urom_vec)

            e[k] = rel_state_error(urom_vec, @view(snapshots.u[:, k]))
            NS.apply_bc_u!(u_fom, t, setup)
            NS.apply_bc_u!(u_rom, t, setup)
            ke_fom[k] = NS.total_kinetic_energy(u_fom, setup)
            ke_rom[k] = NS.total_kinetic_energy(u_rom, setup)
            en_fom[k] = total_enstrophy(u_fom, setup)
            en_rom[k] = total_enstrophy(u_rom, setup)
        end
        errs[item.label] = e
        diagnostics[item.label] = (; times = times[1:K], ke_fom, ke_rom, en_fom, en_rom)

        mul!(urom_vec, fit.rom_basis, @view(traj.a_rom[:, K]))
        copyto!(vec(u), urom_vec)
        NS.apply_bc_u!(u, times[K], setup)
        vort = (NS.vorticity(u, setup) |> Array)[setup.Ip]
        final_vorts[item.short] = Matrix(vort[ix, iy])

        p1_fom = peak_in_window(times[1:K], en_fom, 0.20, 0.45)
        p2_fom = peak_in_window(times[1:K], en_fom, 0.45, 0.80)
        p1_rom = peak_in_window(times[1:K], en_rom, 0.20, 0.45)
        p2_rom = peak_in_window(times[1:K], en_rom, 0.45, 0.80)
        final_fom_vec = vec(vort_fom)
        final_rom_vec = vec((NS.vorticity(u, setup) |> Array)[setup.Ip])
        near_fom = final_fom_vec[near_mask]
        near_rom = final_rom_vec[near_mask]
        near_en_fom = 0.5 * mean(abs2, near_fom)
        near_en_rom = 0.5 * mean(abs2, near_rom)

        push!(rows, [
            item.label,
            string(item.id),
            string(meta.rfit),
            @sprintf("%.8f", rank_energy[Int(meta.rfit)]),
            @sprintf("%.5f", mean(e)),
            @sprintf("%.5f", maximum(e)),
            @sprintf("%.5f", last(e)),
        ])
        push!(diag_rows, [
            item.label,
            string(item.id),
            string(meta.rfit),
            @sprintf("%.8f", rank_energy[Int(meta.rfit)]),
            @sprintf("%.5f", rel_series_error(ke_rom, ke_fom)),
            @sprintf("%.5f", rel_series_error(en_rom, en_fom)),
            @sprintf("%.5f", abs(p1_rom.time - p1_fom.time)),
            @sprintf("%.5f", abs(p2_rom.time - p2_fom.time)),
            @sprintf("%.5f", abs(p1_rom.value - p1_fom.value) / max(abs(p1_fom.value), eps(Float64))),
            @sprintf("%.5f", abs(p2_rom.value - p2_fom.value) / max(abs(p2_fom.value), eps(Float64))),
            @sprintf("%.5f", norm(final_rom_vec - final_fom_vec) / max(norm(final_fom_vec), eps(Float64))),
            @sprintf("%.5f", norm(near_rom - near_fom) / max(norm(near_fom), eps(Float64))),
            @sprintf("%.5f", spatial_correlation(near_rom, near_fom)),
            @sprintf("%.5f", near_en_rom / max(near_en_fom, eps(Float64))),
        ])
    end

    write_csv(
        joinpath(analysis_dir(), "light_method_errors.csv"),
        ["label", "run_id", "rank", "actual_energy_retained", "mean_rel_state_error", "max_rel_state_error", "final_rel_state_error"],
        rows,
    )
    write_csv(
        joinpath(analysis_dir(), "diagnostic_method_errors.csv"),
        [
            "label", "run_id", "rank", "actual_energy_retained",
            "rel_ke_rmse", "rel_enstrophy_rmse",
            "first_peak_time_abs_error", "second_peak_time_abs_error",
            "first_peak_value_rel_error", "second_peak_value_rel_error",
            "final_vorticity_rel_error", "near_wall_vorticity_rel_error",
            "near_wall_vorticity_correlation", "near_wall_enstrophy_ratio",
        ],
        diag_rows,
    )

    colors = Makie.wong_colors()
    styles = Dict(
        "Exact OpInf, eta=0.95" => (color = colors[4], linestyle = :dash, linewidth = 2.0),
        "Exact OpInf, eta=0.99" => (color = colors[1], linestyle = :solid, linewidth = 2.8),
        "Intrusive Galerkin, eta=0.99" => (color = colors[2], linestyle = :dash, linewidth = 1.8),
        "Trajectory OpInf, eta=0.99 request" => (color = colors[3], linestyle = :solid, linewidth = 2.4),
    )
    fig = Figure(size = (950, 560))
    ax = Axis(fig[1, 1]; xlabel = "time", ylabel = "relative state error", yscale = log10)
    for item in RUNS
        lines!(ax, times[1:length(errs[item.label])], errs[item.label]; label = item.short, styles[item.label]...)
    end
    axislegend(ax; position = :lt)
    save(joinpath(figure_dir(), "dipole_rom_relative_error.png"), fig)

    fig = Figure(size = (980, 720))
    ax1 = Axis(fig[1, 1]; xlabel = "time", ylabel = "kinetic energy")
    ax2 = Axis(fig[2, 1]; xlabel = "time", ylabel = "enstrophy")
    ref = diagnostics["Exact OpInf, eta=0.99"]
    lines!(ax1, ref.times, ref.ke_fom; label = "FOM", color = :black, linewidth = 2.4)
    lines!(ax2, ref.times, ref.en_fom; label = "FOM", color = :black, linewidth = 2.4)
    for item in RUNS
        d = diagnostics[item.label]
        lines!(ax1, d.times, d.ke_rom; label = item.short, styles[item.label]...)
        lines!(ax2, d.times, d.en_rom; label = item.short, styles[item.label]...)
    end
    axislegend(ax1; position = :rt)
    axislegend(ax2; position = :rt)
    save(joinpath(figure_dir(), "dipole_rom_energy_enstrophy.png"), fig)

    labels = ["FOM", "Exact OpInf, eta=0.99", "Exact OpInf, eta=0.95", "Trajectory OpInf"]
    clim = maximum(quantile(abs.(vec(final_vorts[label])), 0.99) for label in labels)
    clim = max(clim, eps(Float64))

    fig = Figure(size = (1450, 900), figure_padding = 18)
    hm = nothing
    for (i, label) in enumerate(labels)
        row = i <= 2 ? 1 : 2
        col = i <= 2 ? i : i - 2
        ax = Axis(fig[row, col]; title = label, xlabel = "x", ylabel = "y", aspect = DataAspect())
        plot = heatmap!(ax, xzoom, yzoom, final_vorts[label]; colormap = :RdBu, colorrange = (-clim, clim))
        hm = isnothing(hm) ? plot : hm
    end
    Colorbar(fig[:, 3], hm; label = "vorticity")
    save(joinpath(figure_dir(), "dipole_rom_final_vorticity_zoom.png"), fig)

    println("Wrote light dipole table and figures.")
end

main()
