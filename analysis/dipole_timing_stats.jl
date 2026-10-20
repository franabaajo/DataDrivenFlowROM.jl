using Printf

include(joinpath(@__DIR__, "..", "src", "dipole", "Config.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "Paths.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "Registry.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "Setup.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "FOM.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "POD.jl"))
include(joinpath(@__DIR__, "..", "src", "dipole", "ROM.jl"))

const FOM_ID = "fom_n1536_cfl0p67_st1p2_ns50_nu0p0004_e592f3a8"
const RUNS = [
    (label = "Exact OpInf, eta=0.95", id = "rom_from_e592f3a8_exactopinf_bglob_e0p95_d9547ec6"),
    (label = "Exact OpInf, eta=0.99", id = "rom_from_e592f3a8_exactopinf_bglob_e0p99_87d75a35"),
    (label = "Intrusive Galerkin, eta=0.99", id = "rom_from_e592f3a8_intrusive_bglob_e0p99_tagdirect_intrusive_e12a4370"),
    (label = "Trajectory OpInf", id = "rom_from_e592f3a8_trajopinf_bglob_e0p99_lam0p01_tagtrajopinf_e0p99_2b364ac5"),
]

analysis_dir() = mkpath(joinpath(@__DIR__, "..", "results", "tables", "dipole"))

function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ","))
        for row in rows
            println(io, join(row, ","))
        end
    end
end

function timed_advance(fit, traj, cfg; repeats = 6)
    a0 = copy(traj.a_exact[:, 1])

    # Warm up compilation and allocations before measuring.
    simulate_rom_adaptive(
        a0,
        fit.L,
        fit.Q;
        tsim = cfg.tsim,
        nsave = cfg.nsnapshot,
        safety = 0.03,
        growthfac = 1.05,
    )

    times = Float64[]
    for _ in 1:repeats
        elapsed = @elapsed simulate_rom_adaptive(
            a0,
            fit.L,
            fit.Q;
            tsim = cfg.tsim,
            nsave = cfg.nsnapshot,
            safety = 0.03,
            growthfac = 1.05,
        )
        push!(times, elapsed)
    end
    sort!(times)
    return (; median = times[cld(length(times), 2)], minimum = first(times), maximum = last(times))
end

function main()
    cfg = load_fom_config_from_run(FOM_ID)
    fom_meta = load_fom_metadata(FOM_ID)
    fom_runtime = Float64(fom_meta.runtime_seconds)

    rows = Any[]
    for run in RUNS
        meta = load_rom_metadata(run.id)
        fit = load_rom_fit(run.id)
        traj = load_rom_trajectory(run.id)
        t = timed_advance(fit, traj, cfg)
        push!(rows, [
            run.label,
            run.id,
            string(meta.rfit),
            @sprintf("%.6f", t.median),
            @sprintf("%.6f", t.minimum),
            @sprintf("%.6f", t.maximum),
            @sprintf("%.2f", fom_runtime / t.median),
            @sprintf("%.2f", fom_runtime / t.minimum),
        ])
    end

    write_csv(
        joinpath(analysis_dir(), "timing_comparison.csv"),
        [
            "label", "run_id", "rank", "advance_median_seconds",
            "advance_min_seconds", "advance_max_seconds",
            "speedup_vs_fom_median", "speedup_vs_fom_best",
        ],
        rows,
    )
    println("FOM runtime_seconds = ", fom_runtime)
    println("Wrote: ", joinpath(analysis_dir(), "timing_comparison.csv"))
end

main()
