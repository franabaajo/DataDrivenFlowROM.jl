include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Config.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Paths.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Registry.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Setup.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "FOM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "POD.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "ROM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Plotting.jl"))

function list_from_toml(d::AbstractDict, key::AbstractString, default)
    if !haskey(d, key)
        return default
    end
    v = d[key]
    return v isa AbstractVector ? v : [v]
end

function stretch_from_toml_list(vals)
    return [_stretch_from_toml(v) for v in vals]
end

if length(ARGS) < 2
    error("Usage: julia --project=. experiments/DipoleWallCollision/scripts/sweep_dipole_fom.jl <grid_study.toml> <fom.toml>")
end

sweep = TOML.parsefile(ARGS[1])
basecfg = load_fom_config(ARGS[2])

n_list = Int.(list_from_toml(sweep, "n", [basecfg.n]))
cfl_list = Float64.(list_from_toml(sweep, "cfl", [basecfg.cfl]))
stretch_list = stretch_from_toml_list(list_from_toml(sweep, "stretch", [basecfg.stretch]))
nsnapshot_list = Int.(list_from_toml(sweep, "nsnapshot", [basecfg.nsnapshot]))
viscosity_list = Float64.(list_from_toml(sweep, "viscosity", [basecfg.viscosity]))

for n in n_list, cfl in cfl_list, stretch in stretch_list, ns in nsnapshot_list, nu in viscosity_list
    cfg = FOMConfig(
        viscosity = nu,
        r0 = basecfg.r0,
        p1 = basecfg.p1,
        p2 = basecfg.p2,
        initialenergy = basecfg.initialenergy,
        n = n,
        stretch = stretch,
        lims = basecfg.lims,
        tsim = basecfg.tsim,
        cfl = cfl,
        nsnapshot = ns,
    )

    run_id = run_fom(cfg; skip_if_exists = true)
    println("Finished: ", run_id)
end
