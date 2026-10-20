include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Config.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Paths.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Registry.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Setup.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "FOM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "POD.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "ROM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Plotting.jl"))

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

function limits_from_toml(d::AbstractDict, key::AbstractString, default)
    if !haskey(d, key)
        return default
    end
    vals = d[key]
    length(vals) == 2 || error("$key must contain exactly two values")
    return (Float64(vals[1]), Float64(vals[2]))
end

function grid_pairs_from_toml(sweep::AbstractDict, basecfg::FOMConfig)
    if haskey(sweep, "n")
        if haskey(sweep, "nx") || haskey(sweep, "ny")
            error("Use either square-grid key `n` or rectangular-grid keys `nx`/`ny`, not both.")
        end
        n_list = Int.(list_from_toml(sweep, "n", [basecfg.nx]))
        return [(n, n) for n in n_list]
    end

    nx_list = Int.(list_from_toml(sweep, "nx", [basecfg.nx]))
    ny_list = Int.(list_from_toml(sweep, "ny", [basecfg.ny]))
    return [(nx, ny) for nx in nx_list for ny in ny_list]
end

function stretch_pair_lists_from_toml(sweep::AbstractDict, basecfg::FOMConfig)
    if haskey(sweep, "stretch")
        if haskey(sweep, "stretch_x") || haskey(sweep, "stretch_y")
            error("Use either shared key `stretch` or directional keys `stretch_x`/`stretch_y`, not both.")
        end
        stretch_list = stretch_from_toml_list(list_from_toml(sweep, "stretch", [basecfg.stretch_x]))
        return stretch_list, stretch_list
    end

    stretch_x_list = stretch_from_toml_list(list_from_toml(sweep, "stretch_x", [basecfg.stretch_x]))
    stretch_y_list = stretch_from_toml_list(list_from_toml(sweep, "stretch_y", [basecfg.stretch_y]))
    return stretch_x_list, stretch_y_list
end

if length(ARGS) > 2
    error("Usage: julia --project=. experiments/RayleighBenardConvection/scripts/sweep_rayleigh_fom.jl [grid_study.toml] [fom.toml]")
end

grid_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "configs", "grid_study.toml")
fom_path = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "..", "configs", "fom_base.toml")

sweep = TOML.parsefile(grid_path)
basecfg = load_fom_config(fom_path)

grid_pairs = grid_pairs_from_toml(sweep, basecfg)
xlims = limits_from_toml(sweep, "xlims", basecfg.xlims)
ylims = limits_from_toml(sweep, "ylims", basecfg.ylims)
cfl_list = Float64.(list_from_toml(sweep, "cfl", [basecfg.cfl]))
stretch_x_list, stretch_y_list = stretch_pair_lists_from_toml(sweep, basecfg)
nsnapshot_list = Int.(list_from_toml(sweep, "nsnapshot", [basecfg.nsnapshot]))
viscosity_list = Float64.(list_from_toml(sweep, "viscosity", [basecfg.viscosity]))
conductivity_list = Float64.(list_from_toml(sweep, "conductivity", [basecfg.conductivity]))

for (nx, ny) in grid_pairs, cfl in cfl_list, stx in stretch_x_list, sty in stretch_y_list,
    ns in nsnapshot_list, nu in viscosity_list, kappa in conductivity_list

    cfg = FOMConfig(
        viscosity = nu,
        conductivity = kappa,
        gravity = basecfg.gravity,
        gdir = basecfg.gdir,
        dodissipation = basecfg.dodissipation,
        nx = nx,
        ny = ny,
        xlims = xlims,
        ylims = ylims,
        stretch_x = stx,
        stretch_y = sty,
        temp_hot = basecfg.temp_hot,
        temp_cold = basecfg.temp_cold,
        initial_temperature_profile = basecfg.initial_temperature_profile,
        perturbation_amplitude = basecfg.perturbation_amplitude,
        perturbation_wavenumber = basecfg.perturbation_wavenumber,
        tsim = basecfg.tsim,
        cfl = cfl,
        nsnapshot = ns,
    )

    run_id = run_fom(cfg; skip_if_exists = true)
    println("Finished: ", run_id)
end
