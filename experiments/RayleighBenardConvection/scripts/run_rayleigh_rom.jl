include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Config.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Paths.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Registry.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Setup.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "FOM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "POD.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "ROM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Plotting.jl"))

if length(ARGS) < 2
    error("Usage: julia --project=. experiments/RayleighBenardConvection/scripts/run_rayleigh_rom.jl <fom_run_id> <rom_config.toml>")
end

parent_fom_id = ARGS[1]
romcfg = load_rom_config(ARGS[2])

run_id = run_rom(parent_fom_id, romcfg)
plot_rom_outputs(run_id; make_movie = true)

println("ROM run id: ", run_id)
println("Saved in: ", rom_dir(run_id))
