include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Config.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Paths.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Registry.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Setup.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "FOM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "POD.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "ROM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Plotting.jl"))

if length(ARGS) < 2
    error(
        "Usage: ./dipole.sh run-rom <fom_run_id> <rom_config.toml> " *
        "(or: julia --project=. experiments/DipoleWallCollision/scripts/run_dipole_rom.jl <fom_run_id> <rom_config.toml>)",
    )
end

parent_fom_id = ARGS[1]
romcfg = load_rom_config(ARGS[2])

run_id = run_rom(parent_fom_id, romcfg)
plot_rom_outputs(run_id; make_movie = true)

println("ROM run id: ", run_id)
println("Saved in: ", rom_dir(run_id))
