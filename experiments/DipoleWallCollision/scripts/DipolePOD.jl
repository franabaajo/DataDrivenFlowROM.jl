include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Config.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Paths.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Registry.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Setup.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "FOM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "POD.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "ROM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Plotting.jl"))

function main()
    if length(ARGS) < 2
        error(
            "Usage: ./dipole.sh run-all <fom_config.toml> <rom_config.toml> " *
            "(or: julia --project=. experiments/DipoleWallCollision/scripts/DipolePOD.jl <fom_config.toml> <rom_config.toml>)",
        )
    end
    fomcfg = load_fom_config(ARGS[1])
    romcfg = load_rom_config(ARGS[2])

    fom_run_id = run_fom(fomcfg)
    plot_fom_outputs(fom_run_id; make_movie = true)

    rom_run_id = run_rom(fom_run_id, romcfg)
    plot_rom_outputs(rom_run_id; make_movie = true)

    println("FOM run id: ", fom_run_id)
    println("ROM run id: ", rom_run_id)
end

main()
