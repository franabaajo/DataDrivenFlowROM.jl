include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Config.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Paths.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Registry.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Setup.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "FOM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "POD.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "ROM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Plotting.jl"))

function main()
    if length(ARGS) < 2
        error("Usage: julia --project=. experiments/RayleighBenardConvection/scripts/RayleighBenardPOD.jl <fom_config.toml> <rom_config.toml>")
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
