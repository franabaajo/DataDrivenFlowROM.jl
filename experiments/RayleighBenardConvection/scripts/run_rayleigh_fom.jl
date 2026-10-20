include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Config.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Paths.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Registry.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Setup.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "FOM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "POD.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "ROM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Plotting.jl"))

if length(ARGS) < 1
    error("Usage: julia --project=. experiments/RayleighBenardConvection/scripts/run_rayleigh_fom.jl <fom_config.toml>")
end
cfg = load_fom_config(ARGS[1])
run_id = run_fom(cfg)
plot_fom_outputs(run_id; make_movie = true)

println("FOM run id: ", run_id)
println("Saved in: ", fom_dir(run_id))
