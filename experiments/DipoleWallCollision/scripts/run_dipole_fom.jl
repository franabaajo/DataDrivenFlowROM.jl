include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Config.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Paths.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Registry.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Setup.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "FOM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "POD.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "ROM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "dipole", "Plotting.jl"))

if length(ARGS) < 1
    error(
        "Usage: ./dipole.sh run-fom <fom_config.toml> " *
        "(or: julia --project=. experiments/DipoleWallCollision/scripts/run_dipole_fom.jl <fom_config.toml>)",
    )
end
cfg = load_fom_config(ARGS[1])
run_id = run_fom(cfg)
plot_fom_outputs(run_id; make_movie = true)

println("FOM run id: ", run_id)
println("Saved in: ", fom_dir(run_id))
