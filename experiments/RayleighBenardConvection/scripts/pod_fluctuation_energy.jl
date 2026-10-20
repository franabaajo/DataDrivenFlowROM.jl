using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Config.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Paths.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "Setup.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "rayleighbenard", "FOM.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "src", "workflows", "POD.jl"))

function centered_weighted_spectrum(
    X::AbstractMatrix,
    w::AbstractVector;
    chunk_rows::Integer = 25_000,
)
    N, K = size(X)
    length(w) == N || error("Weight length mismatch: length(w) = $(length(w)) but size(X, 1) = $N")
    chunk_rows >= 1 || error("chunk_rows must be >= 1, got $chunk_rows")

    T = promote_type(eltype(X), eltype(w))
    gram = zeros(T, K, K)

    row = 1
    while row <= N
        lastrow = min(row + Int(chunk_rows) - 1, N)
        Xc = Matrix{T}(@view X[row:lastrow, :])
        Xc .-= sum(Xc; dims = 2) ./ K
        Y = sqrt.(T.(@view w[row:lastrow])) .* Xc
        mul!(gram, transpose(Y), Y, one(T), one(T))
        row = lastrow + 1
    end

    eig = eigen(Symmetric(gram))
    order = sortperm(eig.values; rev = true)
    values = max.(eig.values[order], zero(eltype(eig.values)))
    return sqrt.(values)
end

function cumulative_energy_safe(S::AbstractVector)
    e = S .^ 2
    total = sum(e)
    total <= 0 && return zeros(eltype(e), length(e))
    return cumsum(e) ./ total
end

function rank_for_energy(S::AbstractVector, eta::Real)
    ce = cumulative_energy_safe(S)
    idx = findfirst(>=(eta), ce)
    return idx === nothing ? length(S) : idx
end

function main()
    if isempty(ARGS)
        error("Usage: julia --project=. experiments/RayleighBenardConvection/scripts/pod_fluctuation_energy.jl <fom_run_id> [max_r]")
    end

    run_id = ARGS[1]
    max_r = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 30
    max_r >= 1 || error("max_r must be >= 1")

    cfg = load_fom_config_from_run(run_id)
    snapshots = load_fom_snapshots(run_id)
    setup = getsetup(cfg; backend = KernelAbstractions.CPU())
    w = state_mass_vector(setup, NS.vectorfield(setup), NS.scalarfield(setup))

    raw_spectrum = pod_spectrum_from_weighted_gram(snapshots.state, w)
    Sraw = raw_spectrum.S
    Sfluct = centered_weighted_spectrum(snapshots.state, w)

    raw_energy = cumulative_energy_safe(Sraw)
    fluct_energy = cumulative_energy_safe(Sfluct)
    rmax = min(max_r, length(raw_energy), length(fluct_energy))

    println("FOM run id: ", run_id)
    println("Snapshots: ", size(snapshots.state, 2))
    println("State dofs: ", size(snapshots.state, 1))
    println()
    println(" rank    raw_energy    fluctuation_energy")
    for r in 1:rmax
        @printf("%5d  %12.8f  %18.8f\n", r, raw_energy[r], fluct_energy[r])
    end

    println()
    println("Ranks needed by centered fluctuation energy:")
    for eta in (0.99, 0.999, 0.9999)
        @printf("  %.4f -> r = %d\n", eta, rank_for_energy(Sfluct, eta))
    end
end

main()
