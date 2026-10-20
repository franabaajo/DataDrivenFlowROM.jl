@testitem "Space-local POD: W-orthonormal and Single-Block Support (S1)" begin
    using LinearAlgebra
    using SparseArrays

    include(joinpath(@__DIR__, "..", "src", "workflows", "POD.jl"))

    N, K = 60, 15
    X = randn(N, K)
    w = rand(N) .+ 0.1 # strictly positive weights

    blocks = [collect(1:20), collect(21:45), collect(46:N)]

    sl = build_space_local_basis(
        X,
        blocks;
        requested_total_rank = 9,
        row_weights = w,
        strategy = :S1,
        sparse_basis = true,
    )

    @test sl.rtot == 9
    @test sum(sl.ranks) == sl.rtot
    @test sl.rom_basis isa SparseMatrixCSC
    @test size(sl.rom_basis) == (N, sl.rtot)
    @test size(sl.projector) == (sl.rtot, N)

    # Check W-orthonormality: Phi' * W * Phi ≈ I.
    Phi = sl.rom_basis
    G = Matrix(Phi' * Diagonal(w) * Phi)
    @test norm(G - I(sl.rtot)) ≤ 1e-8

    # Each space-local basis vector must live in exactly one block.
    bid = zeros(Int, N)
    for (i, rows) in enumerate(blocks)
        for r in rows
            bid[r] = i
        end
    end
    for j in 1:size(Phi, 2)
        rows = rowvals(Phi)[nzrange(Phi, j)]
        @test length(unique(bid[rows])) == 1
    end
end

@testitem "Space-local POD: local_ranks Overrides Energy (S3)" begin
    using LinearAlgebra
    using SparseArrays

    include(joinpath(@__DIR__, "..", "src", "workflows", "POD.jl"))

    N, K = 40, 10
    X = randn(N, K)
    w = rand(N) .+ 0.1
    blocks = [collect(1:10), collect(11:N)]

    sl = build_space_local_basis(
        X,
        blocks;
        row_weights = w,
        local_ranks = [2, 3],
        energy_fraction = 0.9999, # should be ignored because local_ranks is explicit
        strategy = :S3,
        sparse_basis = false,
    )

    @test sl.ranks == [2, 3]
    @test sl.rtot == 5
    @test sl.rom_basis isa Matrix
    @test size(sl.rom_basis) == (N, sl.rtot)
    @test size(sl.projector) == (sl.rtot, N)
end

