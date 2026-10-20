@testitem "Public package API" begin
    using DataDrivenFlowROM
    using LinearAlgebra

    r = 3
    @test n_unique_quadratic(r) == 6
    @test n_features_affine_quadratic(r) == 10

    A = rank_ensuring_states_quadratic(r)
    @test size(A) == (r, r + n_unique_quadratic(r))

    features, pairs = quadratic_features_symmetric(A)
    @test size(features) == (size(A, 2), r + n_unique_quadratic(r))
    @test length(pairs) == n_unique_quadratic(r)
    @test rank(features) == size(features, 2)

    X = [1.0 2.0 3.0; 0.0 1.0 0.0; 2.0 0.0 1.0]
    spectrum = pod_spectrum_from_gram(X)
    @test length(spectrum.S) == size(X, 2)
    @test issorted(spectrum.S; rev = true)
end

