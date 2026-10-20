if !isdefined(@__MODULE__, :_WORKFLOW_POD_LOADED)
const _WORKFLOW_POD_LOADED = true

using LinearAlgebra
using SparseArrays

decompose_snapshots(snapshots) = svd(snapshots.u)

"""
Compute POD singular values and right singular vectors from the snapshot Gram matrix.

This avoids materializing the full left singular vector matrix `U`, which is expensive
for large state dimensions.
"""
function pod_spectrum_from_gram(X::AbstractMatrix)
    gram = Symmetric(Matrix(X' * X))
    eig = eigen(gram)
    order = sortperm(eig.values; rev = true)
    values = eig.values[order]
    values = max.(values, zero(eltype(values)))
    S = sqrt.(values)
    V = eig.vectors[:, order]
    return (; S, V)
end

"""
Compute POD singular values and right singular vectors from the *weighted* snapshot
Gram matrix.

Given snapshots `X` of size `(N, K)` and nonnegative row weights `w` of length `N`,
this constructs the Gram matrix

`G = X' * Diagonal(w) * X`,

which corresponds to POD in the weighted inner product `⟨x, y⟩_w = x' * Diagonal(w) * y`.

The computation is performed in row chunks to avoid materializing `Diagonal(sqrt.(w)) * X`
when `N` is large.
"""
function pod_spectrum_from_weighted_gram(
    X::AbstractMatrix,
    w::AbstractVector;
    chunk_rows::Integer = 25_000,
)
    N, K = size(X)
    length(w) == N || error("Weight length mismatch: length(w) = $(length(w)) but size(X,1) = $N")
    chunk_rows >= 1 || error("chunk_rows must be >= 1, got $chunk_rows")

    T = promote_type(eltype(X), eltype(w))
    gram = zeros(T, K, K)

    i = 1
    while i <= N
        j = min(i + Int(chunk_rows) - 1, N)
        Xc = @view X[i:j, :]
        wc = @view w[i:j]

        sqrtwc = sqrt.(T.(wc))
        Y = sqrtwc .* Xc
        mul!(gram, Y', Y, one(T), one(T))

        i = j + 1
    end

    gram = Symmetric(gram)
    eig = eigen(gram)
    order = sortperm(eig.values; rev = true)
    values = eig.values[order]
    values = max.(values, zero(eltype(values)))
    S = sqrt.(values)
    V = eig.vectors[:, order]
    return (; S, V)
end

"""
Recover one POD mode vector from snapshots and Gram decomposition output.
"""
function pod_mode_from_gram(X::AbstractMatrix, spectrum, imode::Integer)
    1 <= imode <= length(spectrum.S) || error("POD mode index out of bounds: $imode")
    sigma = spectrum.S[imode]
    sigma_max = isempty(spectrum.S) ? zero(sigma) : maximum(spectrum.S)
    tol = sqrt(eps(Float64)) * max(Float64(sigma_max), 1.0)
    if Float64(sigma) <= tol
        return zeros(eltype(X), size(X, 1))
    end
    mode = X * @view(spectrum.V[:, imode])
    mode ./= sigma
    return mode
end

function cumulative_energy(S)
    e = S .^ 2
    return cumsum(e) / sum(e)
end

function rank_from_energy(S; energy_fraction::Float64 = 0.99)
    return findfirst(>=(energy_fraction), cumulative_energy(S))
end

function choose_rank(S; requested::Union{Nothing,Int} = nothing, energy_fraction::Float64 = 0.99)
    r_energy = rank_from_energy(S; energy_fraction)
    return isnothing(requested) ? r_energy : min(requested, r_energy)
end

pod_basis(decomp, r::Int) = decomp.U[:, 1:r]

project_coordinates(basis, X) = basis' * X
project_derivatives(basis, dX) = basis' * dX

"""
Select rank based on cumulative energy, robust to the zero-energy case.
Returns `0` if `S` contains no energy.
"""
function rank_from_energy_safe(S; energy_fraction::Float64 = 0.99)
    isempty(S) && return 0
    e = S .^ 2
    tot = sum(e)
    tot <= 0 && return 0
    idx = findfirst(>=(energy_fraction), cumsum(e) / tot)
    return idx === nothing ? length(S) : idx
end

"""
Partition coordinate values into bins defined by sorted `cuts`.

Returns an integer vector `bin` with values in `1:(length(cuts)+1)`, where
`bin[i] = k` means:
- `k = 1`       if `coords[i] <= cuts[1]`
- `1 < k < end` if `cuts[k-1] < coords[i] <= cuts[k]`
- `k = end`     if `coords[i] > cuts[end]`
"""
function bins_from_cuts(coords::AbstractVector, cuts::AbstractVector)
    cuts_sorted = sort(collect(cuts))
    nb = length(cuts_sorted) + 1
    bin = Vector{Int}(undef, length(coords))
    @inbounds for i in eachindex(coords)
        bin[i] = searchsortedfirst(cuts_sorted, coords[i])
    end
    return bin, cuts_sorted, nb
end

"""
Compute interior cut locations that split the coordinate range of `coords` into
`nblocks` uniform bins.

Returns a sorted vector of length `nblocks-1` (or empty if `nblocks == 1`).
"""
function uniform_cuts(coords::AbstractVector, nblocks::Integer)
    nblocks >= 1 || error("nblocks must be >= 1, got $nblocks")
    nblocks == 1 && return Float64[]
    xmin, xmax = extrema(coords)
    xmin < xmax || return Float64[]
    edges = range(Float64(xmin), Float64(xmax); length = Int(nblocks) + 1)
    return collect(edges[2:end-1])
end

"""
Build non-overlapping row-index blocks from x/y coordinates.

If `ycuts` is `nothing`, returns `nbx` blocks based on x only.
If `ycuts` is provided, returns `nbx*nby` blocks based on a tensor-product grid.
"""
function blocks_from_xycoords(
    xcoords::AbstractVector,
    ycoords::Union{Nothing,AbstractVector};
    xcuts::AbstractVector,
    ycuts::Union{Nothing,AbstractVector} = nothing,
)
    length(xcoords) >= 1 || error("Empty coordinate vector")
    ycoords === nothing && ycuts !== nothing && error("ycuts provided but ycoords is nothing")
    if ycoords !== nothing
        length(ycoords) == length(xcoords) || error("x/y coordinate size mismatch")
    end

    bx, xcuts_sorted, nbx = bins_from_cuts(xcoords, xcuts)
    if ycuts === nothing
        blocks = [Int[] for _ in 1:nbx]
        @inbounds for i in eachindex(bx)
            push!(blocks[bx[i]], i)
        end
        return blocks, (; xcuts = xcuts_sorted, ycuts = nothing, nbx, nby = 1)
    end

    by, ycuts_sorted, nby = bins_from_cuts(ycoords, ycuts)
    blocks = [Int[] for _ in 1:(nbx * nby)]
    @inbounds for i in eachindex(bx)
        id = (by[i] - 1) * nbx + bx[i]
        push!(blocks[id], i)
    end
    return blocks, (; xcuts = xcuts_sorted, ycuts = ycuts_sorted, nbx, nby)
end

"""
Compute x/y coordinate vectors for a flattened velocity field `vec(u)` where
`u = vectorfield(setup)` has size `(setup.N..., setup.dimension())`.
"""
function velocity_xycoords(setup; T::Type = Float64)
    D = setup.dimension()
    D in (2, 3) || error("velocity_xycoords only supports D=2 or D=3, got D=$D")
    N = setup.N
    dims = (N..., D)
    LI = LinearIndices(ntuple(d -> Base.OneTo(dims[d]), length(dims)))

    xcoords = Vector{T}(undef, prod(dims))
    ycoords = Vector{T}(undef, prod(dims))

    if D == 2
        N1, N2 = N
        @inbounds for α in 1:D
            xα = Array(setup.xu[α][1])
            yα = Array(setup.xu[α][2])
            for j in 1:N2, i in 1:N1
                k = LI[i, j, α]
                xcoords[k] = xα[i]
                ycoords[k] = yα[j]
            end
        end
    else
        N1, N2, N3 = N
        # ycoords stores the second coordinate; z can be derived similarly if needed later.
        @inbounds for α in 1:D
            xα = Array(setup.xu[α][1])
            yα = Array(setup.xu[α][2])
            for k3 in 1:N3, j in 1:N2, i in 1:N1
                k = LI[i, j, k3, α]
                xcoords[k] = xα[i]
                ycoords[k] = yα[j]
            end
        end
    end

    return xcoords, ycoords
end

"""
Compute x/y coordinate vectors for a flattened scalar field `vec(p)` where
`p = scalarfield(setup)` has size `setup.N`.
"""
function scalar_xycoords(setup; T::Type = Float64)
    D = setup.dimension()
    D in (2, 3) || error("scalar_xycoords only supports D=2 or D=3, got D=$D")
    N = setup.N
    LI = LinearIndices(ntuple(d -> Base.OneTo(N[d]), D))

    xcoords = Vector{T}(undef, prod(N))
    ycoords = Vector{T}(undef, prod(N))

    if D == 2
        N1, N2 = N
        x = Array(setup.xp[1])
        y = Array(setup.xp[2])
        @inbounds for j in 1:N2, i in 1:N1
            k = LI[i, j]
            xcoords[k] = x[i]
            ycoords[k] = y[j]
        end
    else
        N1, N2, N3 = N
        x = Array(setup.xp[1])
        y = Array(setup.xp[2])
        @inbounds for k3 in 1:N3, j in 1:N2, i in 1:N1
            k = LI[i, j, k3]
            xcoords[k] = x[i]
            ycoords[k] = y[j]
        end
    end

    return xcoords, ycoords
end

"""
Build a space-local POD basis from snapshot matrix `X` by performing separate SVDs
on row blocks and assembling a block-sparse basis.

Arguments:
- `X`: `(N, K)` snapshots (columns are snapshots)
- `blocks`: vector of row-index vectors, non-overlapping and 1-based

Keyword arguments:
- `requested_total_rank`: for strategy `:S1`, the total number of modes to retain
- `energy_fraction`: for strategy `:S3`, per-block cumulative energy target
- `local_ranks`: optional explicit per-block ranks (overrides `energy_fraction`)
- `strategy`: `:S1` (global pooling of block singular values) or `:S3` (per-block energy)
- `sparse_basis`: whether to assemble `rom_basis` as `SparseMatrixCSC`

Returns a named tuple with `rom_basis`, `projector`, and diagnostics.
"""
function build_space_local_basis(
    X::AbstractMatrix,
    blocks::AbstractVector{<:AbstractVector{Int}};
    requested_total_rank::Union{Nothing,Int} = nothing,
    row_weights::Union{Nothing,AbstractVector} = nothing,
    energy_fraction::Float64 = 0.99,
    local_ranks::Union{Nothing,AbstractVector{Int}} = nothing,
    strategy::Symbol = :S3,
    sparse_basis::Bool = true,
)
    N, K = size(X)
    isempty(blocks) && error("No blocks provided")
    row_weights !== nothing && length(row_weights) == N ||
        row_weights === nothing ||
        error("row_weights must have length N = $N")

    # Sort blocks for stable CSC row ordering.
    blocks_sorted = [sort!(collect(b)) for b in blocks]

    # Compute per-block spectra. Avoid materializing per-block `U` factors:
    # for tall snapshot matrices (N ≫ K), `svd(Xb)` returns `U` of size `Nᵢ×K`,
    # which can be multiple GB for large grids.
    spectra = Vector{Any}(undef, length(blocks_sorted))
    Sblocks = Vector{Any}(undef, length(blocks_sorted))
    for (ib, rows) in enumerate(blocks_sorted)
        Xb = @view X[rows, :]
        spec = if row_weights === nothing
            pod_spectrum_from_gram(Xb)
        else
            pod_spectrum_from_weighted_gram(Xb, @view row_weights[rows])
        end
        spectra[ib] = spec
        Sblocks[ib] = spec.S
    end

    # Determine local ranks.
    nb = length(blocks_sorted)
    ranks = zeros(Int, nb)
    if local_ranks !== nothing
        length(local_ranks) == nb || error("local_ranks must have length $nb")
        ranks .= Int.(local_ranks)
    elseif strategy == :S3
        for i in 1:nb
            ranks[i] = rank_from_energy_safe(Sblocks[i]; energy_fraction)
        end
    elseif strategy == :S1
        requested_total_rank === nothing &&
            error("strategy=:S1 requires requested_total_rank (e.g. ROMConfig.npod)")
        # Global pooling over block singular values.
        σ = Float64[]
        bid = Int[]
        for i in 1:nb
            Si = Sblocks[i]
            for j in 1:length(Si)
                push!(σ, Float64(Si[j]))
                push!(bid, i)
            end
        end
        nsel = min(requested_total_rank, length(σ))
        order = sortperm(σ; rev = true)
        selected = view(order, 1:nsel)
        for idx in selected
            ranks[bid[idx]] += 1
        end
    else
        error("Unknown space-local strategy: $strategy (use :S1 or :S3)")
    end

    # Clamp ranks to feasible limits.
    for i in 1:nb
        ranks[i] = clamp(ranks[i], 0, length(Sblocks[i]))
    end

    rtot = sum(ranks)
    if rtot == 0
        rom_basis = sparse_basis ? spzeros(eltype(X), N, 0) : zeros(eltype(X), N, 0)
        projector = sparse_basis ? spzeros(eltype(X), 0, N) : zeros(eltype(X), 0, N)
        return (; rom_basis, projector, ranks, blocks = blocks_sorted, strategy, rtot)
    end

    # Assemble basis as block-sparse CSC (or dense with zeros).
    if sparse_basis
        nnz = sum(length(blocks_sorted[i]) * ranks[i] for i in 1:nb)
        colptr = Vector{Int}(undef, rtot + 1)
        rowval = Vector{Int}(undef, nnz)
        nzval = Vector{eltype(X)}(undef, nnz)
        colptr[1] = 1

        col = 1
        off = 1
        for i in 1:nb
            ri = ranks[i]
            ri == 0 && continue
            rows = blocks_sorted[i]
            Ni = length(rows)
            Xb = @view X[rows, :]
            spec = spectra[i]

            for j in 1:ri
                colptr[col + 1] = colptr[col] + Ni
                copyto!(view(rowval, off:(off + Ni - 1)), rows)
                mode = pod_mode_from_gram(Xb, spec, j)
                copyto!(view(nzval, off:(off + Ni - 1)), mode)
                off += Ni
                col += 1
            end
        end

        rom_basis = SparseMatrixCSC(N, rtot, colptr, rowval, nzval)
    else
        rom_basis = zeros(eltype(X), N, rtot)
        col = 1
        for i in 1:nb
            ri = ranks[i]
            ri == 0 && continue
            rows = blocks_sorted[i]
            Xb = @view X[rows, :]
            spec = spectra[i]
            for j in 1:ri
                mode = pod_mode_from_gram(Xb, spec, j)
                copyto!(view(rom_basis, rows, col), mode)
                col += 1
            end
        end
    end

    projector = row_weights === nothing ? rom_basis' : (rom_basis' * Diagonal(row_weights))

    return (; rom_basis, projector, ranks, blocks = blocks_sorted, strategy, rtot)
end

end
