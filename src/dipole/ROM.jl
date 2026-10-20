using Dates
using JLD2
using KernelAbstractions
using LinearAlgebra
import IncompressibleNavierStokes as NS

include(joinpath(@__DIR__, "..", "workflows", "QuadraticROM.jl"))

function _pod_energy_retained(S, r::Integer)
    isempty(S) && return 0.0
    r < 1 && return 0.0
    r_eff = min(Int(r), length(S))
    tot = sum(abs2, S)
    tot <= 0 && return 0.0
    return sum(abs2, @view(S[1:r_eff])) / tot
end

const WORKFLOW_JLD2_IOTYPE = IOStream

function save_object_atomic(path::AbstractString, obj; iotype::Type = WORKFLOW_JLD2_IOTYPE)
    tmp = path * ".tmp"
    if isfile(tmp)
        rm(tmp; force = true)
    end
    save_object(tmp, obj; iotype = iotype)
    mv(tmp, path; force = true)
    return nothing
end

function load_object_iostream(path::AbstractString; iotype::Type = WORKFLOW_JLD2_IOTYPE)
    return jldopen(path, "r"; iotype = iotype) do file
        all_keys = collect(keys(file))
        isempty(all_keys) && error("File $path does not contain any object")
        length(all_keys) > 1 && error("File $path contains more than one object; use jldopen directly.")
        return file[all_keys[1]]
    end
end

"""
Compute the ROM basis and projector from FOM snapshots.

Returns a named tuple with `rom_basis`, `projector`, singular values, and
selected basis indices.
"""
function build_rom_basis(
    snapshots_u::AbstractMatrix,
    setup,
    romcfg::ROMConfig;
    requested_rank::Int,
)
    Kfull = size(snapshots_u, 2)
    idx_basis = 1:Kfull
    Xbasis = snapshots_u

    Kbasis = size(Xbasis, 2)

    # Weighted inner product that matches the discrete kinetic-energy norm on
    # stretched/staggered grids.
    w = velocity_mass_vector(setup, NS.vectorfield(setup))
    length(w) == size(Xbasis, 1) || error("Velocity mass vector length mismatch with snapshot dimension")

    basis_mode = lowercase(romcfg.basis_mode)
    if basis_mode == "global"
        spectrum = pod_spectrum_from_weighted_gram(Xbasis, w)
        rank_sel_raw = lowercase(strip(romcfg.rank_selection))
        rank_selection = if rank_sel_raw in ("npod", "fixed", "rank")
            "npod"
        elseif rank_sel_raw in ("energy_fraction", "energyfraction", "energy", "fraction", "eta")
            "energy_fraction"
        else
            error("Unsupported rank_selection: $(romcfg.rank_selection). Use \"npod\" or \"energy_fraction\".")
        end

        rfit = if rank_selection == "npod"
            r = min(requested_rank, Kbasis)
            if r < requested_rank
                @warn "Requested npod = $(requested_rank), but only $Kbasis basis snapshots are available. Using rfit = $r."
            end
            r
        else
            ef = romcfg.energy_fraction
            (0.0 < ef <= 1.0) || error("energy_fraction must be in (0, 1] for rank_selection=\"energy_fraction\", got $ef")
            r = rank_from_energy_safe(spectrum.S; energy_fraction = ef)
            r >= 1 || error("Energy-based rank selection returned r = $r; check snapshot energy content.")
            r
        end

        energy_retained = _pod_energy_retained(spectrum.S, rfit)
        rom_basis = Matrix{eltype(Xbasis)}(undef, size(Xbasis, 1), rfit)
        for i in 1:rfit
            copyto!(@view(rom_basis[:, i]), pod_mode_from_gram(Xbasis, spectrum, i))
        end
        projector = rom_basis' * Diagonal(w)
        rank_info = (;
            rank_selection = rank_selection,
            requested_npod = requested_rank,
            energy_fraction = romcfg.energy_fraction,
            energy_retained = energy_retained,
            rfit = rfit,
        )
        return (; rom_basis, projector, singular_values = spectrum.S, idx_basis, rank_info)
    elseif basis_mode == "space_local"
        requested_total_rank = min(requested_rank, Kbasis)
        if requested_total_rank < requested_rank
            @warn "Requested npod = $(requested_rank), but only $Kbasis basis snapshots are available. Using rfit = $requested_total_rank."
        end
        xcoords, ycoords = velocity_xycoords(setup; T = Float64)
        part_method = lowercase(romcfg.basis_partition_method)
        nx = romcfg.basis_partition_nx
        ny = romcfg.basis_partition_ny
        nx >= 1 || error("basis_partition_nx must be >= 1, got $nx")
        ny >= 1 || error("basis_partition_ny must be >= 1, got $ny")

        blocks, partinfo = if part_method == "uniform_x"
            xcuts = uniform_cuts(xcoords, nx)
            blocks_from_xycoords(xcoords, nothing; xcuts = xcuts)
        elseif part_method == "uniform_xy"
            xcuts = uniform_cuts(xcoords, nx)
            ycuts = uniform_cuts(ycoords, ny)
            blocks_from_xycoords(xcoords, ycoords; xcuts = xcuts, ycuts = ycuts)
        elseif part_method == "manual_x"
            xcuts = isnothing(romcfg.basis_partition_xcuts) ? Float64[] : romcfg.basis_partition_xcuts
            blocks_from_xycoords(xcoords, nothing; xcuts = xcuts)
        elseif part_method == "manual_xy"
            xcuts = isnothing(romcfg.basis_partition_xcuts) ? Float64[] : romcfg.basis_partition_xcuts
            ycuts = isnothing(romcfg.basis_partition_ycuts) ? Float64[] : romcfg.basis_partition_ycuts
            blocks_from_xycoords(xcoords, ycoords; xcuts = xcuts, ycuts = ycuts)
        else
            error(
                "Unsupported basis_partition_method: $(romcfg.basis_partition_method). " *
                "Use \"uniform_x\", \"uniform_xy\", \"manual_x\", or \"manual_xy\".",
            )
	        end

	        strategy = Symbol(romcfg.basis_local_strategy)
	        slpod = build_space_local_basis(
	            Xbasis,
	            blocks;
	            requested_total_rank = requested_total_rank,
                row_weights = w,
	            energy_fraction = romcfg.basis_local_energy_fraction,
	            local_ranks = romcfg.basis_local_ranks,
	            strategy = strategy,
	            sparse_basis = romcfg.basis_space_local_sparse,
	        )
	        if strategy == :S3 && romcfg.basis_local_ranks === nothing && slpod.rtot != requested_rank
	            @info "Space-local POD (S3) selected an effective rank from the energy threshold" (
	                requested_npod = requested_rank,
	                rfit_effective = slpod.rtot,
	                basis_local_energy_fraction = romcfg.basis_local_energy_fraction,
	            )
	        end

	        rom_basis = slpod.rom_basis
	        projector = slpod.projector
	        singular_values = nothing
	        space_local_info = (;
	            basis_mode = "space_local",
	            basis_local_strategy = String(romcfg.basis_local_strategy),
	            basis_partition_method = String(romcfg.basis_partition_method),
	            basis_partition = partinfo,
	            basis_local_ranks = slpod.ranks,
	            basis_local_energy_fraction = romcfg.basis_local_energy_fraction,
	            basis_space_local_sparse = romcfg.basis_space_local_sparse,
	            rtot = slpod.rtot,
	        )
	        return (; rom_basis, projector, singular_values, idx_basis, space_local_info)
	    else
	        error("Unsupported basis_mode: $(romcfg.basis_mode). Use \"global\" or \"space_local\".")
    end
end

"""
Generate exact-OpInf training states and derivatives by one-step FOM evolution.

Returns `(Atrain, Dtrain, method_info)` in reduced coordinates.
"""
function generate_exact_opinf_snapshots(
    rom_basis::AbstractMatrix,
    basis_projector::AbstractMatrix,
    setup,
    psolver,
    fomcfg::FOMConfig,
    romcfg::ROMConfig;
    dt_train::Real,
)
    rfit = size(rom_basis, 2)
    T = eltype(rom_basis)
    Atrain = rank_ensuring_states_quadratic(rfit; T = T)
    K = size(Atrain, 2)
    Dtrain = zeros(T, rfit, K)

    # Exact OpInf theory assumes a single explicit-Euler-consistent step.
    method = romcfg.exact_opinf_force_explicit_euler ? :explicit_euler : :lmwray3
    method_info = nothing

    for s in 1:K
        a = @view Atrain[:, s]
        u0 = rom_basis * a
        u1, method_info = one_step_fom_from_state(u0, setup, psolver, fomcfg; dt = dt_train, method = method)
        Dtrain[:, s] = basis_projector * ((u1 - u0) / dt_train)
    end

    return Atrain, Dtrain, method_info
end

function _copy_basis_column_to_velocity!(u, rom_basis::AbstractMatrix, j::Integer)
    mode = Array(rom_basis[:, j])
    copyto!(u, reshape(mode, size(u)))
    return u
end

function _copy_basis_pair_sum_to_velocity!(u, rom_basis::AbstractMatrix, j::Integer, k::Integer)
    mode = Array(rom_basis[:, j])
    mode .+= Array(rom_basis[:, k])
    copyto!(u, reshape(mode, size(u)))
    return u
end

function _project_velocity_rhs(
    force,
    basis_projector::AbstractMatrix,
    setup,
    psolver,
    p,
)
    NS.project!(force, setup; psolver, p)
    return basis_projector * vec(Array(force))
end

function _symmetric_tensor_from_unique_coefficients(B::AbstractMatrix, pairs)
    r = size(B, 1)
    Q = zeros(eltype(B), r, r, r)
    @inbounds for i in 1:r
        for (m, (j, k)) in enumerate(pairs)
            coeff = B[i, m]
            if j == k
                Q[i, j, k] = coeff
            else
                Q[i, j, k] = coeff / 2
                Q[i, k, j] = coeff / 2
            end
        end
    end
    return Q
end

"""
Assemble the intrusive dipole ROM operators by direct reduced tensor assembly.

The full-order polynomial operators are not formed. Instead, each reduced
operator column is assembled by applying the FOM diffusion or convection
operator to POD modes, applying the same velocity projection used by the FOM,
and projecting the result onto the POD space.
"""
function assemble_intrusive_rom_operators(
    rom_basis::AbstractMatrix,
    basis_projector::AbstractMatrix,
    setup,
    psolver,
    fomcfg::FOMConfig,
)
    rfit = size(rom_basis, 2)
    T = eltype(rom_basis)

    pairs = quadratic_pairs(rfit)
    nq = length(pairs)
    L = zeros(T, rfit, rfit)
    B = zeros(T, rfit, nq)

    u = NS.vectorfield(setup)
    usum = similar(u)
    force = similar(u)
    p = NS.scalarfield(setup)

    projected_self_convection = Vector{Vector{T}}(undef, rfit)
    n_mixed_pairs = rfit * (rfit - 1) ÷ 2
    n_projected_evaluations = 2 * rfit + n_mixed_pairs
    @info "Assembling direct intrusive ROM operators" (
        rfit = rfit,
        n_unique_quadratic = nq,
        n_mixed_pairs = n_mixed_pairs,
        n_projected_operator_evaluations = n_projected_evaluations,
    )
    prog = Progress(rfit + nq; desc = "Assembling intrusive ROM")

    for j in 1:rfit
        _copy_basis_column_to_velocity!(u, rom_basis, j)
        NS.apply_bc_u!(u, 0.0, setup)

        fill!(force, zero(eltype(force)))
        NS.diffusion!(force, u, setup, fomcfg.viscosity)
        L[:, j] = _project_velocity_rhs(force, basis_projector, setup, psolver, p)

        fill!(force, zero(eltype(force)))
        NS.convection!(force, u, setup)
        projected_self_convection[j] = _project_velocity_rhs(force, basis_projector, setup, psolver, p)
        next!(prog; showvalues = [("stage", "linear/self"), ("mode", j), ("rank", rfit)])
    end

    for (m, (j, k)) in enumerate(pairs)
        if j == k
            B[:, m] = projected_self_convection[j]
        else
            _copy_basis_pair_sum_to_velocity!(usum, rom_basis, j, k)
            NS.apply_bc_u!(usum, 0.0, setup)

            fill!(force, zero(eltype(force)))
            NS.convection!(force, usum, setup)
            projected_pair_convection = _project_velocity_rhs(force, basis_projector, setup, psolver, p)
            B[:, m] = projected_pair_convection - projected_self_convection[j] - projected_self_convection[k]
        end
        next!(prog; showvalues = [("stage", "quadratic"), ("pair", m), ("total_pairs", nq)])
    end

    Q = _symmetric_tensor_from_unique_coefficients(B, pairs)

    method_info = (
        method = "intrusive_direct_reduced_operator_assembly",
        is_intrusive = true,
        assembly = "direct_projected_diffusion_convection",
        n_linear_columns = rfit,
        n_quadratic_columns = nq,
    )
    diagnostics = (
        maxabsL = maximum(abs, L),
        maxabsQ = maximum(abs, Q),
    )

    return (;
        L,
        Q,
        B,
        pairs,
        σ = T[],
        resnorm = zero(T),
        relres = zero(T),
        diagnostics,
        method_info,
    )
end

"""
Check whether a ROM run directory contains the expected output files.
"""
function is_rom_outputs_complete(dir::AbstractString)
    needed = [
        joinpath(dir, "config.toml"),
        joinpath(dir, "fit.jld2"),
        joinpath(dir, "trajectory.jld2"),
        joinpath(dir, "metadata.jld2"),
    ]
    return all(isfile, needed)
end

"""
Run the dipole ROM workflow for a parent FOM run.

This routine loads FOM data, builds a basis, fits a quadratic ROM (trajectory,
exact OpInf, or intrusive mode), simulates the ROM trajectory, and writes
`fit.jld2`, `trajectory.jld2`, and `metadata.jld2`.
"""
function run_rom(parent_fom_id::String, romcfg::ROMConfig; force::Bool = false, skip_if_exists::Bool = true)
    run_id = rom_id(parent_fom_id, romcfg)
    dir = rom_dir(run_id)

    if isdir(dir)
        if skip_if_exists
            if is_rom_outputs_complete(dir)
                # Keep existing outputs but ensure parent id is recorded.
                open(joinpath(dir, "parent_fom_id.txt"), "w") do io
                    write(io, parent_fom_id)
                    write(io, "\n")
                end
                if !isfile(joinpath(dir, "rom_dim.txt")) && isfile(joinpath(dir, "metadata.jld2"))
                    meta = load_object_iostream(joinpath(dir, "metadata.jld2"))
                    if hasproperty(meta, :rfit)
                        open(joinpath(dir, "rom_dim.txt"), "w") do io
                            write(io, string(meta.rfit))
                            write(io, "\n")
                        end
                    end
                end
                if !isfile(joinpath(dir, "rom_rank_info.txt")) && isfile(joinpath(dir, "metadata.jld2"))
                    meta = load_object_iostream(joinpath(dir, "metadata.jld2"))
                    basis_mode = hasproperty(meta, :basis_mode) ? string(meta.basis_mode) : "unknown"
                    rank_selection = hasproperty(meta, :rank_selection) ? string(meta.rank_selection) : "unknown"
                    requested_npod = hasproperty(meta, :requested_npod) ? meta.requested_npod : nothing
                    energy_fraction = hasproperty(meta, :energy_fraction) ? meta.energy_fraction : nothing
                    energy_retained =
                        hasproperty(meta, :rank_info) && meta.rank_info !== nothing && hasproperty(meta.rank_info, :energy_retained) ?
                        meta.rank_info.energy_retained : nothing
                    rfit_meta = hasproperty(meta, :rfit) ? meta.rfit : nothing
                    open(joinpath(dir, "rom_rank_info.txt"), "w") do io
                        write(io, "basis_mode = $basis_mode\n")
                        write(io, "rank_selection = $rank_selection\n")
                        if requested_npod !== nothing
                            write(io, "requested_npod = $requested_npod\n")
                        end
                        if energy_fraction !== nothing
                            write(io, "energy_fraction = $energy_fraction\n")
                        end
                        if energy_retained !== nothing
                            write(io, "energy_retained = $(round(energy_retained; digits = 10))\n")
                        end
                        if rfit_meta !== nothing
                            write(io, "rfit = $rfit_meta\n")
                        end
                    end
                end
                @info "ROM run already exists, reusing $run_id"
                return run_id
            else
                @warn "ROM run directory exists but outputs missing; rerunning" run_id dir
                rm(dir; recursive = true, force = true)
            end
        elseif !force
            error("ROM directory already exists: $dir")
        else
            rm(dir; recursive = true, force = true)
        end
    end

    ensure_run_dirs!(dir)
    save_rom_config(joinpath(dir, "config.toml"), romcfg)
    open(joinpath(dir, "parent_fom_id.txt"), "w") do io
        write(io, parent_fom_id)
        write(io, "\n")
    end

    fomcfg = load_fom_config_from_run(parent_fom_id)
    snapshots = load_fom_snapshots(parent_fom_id)

    Kfull = size(snapshots.u, 2)
    use_intrusive_rom_effective = romcfg.use_intrusive_rom
    use_exact_opinf_effective = !use_intrusive_rom_effective && romcfg.use_exact_opinf
    if use_exact_opinf_effective && !romcfg.exact_opinf_force_explicit_euler && romcfg.fallback_to_standard_opinf
        @warn "Falling back to trajectory-based OpInf because explicit Euler was not forced." (
            exact_opinf_force_explicit_euler = romcfg.exact_opinf_force_explicit_euler,
            fallback_to_standard_opinf = romcfg.fallback_to_standard_opinf,
        )
        use_exact_opinf_effective = false
    end
    train_mode = use_intrusive_rom_effective ? :intrusive : (use_exact_opinf_effective ? :exact_opinf : :traj_opinf)

    Atrain = nothing
    Dtrain = nothing
    method_info = nothing
    train_dt = nothing
    train_snapshot_count = nothing

    setup_basis = getsetup(fomcfg; backend = KernelAbstractions.CPU())
    basis = build_rom_basis(
        snapshots.u,
        setup_basis,
        romcfg;
        requested_rank = romcfg.npod,
    )

    rom_basis = basis.rom_basis
    basis_projector = basis.projector
    rfit = size(rom_basis, 2)
    basis_space_local_info = hasproperty(basis, :space_local_info) ? basis.space_local_info : nothing
    rank_info = hasproperty(basis, :rank_info) ? basis.rank_info : nothing
    a_exact_full = basis_projector * snapshots.u
    dadt_exact_full = basis_projector * snapshots.dudt

    @info "ROM basis selection" (
        basis_snapshot_start = first(basis.idx_basis),
        basis_snapshot_end = last(basis.idx_basis),
        nbasis = length(basis.idx_basis),
    )

    if !isnothing(rank_info) && hasproperty(rank_info, :energy_retained)
        @info "ROM POD energy retained" (
            rank_selection = rank_info.rank_selection,
            requested_npod = rank_info.requested_npod,
            rfit = rank_info.rfit,
            energy_retained = rank_info.energy_retained,
        )
    end

    if train_mode in (:intrusive, :exact_opinf)
        a_exact = a_exact_full
        dadt_exact = dadt_exact_full

        backend = getbackend()
        @info "ROM full-order operator backend" backend = string(typeof(backend))
        ensure_cuda_psolver!(backend)
        setup = getsetup(fomcfg; backend = backend)
        psolver = NS.default_psolver(setup)

        if train_mode == :intrusive
            fit = assemble_intrusive_rom_operators(
                rom_basis,
                basis_projector,
                setup,
                psolver,
                fomcfg,
            )
            method_info = fit.method_info
            @info "Intrusive ROM training diagnostics" (
                K_train = 0,
                nf = rfit + n_unique_quadratic(rfit),
                n_linear_columns = method_info.n_linear_columns,
                n_quadratic_columns = method_info.n_quadratic_columns,
                resnorm = fit.resnorm,
                relres = fit.relres,
                method = method_info,
            )
        ## Exact OpInf ##
        else
            snapshot_dt = fomcfg.tsim / fomcfg.nsnapshot
            estimated_dt = estimate_exact_opinf_dt_matlab(a_exact, dadt_exact; fallback_dt = snapshot_dt)
            train_dt = min(estimated_dt, snapshot_dt)
            if train_dt < estimated_dt
                @info "Capped exact-OpInf training step at the FOM snapshot spacing" estimated_dt train_dt snapshot_dt
            else
                @info "Estimated exact-OpInf training step" train_dt snapshot_dt
            end
            Atrain, Dtrain, method_info = generate_exact_opinf_snapshots(
                rom_basis,
                basis_projector,
                setup,
                psolver,
                fomcfg,
                romcfg;
                dt_train = train_dt,
            )

            fit = fit_quadratic_rom_exact_symmetric(Atrain, Dtrain)
            if !romcfg.exact_opinf_force_explicit_euler
                @warn "Exact OpInf training did not use explicit Euler; exact recovery assumptions are not fully satisfied." method = method_info
            end

            @info "Exact OpInf training diagnostics" (
                K_train = size(Atrain, 2),
                nf = rfit + n_unique_quadratic(rfit),
                condPhi = fit.σ[1] / max(fit.σ[end], eps(eltype(fit.σ))),
                resnorm = fit.resnorm,
                relres = fit.relres,
                method = method_info,
            )
        end
        train_snapshot_count = size(a_exact, 2)
        ##
    else
        a_exact = a_exact_full
        dadt_exact = dadt_exact_full
        idx_train = basis.idx_basis
        a_train = @view a_exact[:, idx_train]
        dadt_train = @view dadt_exact[:, idx_train]

        K = size(a_train, 2)
        rmax = max_quadratic_rank_symmetric(K)
        rfit = min(romcfg.npod, rmax)
        if rfit < romcfg.npod
            @warn "Requested npod = $(romcfg.npod), but with K = $K only up to r = $rmax is used. Using rfit = $rfit."
        end

        fit = fit_quadratic_rom_symmetric(
            @view(a_train[1:rfit, :]),
            @view(dadt_train[1:rfit, :]);
            lambda = romcfg.lambda,
        )

        rom_basis = @view rom_basis[:, 1:rfit]
        basis_projector = @view basis_projector[1:rfit, :]
        a_exact = @view a_exact[1:rfit, :]
        dadt_exact = @view dadt_exact[1:rfit, :]
        train_snapshot_count = size(a_train, 2)
    end
    ## General ROM ##
    a0 = copy(a_exact[:, 1])
    a_rom = simulate_rom_adaptive(
        a0,
        fit.L,
        fit.Q;
        tsim = fomcfg.tsim,
        nsave = fomcfg.nsnapshot,
        safety = romcfg.safety,
        growthfac = romcfg.growthfac,
    )

    condΦ = isempty(fit.σ) ? convert(eltype(fit.L), NaN) :
        fit.σ[1] / max(fit.σ[end], eps(eltype(fit.σ)))
    diagnostics = (
        rfit = rfit,
        condPhi = condΦ,
        maxabsL = maximum(abs, fit.L),
        maxabsQ = maximum(abs, fit.Q),
    )

    fit_data = (
        rom_basis = rom_basis,
        basis_projector = basis_projector,
        basis_idx = basis.idx_basis,
        L = fit.L,
        Q = fit.Q,
        B = fit.B,
        σ = fit.σ,
        pairs = fit.pairs,
        train_mode = String(train_mode),
        run_tag = romcfg.run_tag,
        use_intrusive_rom = use_intrusive_rom_effective,
        use_intrusive_rom_config = romcfg.use_intrusive_rom,
        use_exact_opinf = use_exact_opinf_effective,
        use_exact_opinf_config = romcfg.use_exact_opinf,
        train_dt = use_exact_opinf_effective ? train_dt : nothing,
        train_residual_norm = hasproperty(fit, :resnorm) ? fit.resnorm : nothing,
        train_residual_rel = hasproperty(fit, :relres) ? fit.relres : nothing,
        Atrain = Atrain,
        Dtrain = Dtrain,
        train_method = method_info,
        diagnostics = diagnostics,
    )
    traj_data = (
        a_exact = a_exact,
        dadt_exact = dadt_exact,
        a_rom = a_rom,
    )
    metadata = (
        created_at = string(Dates.now()),
        parent_fom_id = parent_fom_id,
        run_tag = romcfg.run_tag,
        requested_npod = romcfg.npod,
        rfit = rfit,
        basis_mode = lowercase(romcfg.basis_mode),
        rank_selection = isnothing(rank_info) ? "unknown" : rank_info.rank_selection,
        rank_info = rank_info,
        basis_space_local_info = basis_space_local_info,
        diagnostics = diagnostics,
        train_mode = String(train_mode),
        use_intrusive_rom = use_intrusive_rom_effective,
        use_intrusive_rom_config = romcfg.use_intrusive_rom,
        use_exact_opinf = use_exact_opinf_effective,
        use_exact_opinf_config = romcfg.use_exact_opinf,
        exact_opinf_force_explicit_euler = romcfg.exact_opinf_force_explicit_euler,
        exact_opinf_assumptions_satisfied = use_exact_opinf_effective && romcfg.exact_opinf_force_explicit_euler,
        dt_train = use_exact_opinf_effective ? train_dt : nothing,
        K_train = train_mode == :intrusive ? 0 :
            (train_mode == :exact_opinf ? (rfit + n_unique_quadratic(rfit)) : train_snapshot_count),
        basis_snapshot_start = first(basis.idx_basis),
        basis_snapshot_end = last(basis.idx_basis),
        basis_nsnapshot = length(basis.idx_basis),
    )

    open(joinpath(dir, "rom_dim.txt"), "w") do io
        # Effective ROM dimension (can differ from requested npod for space-local S3).
        write(io, string(rfit))
        write(io, "\n")
    end
    open(joinpath(dir, "rom_rank_info.txt"), "w") do io
        write(io, "basis_mode = $(lowercase(romcfg.basis_mode))\n")
        if !isnothing(romcfg.run_tag)
            write(io, "run_tag = $(romcfg.run_tag)\n")
        end
        if !isnothing(rank_info)
            write(io, "rank_selection = $(rank_info.rank_selection)\n")
        else
            write(io, "rank_selection = unknown\n")
        end
        write(io, "requested_npod = $(romcfg.npod)\n")
        write(io, "energy_fraction = $(romcfg.energy_fraction)\n")
        if !isnothing(rank_info) && hasproperty(rank_info, :energy_retained)
            write(io, "energy_retained = $(round(rank_info.energy_retained; digits = 10))\n")
        end
        write(io, "rfit = $(rfit)\n")
    end

    save_object_atomic(joinpath(dir, "fit.jld2"), fit_data)
    save_object_atomic(joinpath(dir, "trajectory.jld2"), traj_data)
    save_object_atomic(joinpath(dir, "metadata.jld2"), metadata)

    register_rom!(run_id, parent_fom_id, romcfg, fomcfg; dt_train = train_dt)
    @info "Saved ROM run to $dir"

    return run_id
end

"""Load fitted ROM operators and diagnostics for `run_id`."""
load_rom_fit(run_id::String) = load_object_iostream(joinpath(rom_dir(run_id), "fit.jld2"))
"""Load projected and simulated ROM trajectories for `run_id`."""
load_rom_trajectory(run_id::String) = load_object_iostream(joinpath(rom_dir(run_id), "trajectory.jld2"))
"""Load ROM metadata for `run_id`."""
load_rom_metadata(run_id::String) = load_object_iostream(joinpath(rom_dir(run_id), "metadata.jld2"))
"""Load the persisted ROM config TOML for `run_id`."""
load_rom_config_from_run(run_id::String) = load_rom_config(joinpath(rom_dir(run_id), "config.toml"))
