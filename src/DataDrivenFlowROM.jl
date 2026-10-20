"""
    DataDrivenFlowROM

Reusable algorithms for POD and affine-quadratic reduced-order modelling.
The case-study orchestration lives under `experiments/` and is intentionally
kept separate from this small public API.
"""
module DataDrivenFlowROM

include("workflows/FOM.jl")
include("workflows/POD.jl")
include("workflows/QuadraticROM.jl")

export snapshot_stream_paths,
    has_streamed_snapshots,
    load_streamed_snapshots,
    one_step_fom_from_state,
    pod_spectrum_from_gram,
    pod_spectrum_from_weighted_gram,
    pod_mode_from_gram,
    cumulative_energy,
    rank_from_energy,
    choose_rank,
    rank_from_energy_safe,
    build_space_local_basis,
    uniform_cuts,
    bins_from_cuts,
    blocks_from_xycoords,
    velocity_xycoords,
    scalar_xycoords,
    n_unique_quadratic,
    n_features_affine_quadratic,
    max_quadratic_rank_symmetric,
    max_affine_quadratic_rank_symmetric,
    quadratic_pairs,
    quadratic_features_symmetric,
    quadratic_features_affine_symmetric,
    fit_quadratic_rom_symmetric,
    fit_affine_quadratic_rom_symmetric,
    rank_ensuring_states_quadratic,
    rank_ensuring_states_affine_quadratic,
    fit_quadratic_rom_exact_symmetric,
    fit_affine_quadratic_rom_exact_symmetric,
    simulate_rom_adaptive

end

