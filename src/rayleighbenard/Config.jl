using TOML

Base.@kwdef struct FOMConfig
    viscosity::Float64 = 2.5e-4
    conductivity::Float64 = 2.5e-4
    gravity::Float64 = 1.0
    gdir::Int = 2
    dodissipation::Bool = true

    nx::Int = 256
    ny::Int = 128
    xlims::NTuple{2,Float64} = (0.0, 2.0)
    ylims::NTuple{2,Float64} = (0.0, 1.0)
    stretch_x::Union{Nothing,Float64} = nothing
    stretch_y::Union{Nothing,Float64} = nothing

    temp_hot::Float64 = 1.0
    temp_cold::Float64 = 0.0
    # Initial temperature profile:
    # - "conduction": linear conduction state (matches wall temperatures)
    # - "uniform": uniform mean temperature (walls impose gradient at t=0)
    initial_temperature_profile::String = "uniform"
    perturbation_amplitude::Float64 = 1e-2
    perturbation_wavenumber::Float64 = 20.0

    tsim::Float64 = 20.0
    cfl::Float64 = 0.9
    nsnapshot::Int = 200
end

Base.@kwdef struct ROMConfig
    # Optional non-behavioral tag used only to create a distinct ROM run id.
    run_tag::Union{Nothing,String} = nothing
    npod::Int = 20
    # Rank selection strategy for global POD basis:
    # - "npod": use `npod` (fixed rank)
    # - "energy_fraction": choose the smallest rank capturing `energy_fraction`
    rank_selection::String = "npod"
    # Build global POD modes from mean-subtracted snapshots. This is independent
    # from the rank-selection strategy.
    basis_centered::Bool = true
    lambda::Float64 = 1e-2
    safety::Float64 = 0.03
    growthfac::Float64 = 1.05
    energy_fraction::Float64 = 0.99
    # Ignore the initial FOM transient when constructing/training the ROM.
    # The first saved FOM snapshot at or after this time is used.
    training_start_time::Float64 = 0.0
    # Optional end of the FOM training window. The last saved snapshot at or
    # before this time is used. `nothing` means train through the final FOM
    # snapshot.
    training_end_time::Union{Nothing,Float64} = nothing
    use_intrusive_rom::Bool = false
    use_exact_opinf::Bool = true
    exact_opinf_force_explicit_euler::Bool = true
    fallback_to_standard_opinf::Bool = true
    # Include an affine (constant) term `c` in the quadratic ROM:
    #   da/dt = c + L*a + q(a).
    #
    # Rayleigh-Bénard simulations with Dirichlet temperature boundary
    # conditions generally induce an affine term in reduced coordinates unless
    # the state is lifted/centered around an equilibrium.
    include_constant_term::Bool = true

    # POD basis mode:
    # - "global": standard POD (current default)
    # - "space_local": space-local POD on row blocks
    basis_mode::String = "global"
    basis_local_strategy::String = "S3" # "S1" (fixed total rank) or "S3" (per-block energy)
    basis_partition_method::String = "uniform_x" # "uniform_x", "uniform_xy", "manual_x", or "manual_xy"
    basis_partition_nx::Int = 1
    basis_partition_ny::Int = 1
    basis_partition_xcuts::Union{Nothing,Vector{Float64}} = nothing
    basis_partition_ycuts::Union{Nothing,Vector{Float64}} = nothing
    basis_local_energy_fraction::Float64 = 0.99
    basis_local_ranks::Union{Nothing,Vector{Int}} = nothing
    basis_space_local_sparse::Bool = true
end

Base.@kwdef struct SwitchedROMConfig
    # Optional non-behavioral tag used only to create a distinct switched ROM id.
    run_tag::Union{Nothing,String} = nothing
    # Existing local ROM run ids ordered as [pre, transition, post].
    local_rom_ids::Vector{String} = String[]
    regime_names::Vector{String} = ["pre", "transition", "post"]
    start_time::Float64 = 0.0
    # `nothing` means predict through the final FOM snapshot.
    end_time::Union{Nothing,Float64} = nothing
    initial_regime::String = "pre"

    # Scheduled switching times, usually derived from the reference FOM.
    switching_policy::String = "scheduled"
    min_switch_time::Float64 = 20.0
    pre_to_transition_threshold::Float64 = 0.15
    transition_to_post_threshold::Float64 = -0.30
    transition_to_pre_threshold::Float64 = 0.20
    post_to_transition_threshold::Float64 = -0.20
    allow_backward_switches::Bool = false
    pre_to_transition_time::Union{Nothing,Float64} = nothing
    transition_to_post_time::Union{Nothing,Float64} = nothing

    # Allow local ROMs trained on another FOM run if the FOM configurations
    # match except for tsim/nsnapshot. This is used for out-of-training-horizon
    # prediction tests against a longer reference FOM.
    allow_parent_fom_mismatch::Bool = false

    safety::Float64 = 0.03
    growthfac::Float64 = 1.05
end

function rom_training_window(
    fomcfg::FOMConfig,
    nsaved::Integer,
    requested_start_time::Real,
    requested_end_time = nothing,
)
    expected = fomcfg.nsnapshot + 1
    nsaved == expected || error("Expected $expected FOM snapshots from the config, found $nsaved")

    t_requested = Float64(requested_start_time)
    isfinite(t_requested) || error("training_start_time must be finite, got $requested_start_time")
    0.0 <= t_requested < fomcfg.tsim ||
        error("training_start_time must satisfy 0 <= training_start_time < tsim=$(fomcfg.tsim), got $t_requested")

    t_end_requested = if isnothing(requested_end_time)
        fomcfg.tsim
    else
        t_end = Float64(requested_end_time)
        isfinite(t_end) || error("training_end_time must be finite, got $requested_end_time")
        t_requested < t_end <= fomcfg.tsim ||
            error("training_end_time must satisfy training_start_time < training_end_time <= tsim=$(fomcfg.tsim), got $t_end")
        t_end
    end

    all_times = collect(range(0.0, fomcfg.tsim; length = expected))
    first_snapshot = searchsortedfirst(all_times, t_requested)
    first_snapshot <= nsaved || error("No FOM snapshot is available at or after training_start_time=$t_requested")
    last_snapshot = searchsortedlast(all_times, t_end_requested)
    last_snapshot >= first_snapshot ||
        error("No FOM snapshot interval is available between training_start_time=$t_requested and training_end_time=$t_end_requested")

    snapshot_indices = first_snapshot:last_snapshot
    length(snapshot_indices) >= 2 ||
        error("The ROM training window must contain at least two snapshots; got $(length(snapshot_indices))")

    times = all_times[snapshot_indices]
    actual_start_time = first(times)
    if actual_start_time > t_requested
        @info "Rounded ROM training start forward to the next stored FOM snapshot" (
            requested_start_time = t_requested,
            actual_start_time = actual_start_time,
        )
    end
    actual_end_time = last(times)
    if actual_end_time < t_end_requested
        @info "Rounded ROM training end backward to the previous stored FOM snapshot" (
            requested_end_time = t_end_requested,
            actual_end_time = actual_end_time,
        )
    end

    return (;
        requested_start_time = t_requested,
        requested_end_time = t_end_requested,
        actual_start_time,
        end_time = actual_end_time,
        duration = actual_end_time - actual_start_time,
        snapshot_indices,
        times,
    )
end

function validate_switched_rom_config(cfg::SwitchedROMConfig)
    length(cfg.local_rom_ids) == 3 ||
        error("SwitchedROMConfig currently requires exactly three local_rom_ids: pre, transition, post")
    length(cfg.regime_names) == 3 ||
        error("SwitchedROMConfig currently requires exactly three regime_names")
    length(unique(cfg.regime_names)) == 3 ||
        error("SwitchedROMConfig regime_names must be unique")
    cfg.initial_regime in cfg.regime_names ||
        error("initial_regime=$(cfg.initial_regime) is not in regime_names=$(cfg.regime_names)")
    isfinite(cfg.start_time) || error("start_time must be finite, got $(cfg.start_time)")
    cfg.start_time >= 0.0 || error("start_time must be nonnegative, got $(cfg.start_time)")
    if !isnothing(cfg.end_time)
        isfinite(cfg.end_time) || error("end_time must be finite, got $(cfg.end_time)")
        cfg.end_time > cfg.start_time || error("end_time must be greater than start_time")
    end
    cfg.safety > 0.0 || error("safety must be positive, got $(cfg.safety)")
    cfg.growthfac > 1.0 || error("growthfac must be greater than 1, got $(cfg.growthfac)")
    cfg.switching_policy == "scheduled" ||
        error("switching_policy=\"$(cfg.switching_policy)\" is disabled; use switching_policy=\"scheduled\"")
    cfg.transition_to_pre_threshold > cfg.pre_to_transition_threshold ||
        error("transition_to_pre_threshold should exceed pre_to_transition_threshold to provide hysteresis")
    cfg.post_to_transition_threshold > cfg.transition_to_post_threshold ||
        error("post_to_transition_threshold should exceed transition_to_post_threshold to provide hysteresis")
    isnothing(cfg.pre_to_transition_time) &&
        error("pre_to_transition_time is required when switching_policy=\"scheduled\"")
    isnothing(cfg.transition_to_post_time) &&
        error("transition_to_post_time is required when switching_policy=\"scheduled\"")
    isfinite(cfg.pre_to_transition_time) ||
        error("pre_to_transition_time must be finite, got $(cfg.pre_to_transition_time)")
    isfinite(cfg.transition_to_post_time) ||
        error("transition_to_post_time must be finite, got $(cfg.transition_to_post_time)")
    cfg.pre_to_transition_time >= cfg.start_time ||
        error("pre_to_transition_time must be at or after start_time")
    cfg.transition_to_post_time > cfg.pre_to_transition_time ||
        error("transition_to_post_time must be after pre_to_transition_time")
    if !isnothing(cfg.end_time)
        cfg.transition_to_post_time <= cfg.end_time ||
            error("transition_to_post_time must be at or before end_time")
    end
    return cfg
end

_stretch_to_toml(x) = isnothing(x) ? "uniform" : x
function _stretch_from_toml(x)
    if x isa AbstractString
        xl = lowercase(x)
        if xl in ("uniform", "none", "nothing")
            return nothing
        end
    end
    return Float64(x)
end

function fom_to_dict(cfg::FOMConfig)
    Dict(
        "viscosity" => cfg.viscosity,
        "conductivity" => cfg.conductivity,
        "gravity" => cfg.gravity,
        "gdir" => cfg.gdir,
        "dodissipation" => cfg.dodissipation,
        "nx" => cfg.nx,
        "ny" => cfg.ny,
        "xlims" => collect(cfg.xlims),
        "ylims" => collect(cfg.ylims),
        "stretch_x" => _stretch_to_toml(cfg.stretch_x),
        "stretch_y" => _stretch_to_toml(cfg.stretch_y),
        "temp_hot" => cfg.temp_hot,
        "temp_cold" => cfg.temp_cold,
        "initial_temperature_profile" => cfg.initial_temperature_profile,
        "perturbation_amplitude" => cfg.perturbation_amplitude,
        "perturbation_wavenumber" => cfg.perturbation_wavenumber,
        "tsim" => cfg.tsim,
        "cfl" => cfg.cfl,
        "nsnapshot" => cfg.nsnapshot,
    )
end

function rom_to_dict(cfg::ROMConfig)
    d = Dict{String,Any}(
        "npod" => cfg.npod,
        "rank_selection" => cfg.rank_selection,
        "basis_centered" => cfg.basis_centered,
        "lambda" => cfg.lambda,
        "safety" => cfg.safety,
        "growthfac" => cfg.growthfac,
        "energy_fraction" => cfg.energy_fraction,
        "training_start_time" => cfg.training_start_time,
        "use_intrusive_rom" => cfg.use_intrusive_rom,
        "use_exact_opinf" => cfg.use_exact_opinf,
        "exact_opinf_force_explicit_euler" => cfg.exact_opinf_force_explicit_euler,
        "fallback_to_standard_opinf" => cfg.fallback_to_standard_opinf,
        "include_constant_term" => cfg.include_constant_term,

        "basis_mode" => cfg.basis_mode,
        "basis_local_strategy" => cfg.basis_local_strategy,
        "basis_partition_method" => cfg.basis_partition_method,
        "basis_partition_nx" => cfg.basis_partition_nx,
        "basis_partition_ny" => cfg.basis_partition_ny,
        "basis_local_energy_fraction" => cfg.basis_local_energy_fraction,
        "basis_space_local_sparse" => cfg.basis_space_local_sparse,
    )

    # TOML does not support null. Only write optional fields if present.
    if !isnothing(cfg.basis_partition_xcuts)
        d["basis_partition_xcuts"] = cfg.basis_partition_xcuts
    end
    if !isnothing(cfg.basis_partition_ycuts)
        d["basis_partition_ycuts"] = cfg.basis_partition_ycuts
    end
    if !isnothing(cfg.basis_local_ranks)
        d["basis_local_ranks"] = cfg.basis_local_ranks
    end
    if !isnothing(cfg.run_tag)
        d["run_tag"] = cfg.run_tag
    end
    if !isnothing(cfg.training_end_time)
        d["training_end_time"] = cfg.training_end_time
    end

    return d
end

function switched_rom_to_dict(cfg::SwitchedROMConfig)
    d = Dict{String,Any}(
        "local_rom_ids" => cfg.local_rom_ids,
        "regime_names" => cfg.regime_names,
        "start_time" => cfg.start_time,
        "initial_regime" => cfg.initial_regime,
        "switching_policy" => cfg.switching_policy,
        "min_switch_time" => cfg.min_switch_time,
        "pre_to_transition_threshold" => cfg.pre_to_transition_threshold,
        "transition_to_post_threshold" => cfg.transition_to_post_threshold,
        "transition_to_pre_threshold" => cfg.transition_to_pre_threshold,
        "post_to_transition_threshold" => cfg.post_to_transition_threshold,
        "allow_backward_switches" => cfg.allow_backward_switches,
        "allow_parent_fom_mismatch" => cfg.allow_parent_fom_mismatch,
        "safety" => cfg.safety,
        "growthfac" => cfg.growthfac,
    )
    isnothing(cfg.run_tag) || (d["run_tag"] = cfg.run_tag)
    isnothing(cfg.end_time) || (d["end_time"] = cfg.end_time)
    isnothing(cfg.pre_to_transition_time) || (d["pre_to_transition_time"] = cfg.pre_to_transition_time)
    isnothing(cfg.transition_to_post_time) || (d["transition_to_post_time"] = cfg.transition_to_post_time)
    return d
end

function dict_to_fom(d::AbstractDict)
    FOMConfig(
        viscosity = Float64(d["viscosity"]),
        conductivity = Float64(d["conductivity"]),
        gravity = Float64(d["gravity"]),
        gdir = Int(d["gdir"]),
        dodissipation = haskey(d, "dodissipation") ? Bool(d["dodissipation"]) : true,
        nx = Int(d["nx"]),
        ny = haskey(d, "ny") ? Int(d["ny"]) : Int(d["n"]),
        xlims = Tuple(Float64.(d["xlims"])),
        ylims = Tuple(Float64.(d["ylims"])),
        stretch_x = haskey(d, "stretch_x") ? _stretch_from_toml(d["stretch_x"]) : nothing,
        stretch_y = haskey(d, "stretch_y") ? _stretch_from_toml(d["stretch_y"]) : nothing,
        temp_hot = Float64(d["temp_hot"]),
        temp_cold = Float64(d["temp_cold"]),
        initial_temperature_profile = haskey(d, "initial_temperature_profile") ?
            String(d["initial_temperature_profile"]) : "uniform",
        perturbation_amplitude = Float64(d["perturbation_amplitude"]),
        perturbation_wavenumber = Float64(d["perturbation_wavenumber"]),
        tsim = Float64(d["tsim"]),
        cfl = Float64(d["cfl"]),
        nsnapshot = Int(d["nsnapshot"]),
    )
end

function dict_to_rom(d::AbstractDict)
    ROMConfig(
        run_tag = haskey(d, "run_tag") ? String(d["run_tag"]) : nothing,
        npod = Int(d["npod"]),
        rank_selection = haskey(d, "rank_selection") ? String(d["rank_selection"]) : "npod",
        basis_centered = haskey(d, "basis_centered") ? Bool(d["basis_centered"]) :
            lowercase(String(get(d, "basis_mode", "global"))) == "global",
        lambda = Float64(d["lambda"]),
        safety = Float64(d["safety"]),
        growthfac = Float64(d["growthfac"]),
        energy_fraction = haskey(d, "energy_fraction") ? Float64(d["energy_fraction"]) : 0.99,
        training_start_time = haskey(d, "training_start_time") ? Float64(d["training_start_time"]) : 0.0,
        training_end_time = haskey(d, "training_end_time") ? Float64(d["training_end_time"]) : nothing,
        use_intrusive_rom = haskey(d, "use_intrusive_rom") ? Bool(d["use_intrusive_rom"]) : false,
        use_exact_opinf = haskey(d, "use_exact_opinf") ? Bool(d["use_exact_opinf"]) : true,
        exact_opinf_force_explicit_euler = haskey(d, "exact_opinf_force_explicit_euler") ?
            Bool(d["exact_opinf_force_explicit_euler"]) : true,
        fallback_to_standard_opinf = haskey(d, "fallback_to_standard_opinf") ?
            Bool(d["fallback_to_standard_opinf"]) : true,
        include_constant_term = haskey(d, "include_constant_term") ? Bool(d["include_constant_term"]) : true,

        basis_mode = haskey(d, "basis_mode") ? String(d["basis_mode"]) : "global",
        basis_local_strategy = haskey(d, "basis_local_strategy") ? String(d["basis_local_strategy"]) : "S3",
        basis_partition_method = haskey(d, "basis_partition_method") ? String(d["basis_partition_method"]) : "uniform_x",
        basis_partition_nx = haskey(d, "basis_partition_nx") ? Int(d["basis_partition_nx"]) : 1,
        basis_partition_ny = haskey(d, "basis_partition_ny") ? Int(d["basis_partition_ny"]) : 1,
        basis_partition_xcuts = haskey(d, "basis_partition_xcuts") && !isnothing(d["basis_partition_xcuts"]) ?
            Float64.(d["basis_partition_xcuts"]) : nothing,
        basis_partition_ycuts = haskey(d, "basis_partition_ycuts") && !isnothing(d["basis_partition_ycuts"]) ?
            Float64.(d["basis_partition_ycuts"]) : nothing,
        basis_local_energy_fraction = haskey(d, "basis_local_energy_fraction") ? Float64(d["basis_local_energy_fraction"]) : 0.99,
        basis_local_ranks = haskey(d, "basis_local_ranks") && !isnothing(d["basis_local_ranks"]) ?
            Int.(d["basis_local_ranks"]) : nothing,
        basis_space_local_sparse = haskey(d, "basis_space_local_sparse") ? Bool(d["basis_space_local_sparse"]) : true,
    )
end

function dict_to_switched_rom(d::AbstractDict)
    cfg = SwitchedROMConfig(
        run_tag = haskey(d, "run_tag") ? String(d["run_tag"]) : nothing,
        local_rom_ids = haskey(d, "local_rom_ids") ? String.(d["local_rom_ids"]) : String[],
        regime_names = haskey(d, "regime_names") ? String.(d["regime_names"]) : ["pre", "transition", "post"],
        start_time = haskey(d, "start_time") ? Float64(d["start_time"]) : 0.0,
        end_time = haskey(d, "end_time") ? Float64(d["end_time"]) : nothing,
        initial_regime = haskey(d, "initial_regime") ? String(d["initial_regime"]) : "pre",
        switching_policy = haskey(d, "switching_policy") ? String(d["switching_policy"]) : "scheduled",
        min_switch_time = haskey(d, "min_switch_time") ? Float64(d["min_switch_time"]) : 20.0,
        pre_to_transition_threshold = haskey(d, "pre_to_transition_threshold") ?
            Float64(d["pre_to_transition_threshold"]) : 0.15,
        transition_to_post_threshold = haskey(d, "transition_to_post_threshold") ?
            Float64(d["transition_to_post_threshold"]) : -0.30,
        transition_to_pre_threshold = haskey(d, "transition_to_pre_threshold") ?
            Float64(d["transition_to_pre_threshold"]) : 0.20,
        post_to_transition_threshold = haskey(d, "post_to_transition_threshold") ?
            Float64(d["post_to_transition_threshold"]) : -0.20,
        allow_backward_switches = haskey(d, "allow_backward_switches") ?
            Bool(d["allow_backward_switches"]) : false,
        pre_to_transition_time = haskey(d, "pre_to_transition_time") ?
            Float64(d["pre_to_transition_time"]) : nothing,
        transition_to_post_time = haskey(d, "transition_to_post_time") ?
            Float64(d["transition_to_post_time"]) : nothing,
        allow_parent_fom_mismatch = haskey(d, "allow_parent_fom_mismatch") ?
            Bool(d["allow_parent_fom_mismatch"]) : false,
        safety = haskey(d, "safety") ? Float64(d["safety"]) : 0.03,
        growthfac = haskey(d, "growthfac") ? Float64(d["growthfac"]) : 1.05,
    )
    return validate_switched_rom_config(cfg)
end

function save_fom_config(path::AbstractString, cfg::FOMConfig)
    open(path, "w") do io
        TOML.print(io, fom_to_dict(cfg))
    end
end

function save_rom_config(path::AbstractString, cfg::ROMConfig)
    open(path, "w") do io
        TOML.print(io, rom_to_dict(cfg))
    end
end

function load_switched_rom_config(path::AbstractString)
    dict_to_switched_rom(TOML.parsefile(path))
end

function save_switched_rom_config(path::AbstractString, cfg::SwitchedROMConfig)
    open(path, "w") do io
        TOML.print(io, switched_rom_to_dict(cfg); sorted = true)
    end
end

load_fom_config(path::AbstractString) = dict_to_fom(TOML.parsefile(path))
load_rom_config(path::AbstractString) = dict_to_rom(TOML.parsefile(path))
