@testitem "Rayleigh-Benard ROM training window" begin
    include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Config.jl"))
    include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Paths.jl"))

    fomcfg = FOMConfig(tsim = 200.0, nsnapshot = 500)

    window = rom_training_window(fomcfg, 501, 20.0)
    @test window.snapshot_indices == 51:501
    @test window.actual_start_time == 20.0
    @test window.duration == 180.0
    @test length(window.times) == 451

    bounded = rom_training_window(fomcfg, 501, 20.0, 80.0)
    @test bounded.snapshot_indices == 51:201
    @test bounded.actual_start_time == 20.0
    @test bounded.end_time == 80.0
    @test bounded.duration == 60.0
    @test length(bounded.times) == 151

    rounded = rom_training_window(fomcfg, 501, 20.1)
    @test first(rounded.snapshot_indices) == 52
    @test rounded.actual_start_time == 20.4

    rounded_end = rom_training_window(fomcfg, 501, 20.0, 80.1)
    @test last(rounded_end.snapshot_indices) == 201
    @test rounded_end.end_time == 80.0

    @test_throws Exception rom_training_window(fomcfg, 501, -1.0)
    @test_throws Exception rom_training_window(fomcfg, 501, 200.0)
    @test_throws Exception rom_training_window(fomcfg, 501, 20.0, 20.0)
    @test_throws Exception rom_training_window(fomcfg, 501, 20.0, 201.0)
    @test_throws Exception rom_training_window(fomcfg, 500, 20.0)

    parent = "fom_test_12345678"
    cfg_full = ROMConfig()
    cfg_windowed = ROMConfig(training_start_time = 20.0)
    cfg_bounded = ROMConfig(training_start_time = 20.0, training_end_time = 80.0)
    @test rom_id(parent, cfg_full) != rom_id(parent, cfg_windowed)
    @test occursin("_ts20p0_", rom_id(parent, cfg_windowed))
    @test occursin("_te80p0_", rom_id(parent, cfg_bounded))
    @test dict_to_rom(rom_to_dict(cfg_bounded)).training_end_time == 80.0

    switchcfg = SwitchedROMConfig(
        local_rom_ids = ["rom_pre", "rom_transition", "rom_post"],
        start_time = 20.0,
        end_time = 80.0,
        pre_to_transition_time = 35.0,
        transition_to_post_time = 60.0,
        run_tag = "state_switch",
    )
    @test dict_to_switched_rom(switched_rom_to_dict(switchcfg)).end_time == 80.0
    @test occursin("_switched_lz_", switched_rom_id(parent, switchcfg))
    @test occursin("_te80p0_", switched_rom_id(parent, switchcfg))
    @test_throws Exception validate_switched_rom_config(SwitchedROMConfig(local_rom_ids = ["only_one"]))

    @test cfg_full.basis_centered == true
    cfg_raw = ROMConfig(basis_centered = false)
    @test rom_id(parent, cfg_full) != rom_id(parent, cfg_raw)

    raw_energy_rank = ROMConfig(rank_selection = "energy_fraction", basis_centered = false)
    @test dict_to_rom(rom_to_dict(raw_energy_rank)).basis_centered == false

    centered_fixed_rank = ROMConfig(rank_selection = "npod", basis_centered = true)
    @test dict_to_rom(rom_to_dict(centered_fixed_rank)).basis_centered == true

    old_dict = rom_to_dict(cfg_windowed)
    delete!(old_dict, "training_start_time")
    @test dict_to_rom(old_dict).training_start_time == 0.0

    legacy_rom_dict = rom_to_dict(cfg_raw)
    delete!(legacy_rom_dict, "basis_centered")
    @test dict_to_rom(legacy_rom_dict).basis_centered == true

    legacy_space_local_dict = rom_to_dict(ROMConfig(basis_mode = "space_local", basis_centered = false))
    delete!(legacy_space_local_dict, "basis_centered")
    @test dict_to_rom(legacy_space_local_dict).basis_centered == false
end
