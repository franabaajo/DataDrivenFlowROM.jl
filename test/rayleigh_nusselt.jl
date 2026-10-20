@testitem "Rayleigh-Benard Nusselt conduction sanity" begin
    using IncompressibleNavierStokes

    include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "Config.jl"))
    include(joinpath(@__DIR__, "..", "src", "rayleighbenard", "ThermalDiagnostics.jl"))

    for stretch in (nothing, 1.2)
        cfg = FOMConfig(
            nx = 24,
            ny = 20,
            xlims = (0.0, 2.0),
            ylims = (-0.5, 1.5),
            stretch_x = stretch,
            stretch_y = stretch,
            temp_hot = 3.0,
            temp_cold = 1.0,
            conductivity = 0.2,
            dodissipation = false,
        )
        x = isnothing(stretch) ?
            range(cfg.xlims..., cfg.nx + 1) :
            tanh_grid(cfg.xlims..., cfg.nx, stretch)
        y = isnothing(stretch) ?
            range(cfg.ylims..., cfg.ny + 1) :
            tanh_grid(cfg.ylims..., cfg.ny, stretch)
        setup = Setup(
            x = (x, y),
            boundary_conditions = (
                u = ((DirichletBC(), DirichletBC()), (DirichletBC(), DirichletBC())),
                temp = (
                    (SymmetricBC(), SymmetricBC()),
                    (DirichletBC(cfg.temp_hot), DirichletBC(cfg.temp_cold)),
                ),
            ),
        )

        H = cfg.ylims[2] - cfg.ylims[1]
        ΔT = cfg.temp_hot - cfg.temp_cold
        temp = temperaturefield(
            setup,
            (x, y) -> cfg.temp_hot - ΔT * (y - cfg.ylims[1]) / H,
        )
        u = vectorfield(setup)
        fill!(u, 0)
        u = apply_bc_u(u, 0.0, setup)

        plate = plate_nusselt_numbers(temp, setup, cfg)
        @test plate.bottom ≈ 1.0 atol = 5e-13
        @test plate.top ≈ 1.0 atol = 5e-13
        @test volume_nusselt_number(u, temp, setup, cfg) ≈ 1.0 atol = 5e-13

        budget = temperature_heat_budget(u, temp, zero(temp), setup, cfg)
        @test budget.rate ≈ 0.0 atol = 5e-13
        @test budget.boundary ≈ 0.0 atol = 5e-13
        @test budget.source ≈ 0.0 atol = 5e-13
        @test budget.residual ≈ 0.0 atol = 5e-13
    end
end
