if !isdefined(@__MODULE__, :_RAYLEIGH_THERMAL_DIAGNOSTICS_LOADED)
const _RAYLEIGH_THERMAL_DIAGNOSTICS_LOADED = true

import IncompressibleNavierStokes as NS

function _rayleigh_thermal_scales(cfg::FOMConfig)
    Lx = cfg.xlims[2] - cfg.xlims[1]
    H = cfg.ylims[2] - cfg.ylims[1]
    ΔT = cfg.temp_hot - cfg.temp_cold
    Lx > 0 || error("Rayleigh-Benard domain width must be positive, got $Lx")
    H > 0 || error("Rayleigh-Benard domain height must be positive, got $H")
    ΔT != 0 || error("Rayleigh-Benard plate temperature difference must be nonzero")
    return (; Lx, H, ΔT, volume = Lx * H)
end

"""
Return bottom and top plate Nusselt numbers using the wall value and adjacent
interior cell-center value:

`Nu = -(H/ΔT) ⟨∂T/∂y⟩ₓ`.

Both values are positive and equal to one for the conductive state
`T = T_hot - ΔT * (y-ymin)/H`.
"""
function plate_nusselt_numbers(temp, setup, cfg::FOMConfig)
    setup.dimension() == 2 || error("plate_nusselt_numbers currently supports 2D setups")
    (; Lx, H, ΔT) = _rayleigh_thermal_scales(cfg)

    temp_cpu = Array(temp)
    Δx = Array(setup.Δ[1])
    Δyu = Array(setup.Δu[2])
    ix, iy = setup.Ip.indices

    jbottom = first(iy) - 1
    jfirst = first(iy)
    jlast = last(iy)
    jtop = last(iy) + 1

    grad_bottom = (@view(temp_cpu[ix, jfirst])) .- (@view(temp_cpu[ix, jbottom]))
    grad_bottom ./= Δyu[jbottom]
    grad_top = (@view(temp_cpu[ix, jtop])) .- (@view(temp_cpu[ix, jlast]))
    grad_top ./= Δyu[jlast]

    scale = -H / ΔT / Lx
    bottom = scale * sum(@view(Δx[ix]) .* grad_bottom)
    top = scale * sum(@view(Δx[ix]) .* grad_top)
    return (; bottom, top)
end

"""Volume average of a cell-centered scalar over the physical domain."""
function scalar_volume_average(field, setup, cfg::FOMConfig)
    setup.dimension() == 2 || error("scalar_volume_average currently supports 2D setups")
    (; volume) = _rayleigh_thermal_scales(cfg)

    field_cpu = Array(field)
    Δx = Array(setup.Δ[1])
    Δy = Array(setup.Δ[2])
    ix, iy = setup.Ip.indices
    weights = reshape(@view(Δx[ix]), :, 1) .* reshape(@view(Δy[iy]), 1, :)
    return sum(@view(field_cpu[ix, iy]) .* weights) / volume
end

"""
Return the volume/global Nusselt number

`Nu_V = H/(κ ΔT) ⟨vT - κ ∂T/∂y⟩_V`.

Since this workflow evolves the full temperature `T`, the conductive part
integrates to one and this formula is equivalent to
`1 + H/(κ ΔT) ⟨vT⟩_V`.
"""
function volume_nusselt_number(u, temp, setup, cfg::FOMConfig)
    setup.dimension() == 2 || error("volume_nusselt_number currently supports 2D setups")
    cfg.conductivity > 0 || error("conductivity must be positive")
    (; Lx, H, ΔT) = _rayleigh_thermal_scales(cfg)

    u_cpu = Array(u)
    temp_cpu = Array(temp)
    Δx = Array(setup.Δ[1])
    Δy = Array(setup.Δ[2])
    Δyu = Array(setup.Δu[2])
    ix, iy = setup.Ip.indices

    integrated_flux = 0.0
    for j in (first(iy)-1):last(iy)
        denominator = Δy[j] + Δy[j+1]
        temp_face =
            (Δy[j+1] .* @view(temp_cpu[ix, j]) .+
             Δy[j] .* @view(temp_cpu[ix, j+1])) ./ denominator
        grad_face = ((@view(temp_cpu[ix, j+1])) .- (@view(temp_cpu[ix, j]))) ./ Δyu[j]
        flux_face =
            sum(@view(Δx[ix]) .* ((@view(u_cpu[ix, j, 2])) .* temp_face .-
                                  cfg.conductivity .* grad_face)) / Lx
        integrated_flux += flux_face * Δyu[j]
    end

    mean_vertical_flux = integrated_flux / H
    return H * mean_vertical_flux / (cfg.conductivity * ΔT)
end

"""Volume-averaged viscous-heating source used by the configured temperature equation."""
function temperature_source_average(u, setup, cfg::FOMConfig)
    cfg.dodissipation || return 0.0
    source = NS.dissipation(u, setup, cfg.viscosity)
    return scalar_volume_average(source, setup, cfg)
end

"""
Evaluate the discrete mean-temperature heat budget:

`d⟨T⟩/dt = κ ΔT/H² (Nu_bottom - Nu_top) + ⟨S_T⟩`.
"""
function temperature_heat_budget(u, temp, dtemp, setup, cfg::FOMConfig)
    (; H, ΔT) = _rayleigh_thermal_scales(cfg)
    plate = plate_nusselt_numbers(temp, setup, cfg)
    rate = scalar_volume_average(dtemp, setup, cfg)
    source = temperature_source_average(u, setup, cfg)
    boundary = cfg.conductivity * ΔT / H^2 * (plate.bottom - plate.top)
    return (;
        plate...,
        rate,
        boundary,
        source,
        residual = rate - boundary - source,
    )
end

end
