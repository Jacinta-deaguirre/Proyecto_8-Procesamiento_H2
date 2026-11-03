module HydrogenAdsorption
export generate_profiles_plot, heat_equation

using SearchSortedNearest
using Reexport
using Plots

include("coefficient_matrix.jl")
export coefficient_matrix

include("equations_of_state.jl")
export ideal_gas_equation, cyl_avg

include("adsorption_systems.jl")
export IsothermParameters, MDAParameters, DAParameters, MaterialProperties, GeometricParameters, OperationalParameters, AdsorptionParameters,
    adsorption_isotherm, isosteric_heat_of_adsorption, MDA_adsorption!, dae_setup

"""
Energy balance of the packed bed (enthalpy form), per unit bed volume:

    (ρ_s(1-ε_b) c_ps + ρ_g ε_b c_pg + n_a ρ_s(1-ε_b) M_H2 c_vg) ∂T/∂t
        = k_eff/r ∂/∂r(r ∂T/∂r)                 (Q_conductive, discretised as A*T - b)
        + ρ_s(1-ε_b) ΔH ∂n_a/∂t                 (Q_adsorptive)
        + max(ṁ,0)/V c_pg (T_H2 - T)            (Q_mass_flow: inflow enthalpy; the outflow
                                                 leaves at the local bed temperature, so the
                                                 term vanishes during discharge)
        + ε_b dP/dt                             (Q_compression: compression/expansion work of
                                                 the gas phase, uniform because P is uniform)

The compression term is evaluated from the differential variables through the ideal-gas
closure P = ρ̄ R_H2 T_H, with T_H the harmonic volume mean of T (see `MDA_adsorption!`).

# Inputs:
- `materialProps`: Material properties struct containing ρₛ, cₛ, mₛ, ε_b, cₚ, cᵥ, M_H2, R
- `geometricParams`: Geometric parameters struct containing A, b, V, n_r, r_span, R_T
- `operationalParams`: Operational parameters (m_in, T_gas)
- `u`: State variables vector (T, n_a, ρ_avg, P, ρ, T_wall)
- `du`: Time derivative of the state variables
- `dH`: Isosteric heat of adsorption vector / J/mol

# Outputs:
- `dT`: Rate of change of temperature vector
"""
function heat_equation(materialProps::MaterialProperties, geometricParams::GeometricParameters, operationalParams::OperationalParameters, u, du, dH)
    # Material properties
    cₛ = materialProps.cₛ
    ε_b = materialProps.ε_b
    cₚ = materialProps.cₚ
    cᵥ = materialProps.cᵥ
    mₛ = materialProps.mₛ
    M_H2 = materialProps.M_H2
    R = materialProps.R

    # Geometric parameters
    A = geometricParams.A
    b = geometricParams.b
    V = geometricParams.V
    n_r = geometricParams.n_r
    r_span = geometricParams.r_span
    R_T = geometricParams.R_T

    # Operational parameters
    m_in = operationalParams.m_in
    T_H2 = operationalParams.T_gas

    # Unpack state variables
    T = u[1:n_r]
    nₐ = u[n_r+1:2*n_r]
    ρ_avg = u[2*n_r+1]
    ρ = u[2*n_r+3:3*n_r+2]
    dTdt = du[1:n_r]
    dρ_avg = du[2*n_r+1]

    # Effective volumetric heat capacity (solid + gas at constant pressure + adsorbed phase)
    dT_cte_term = mₛ .* cₛ ./ V .+ ρ .* ε_b .* cₚ .+ nₐ .* mₛ .* M_H2 ./ V .* cᵥ

    Q_conductive = A * T .- b                          # Radial conduction
    Q_adsorptive = mₛ / V .* dH .* du[n_r+1:2*n_r]     # Heat of adsorption
    Q_mass_flow = max(m_in, 0.0) / V .* cₚ .* (T_H2 .- T) # Inlet enthalpy (charging only)

    # Compression work ε_b dP/dt, with P = ρ̄ R_H2 T_H and T_H = 1/<1/T>
    T_H = 1 / cyl_avg(1 ./ T, r_span, R_T)
    dT_H = T_H^2 * cyl_avg(dTdt ./ T .^ 2, r_span, R_T)
    dPdt = R / M_H2 * (dρ_avg * T_H + ρ_avg * dT_H)
    Q_compression = ε_b * dPdt

    dT = (Q_conductive .+ Q_adsorptive .+ Q_mass_flow .+ Q_compression) ./ dT_cte_term

    return dT
end

"""
Function to generate the profiles plot for temperature, adsorption, and gas density
at specified time steps.

Inputs:
- `t`: Vector of time values / s
- `r_span`: Radial span of the tank / Vector of radial positions at each node / m
- `T_profile`: Matrix of temperature profiles at each time step / K
- `Adsorption_profile`: Matrix of adsorption profiles at each time step / mol/kg_ads
- `P_profile`: Vector of pressure profiles at each time step / Pa
- `Density_profile`: Matrix of gas density profiles at each time step / kg/m³
- `time_step`: Time step interval for plotting / s
- `colour_map`: Colour map for the plots

Outputs:
- A 2x2 grid plot with temperature, adsorption, gas density, and pressure profiles.
"""
function generate_profiles_plot(t, r_span, T_profile, Adsorption_profile, P_profile, Density_profile, time_step, colour_map)
    plt = plot(layout=(2, 2), palette=colour_map) # This layout displays in a 2x2 grid

    idx = [] # Empty array to store the indexes of the time steps

    for i in range(0, stop=t[end], step=time_step)
        push!(idx, searchsortednearest(t, i))
    end

    unique!(idx) # Remove duplicates

    for i in range(1, stop=length(idx), step=1)
        plot!(plt[1], r_span, T_profile[idx[i], :], xlabel="Radius (m)", ylabel="Temperature (K)", label="t = $(round(t[idx[i]], digits=2)) s", legend=:outerbottomright)
        plot!(plt[2], r_span, Adsorption_profile[idx[i], :], xlabel="Radius (m)", ylabel="q_a / mol/kg_ads", label="t = $(round(t[idx[i]], digits=2)) s", legend=false)
        plot!(plt[4], r_span, Density_profile[idx[i], :], xlabel="Radius (m)", ylabel="Gas density (kg/m³)", label="t = $(round(t[idx[i]], digits=2)) s", legend=false)
    end

    plot!(plt[3], t, P_profile .* 1e-6, xlabel="Time (s)", ylabel="Pressure (MPa)", label="Pressure vs Time", legend=false, color=:black)
    plot!(size=(1600, 760))
    plot!(margin=Plots.cm)
    display(plt)
end

end
