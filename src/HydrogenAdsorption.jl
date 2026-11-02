module HydrogenAdsorption
export generate_profiles_plot, heat_equation

using SearchSortedNearest
using Reexport
using Plots

include("coefficient_matrix.jl")
export coefficient_matrix

include("adsorption_systems.jl")
export IsothermParameters, MDAParameters, DAParameters, MaterialProperties, GeometricParameters, OperationalParameters, AdsorptionParameters,
    adsorption_isotherm, isosteric_heat_of_adsorption, MDA_adsorption!, dae_setup

include("equations_of_state.jl")
export ideal_gas_equation

"""
Function to compute the heat equations for the adsorption system.

# Inputs:
- `materialProps`: Material properties struct containing ρₛ, cₛ, ε_b, cₚ
- `geometricParams`: Geometric parameters struct containing A, V, n_r
- `u`: State variables vector containing temperature T, average density ρ_avg
- `du`: Derivative of state variables vector
- `dH`: Isosteric heat of adsorption vector

# Outputs:
- `dT`: Rate of change of temperature vector
"""
function heat_equation(materialProps::MaterialProperties, geometricParams::GeometricParameters, operationalParams::OperationalParameters, u, du, dH)
    # Unpack parameters
    # Material properties
    ρₛ = materialProps.ρₛ
    cₛ = materialProps.cₛ
    ε_b = materialProps.ε_b
    cₚ = materialProps.cₚ
    cᵥ = materialProps.cᵥ
    mₛ = materialProps.mₛ
    M_H2 = materialProps.M_H2

    # Geometric parameters
    A = geometricParams.A
    V = geometricParams.V
    n_r = geometricParams.n_r
    b = geometricParams.b

    # Operational parameters
    m_in = operationalParams.m_in
    T_H2 = operationalParams.T_gas

    # Unpack state variables
    T = u[1:n_r]
    nₐ = u[n_r+1:2*n_r]
    ρ_avg = u[2*n_r+1]
    P = u[2*n_r+2]
    ρ = u[2*n_r+3:3*n_r+2]
    T_wall = u[3*n_r+3] # Tank wall temperature

    # Prueba dividiendo por dV
    dT_cte_term = 1 .* (mₛ .* cₛ ./ V .+ ρ .* ε_b .* cᵥ .+ nₐ .* mₛ .* M_H2 ./ V .* cᵥ)

    # Q_conductive = 2 .* π .* L .* dr ./ R_T .* A * T # Conductive heat transfer term
    Q_conductive = A * T .- b # Conductive heat transfer term
    Q_adsorptive = mₛ / V .* dH .* du[n_r+1:2*n_r] # Adsorptive heat transfer term
    # Q_mass_flow = m_in .* cₚ .* (T_H2 .- T) # Mass flow heat transfer term
    Q_mass_flow = m_in / V .* cₚ .* (T_H2 .- T) # Mass flow heat transfer term
    # avg_Q_mass_flow = mean(Q_mass_flow)
    #println("Q_mass_flow: $(avg_Q_mass_flow) W/m³")
    #println("Q_conductive: $(mean(Q_conductive)) W/m³")
    #println("Q_adsorptive: $(mean(Q_adsorptive)) W/m³") 

    dT = 1 ./ (dT_cte_term) .* (1 * Q_conductive .+ 1 * Q_adsorptive .+ 1 * Q_mass_flow)

    #println("Denominador: $dT_cte_term")
    #println("Numerador: $(A * T + mₛ .* dH .* du[n_r+1:2*n_r] .+ m_in .* cₚ .* (T_H2 .- T))")

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
        # This loop is used to find the indexes of the time steps nearest to 
        # t = 100s, 200s, 300s, etc.
        # serachsortednearest is used to find the index of the nearest value
        # push! is used to append the index to the idx array
        push!(idx, searchsortednearest(t, i))
    end

    unique!(idx) # Remove duplicates

    for i in range(1, stop=length(idx), step=1)
        # This loop is used to plot the temperature, adsorption, and gas density
        # property[idx[i]] is used to access the values of the properties at 
        # t = 100s, 200s, 300s, etc.
        # t[idx[i]] the time at the specific index

        # Temperature plot
        plot!(plt[1], r_span, T_profile[idx[i], :], xlabel="Radius (m)", ylabel="Temperature (K)", label="t = $(round(t[idx[i]], digits=2)) s", legend=:outerbottomright)

        # Adsorption plot
        plot!(plt[2], r_span, Adsorption_profile[idx[i], :], xlabel="Radius (m)", ylabel="q_a / mol/kg_ads", label="t = $(round(t[idx[i]], digits=2)) s", legend=false)

        # Gas density plot
        plot!(plt[4], r_span, Density_profile[idx[i], :], xlabel="Radius (m)", ylabel="Gas density (kg/m³)", label="t = $(round(t[idx[i]], digits=2)) s", legend=false)
    end

    # Pressure plot
    plot!(plt[3], t, P_profile .* 1e-6, xlabel="Time (s)", ylabel="Pressure (MPa)", label="Pressure vs Time", legend=false, color=:black)

    # Plot size for better visualization
    plot!(size=(1600, 760))
    plot!(margin=Plots.cm)

    # Display the plot
    display(plt)
    # savefig(plt, "1D_cylindrical/figs/1D_cylindrical_model_$(time_step)s.svg")
end

end
