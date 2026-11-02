using Statistics

"""
Abstract type for adsorption isotherm parameters.
This type serves as a base for different isotherm models, such as the Modified Dubinin-Astakov (MDA) and Dubinin-Astakov (DA) models.
It allows for multiple dispatch in functions that handle different isotherm models, enabling the use of different parameter sets without changing the function signatures.
"""
abstract type IsothermParameters end

"""
Struct used to hold the adsorption isotherm parameters for the Modified Dubinin-Astakov model.
"""
struct MDAParameters <: IsothermParameters
    n₀::Float64     # Limit adsorption / mol/kg
    p₀::Float64     # Saturation pressure / Pa
    α::Float64      # Enthalpic Factor / J/mol
    β::Float64      # Entropic Factor / J/mol K
    m::Float64      # Exponent in the isotherm equation
end

"""
Struct used to hold the adsorption isotherm parameters for the Dubinin-Astakov model.
"""
struct DAParameters <: IsothermParameters
    P_limit::Float64    # Limit pressure / Pa
    # Limiting adsorption parameters
    ψ::Float64          # mmol/g
    β::Float64          # mol/(kg K)
    # Denominator parameters
    κ::Float64          # J mol-1
    γ::Float64          # J mol-1 K-1
    m::Float64          # Exponent in the isotherm equation
end

"""
Struct to hold the material properties of the hydrogen tank and the adsorbent material.
This struct includes properties for both the activated carbon adsorbent and the hydrogen gas, as well as effective properties that are derived from these materials.
"""
struct MaterialProperties
    # Activated carbon properties
    ρₛ::Float64        # Density of the solid phase / kg/m³
    cₛ::Float64        # Specific heat capacity of the solid phase / J/(kg·K)
    mₛ::Float64        # Mass of the solid phase / kg
    kₛ::Float64        # Thermal conductivity of the solid phase / W/(m·K)
    ε_b::Float64       # Bulk porosity of the material

    # Hydrogen properties
    cₚ::Float64        # Specific heat capacity of hydrogen / J/(kg·K)
    cᵥ::Float64        # Specific heat capacity at constant volume / J/(kg·K)
    M_H2::Float64      # Molar mass of hydrogen / kg/mol
    R::Float64         # Universal gas constant / J/(mol·K)
    k_g::Float64        # Thermal conductivity of hydrogen / W/(m·K)

    # Effective properties
    k_eff::Float64     # Effective thermal conductivity / W/(m·K)
end

"""
Struct to hold the geometric parameters of the hydrogen tank.
This struct includes parameters such as the number of radial nodes, radial step size, volume of the tank, coefficient matrix for spatial discretization, radial span of the tank, and the radius of the tank.
It also includes the steel thermophysical properties.
"""
struct GeometricParameters
    n_r::Int                # Number of radial nodes
    dr::Float64             # Radial step size / m
    V::Float64              # Volume of the tank / m³
    L::Float64              # Length of the tank / m
    A::Matrix{Float64}      # Coefficient matrix for the spatial discretization
    b::Vector{Float64}      # Right-hand side vector for the spatial discretization
    r_span::Vector{Float64} # Radial span of the tank / m
    R_T::Float64            # Radius of the tank / m
    e::Float64              # Thickness of the tank wall / m
    Ao::Float64             # Outer surface area of the tank / m²
    Ai::Float64             # Inner surface area of the tank / m²

    # Wall properties
    c_wall::Float64         # Specific heat capacity of the tank wall / J/(kg·K)
    m_wall::Float64         # Mass of the tank wall / kg
    k_wall::Float64         # Thermal conductivity of the tank wall / W/(m·K)
end

"""
Struct to hold the operational parameters of the hydrogen tank system.
This struct includes parameters such as the heat transfer coefficient and the mass flow rate of hydrogen into the tank.
"""
struct OperationalParameters
    U::Float64      # Heat transfer coefficient / W/(m²·K)
    T_air::Float64  # Ambient temperature / K
    m_in            # Mass flow rate of hydrogen into the tank / kg/s
    T_gas::Float64  # Temperature of the incoming hydrogen gas / K
end

struct AdsorptionParameters
    isotherm::IsothermParameters
    material::MaterialProperties
    geometric::GeometricParameters
    operational::OperationalParameters
end

"""
===MODIFIED DUBININ-ASTAKOV DISPATCH===

Function to compute the adsorption isotherm using the Modified Dubinin-Astakov model.
This function calculates the amount of hydrogen adsorbed based on the pressure and temperature, using the parameters defined in the `MDAParameters` struct.

Inputs:
- `params`: An instance of `MDAParameters` containing the isotherm parameters.
- `P`: Pressure in the tank / Pa
- `T`: Temperature in the tank / Vector of temperatures at each radial node / K

Outputs:
- `nₐ`: Amount of hydrogen adsorbed at each radial node / mol/kg
"""
function adsorption_isotherm(params::MDAParameters, P, T)
    # Unpacking parameters
    n₀ = params.n₀      # Limit adsorption / mol/kg
    p₀ = params.p₀      # Saturation pressure / Pa
    α = params.α        # Enthalpic Factor / J/mol
    β = params.β        # Entropic Factor / J/mol K
    m = params.m        # Exponent in the isotherm equation

    # Other constants
    R = 8.314 # Universal gas constant / J/(mol·K)

    # Modified Dubinin-Astakov isotherm equation
    nₐ = n₀ .* exp.(-(R .* T ./ (α .+ β .* T)) .^ m .* (log.(p₀ ./ P)) .^ m)

    return nₐ
end

"""
===DUBININ-ASTAKOV DISPATCH===

Function to compute the adsorption isotherm using the Dubinin-Astakov model.
This function calculates the amount of hydrogen adsorbed based on the pressure and temperature, using the parameters defined in the `DAParameters` struct.

Inputs:
- `params`: An instance of `DAParameters` containing the isotherm parameters.
- `P`: Pressure in the tank / Pa
- `T`: Temperature in the tank / Vector of temperatures at each radial node / K

Outputs:
- `nₐ`: Amount of hydrogen adsorbed at each radial node / mol/kg
"""
function adsorption_isotherm(params::DAParameters, P, T)
    # Unpacking parameters
    P_lim = params.P_limit  # Limit pressure / Pa
    ψ = params.ψ            # Limiting adsorption / mmol/g
    β = params.β
    κ = params.κ
    γ = params.γ
    m = params.m            # Exponent in the isotherm equation

    # Other constants
    R = 8.314 # Universal gas constant / J/(mol·K)

    # Calculate Parameters
    n_0 = ψ .+ β .* T
    E = κ .+ γ .* T
    A = R .* T .* log.(P_lim ./ P)

    # Dubinin-Astakov isotherm equation
    nₐ = n_0 .* exp.(-(A ./ E) .^ m)
    return nₐ
end

"""
    Function to compute the isosteric heat of adsorption based on the isotherm parameters.
    This function implements the Dubinin Astakov equation.
    
    Inputs:
    - `params`: An instance of `MDAParameters` containing the isotherm parameters.
    - `P`: Pressure in the tank / Pa
    - `T`: Temperature in the tank / K
    
    Outputs:
    - `dH`: Isosteric heat of adsorption / J/mol
    """
function isosteric_heat_of_adsorption(params::MDAParameters, P, T)
    # Unpacking parameters
    n₀ = params.n₀      # Limit adsorption / mol/kg
    α = params.α        # Enthalpic Factor / J/mol
    m = params.m        # Exponent in the isotherm equation

    # Calculate the isosteric heat of adsorption
    dH = α .* ((log.(n₀ ./ adsorption_isotherm(params, P, T))) .^ (1 / m))

    return dH
end

"""
System of differential equations for the adsorption process in a hydrogen tank. 
This function computes the time derivative of the state variables, which include temperature, adsorption amount, average density, pressure, and local densities.
It uses the method of lines to discretize the spatial domain and solve the system of equations.

This system in particular uses the Modified Dubinin-Astakov isotherm for adsorption, which is suitable for hydrogen storage applications.

Inputs
- `out`: Output array
- `du`: Array containing the current time derivatives of the state variables.
- `u`: Current state of the system, including temperature, adsorption amount, average density, pressure, and local densities.
- `p`: Parameters for the system
- `t`: Current time (not used in this function, but included for compatibility with ODE solvers).

Outputs: The function modifies the `out` array in place

Description of variables:
- `T`: Array of dimension `n_r` representing the temperature at each radial node
- `nₐ`: Array of dimension `n_r` representing the adsorption amount at each radial node
- `ρ_avg`: Average density of hydrogen in the tank
- `P`: Pressure in the tank
- `ρ`: Array of dimension `n_r` representing the local densities of hydrogen at each radial node
"""
function MDA_adsorption!(out, du, u, p, t)
    # Unpacking variables
    T = u[1:n_r]
    nₐ = u[n_r+1:2*n_r]
    ρ_avg = u[2*n_r+1]
    P = u[2*n_r+2]
    ρ = u[2*n_r+3:end]

    dH = α .* ((log.(n₀ ./ nₐ)) .^ (1 / m))

    # Heat diffusion, implicit form
    # f(x) = 0 = dT/dt - A*T
    # A is the discretisation matrix
    # b should be zero as it is recalculated in 63.
    out[1:n_r] .= du[1:n_r] .- 1 / (ρₛ * cₛ * (1 - ε_b) + ρ_avg * cₚ * ε_b) * (A * T .+ dH .* du[n_r+1:2*n_r] ./ V .+ du[2*n_r+2])

    # Neumann BC time derivative for the tank centre 
    out[1] = du[1] - (4 * du[2] - du[3]) / 3

    # Robin BC time derivative to apply method of lines
    out[n_r] = du[n_r] - (4 * du[n_r-1] - du[n_r-2]) / (3 + 2 * U * dr / k_eff)

    # Adsorption equation
    out[n_r+1:2*n_r] .= nₐ .- n₀ .* exp.(-(R .* T ./ (α .+ β .* T)) .^ m .* (log.(p₀ / P)) .^ m)

    # Macroscopic mass balance 
    # mean(n_a .* r_span) / R computes the average adsorption of H2 
    out[2*n_r+1] = du[2*n_r+1] - (m_in / (V * ε_b) - ρₛ * (1 - ε_b) * M_H2 / ε_b * mean(du[n_r+1:2*n_r] .* r_span) / R_T)

    # Ideal gas pressure equation
    # mean(T .* r_span) / R_T computes the pressure inside the tank
    # out[2*n_r+2] = P - ρ_avg * R .* (mean(T .* r_span) / R_T) / M_H2
    # out[2*n_r+2] = du[2*n_r+2] - dPdt_fixed
    out[2*n_r+2] = du[2*n_r+2] - R / M_H2 * (du[2*n_r+1] * (mean(T .* r_span) / R_T) + ρ_avg * (mean(du[1:n_r] .* r_span) / R_T))

    # Pressure equation used to compute ρ in each node
    out[2*n_r+3:end] = P .- ρ .* R .* T ./ M_H2
end

function dae_setup(isotermParams::IsothermParameters,
    materialProps::MaterialProperties,
    geometricParams::GeometricParameters,
    operationalParams::OperationalParameters,
    Pᵢ::Float64,
    T_0::Float64,
    R=8.314)

    # Unpack parameters
    # Relevant material properties
    M_H2, k_eff, ε_b, ρₛ = materialProps.M_H2, materialProps.k_eff, materialProps.ε_b, materialProps.ρₛ

    # Relevant material parameters
    U, T_air, m_in = operationalParams.U, operationalParams.T_air, operationalParams.m_in

    # Relevant geometric parameters
    n_r, dr, V, R_T, r_span = geometricParams.n_r, geometricParams.dr, geometricParams.V, geometricParams.R_T, geometricParams.r_span
    e = geometricParams.e
    k_wall = geometricParams.k_wall


    # Calculate initial temperature vector
    Tᵢ = ones(n_r) * T_0

    # Calculate initial adsorption vector
    nₐᵢ = adsorption_isotherm(isotermParams, Pᵢ, Tᵢ)

    # Calculate initial density vector
    ρᵢ = ideal_gas_equation(Tᵢ, R, M_H2, P=Pᵢ)

    # Initial state vector
    u₀ = vcat(Tᵢ, nₐᵢ, ρᵢ[1], Pᵢ, ρᵢ, Tᵢ[end]) # Last term is tank wall temperature

    # Update Boundary Conditions
    # Robin BC for tank exterior
    u₀[n_r] = (4 * u₀[n_r-1] - u₀[n_r-2] + 4 * k_wall * dr / k_eff / e * T_0) / (3 + 4 * k_wall * dr / k_eff / e)

    # Update adsorption amount at the last node
    u₀[2*n_r] = adsorption_isotherm(isotermParams, Pᵢ, T_0)

    # Update density at the last node
    u₀[3*n_r+2] = ideal_gas_equation(u₀[n_r], R, M_H2, P=Pᵢ)

    # Initial time derivative vector
    du₀ = ones(length(u₀)) * 1e-5

    # Initial time derivative for the last node
    du₀[n_r] = (4 * du₀[n_r-1] - du₀[n_r-2]) / (3 + 2 * U * dr / k_eff)

    # Initial average density derivative for the last node
    du₀[2*n_r+1] = m_in / V / ε_b - ρₛ * (1 - ε_b) * M_H2 / ε_b * mean(du₀[n_r+1:2*n_r] .* r_span) / R_T

    # Explicit differential variables
    differential_vars = [true for _ in 1:(3*n_r+3)]
    differential_vars[2*n_r+3:3*n_r+2] .= false
    differential_vars[n_r] = false # Algebraic equation aproach for T_N

    return u₀, du₀, differential_vars
end
