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
`mₛ` must be the adsorbent mass implied by the bed, ρₛ (1 - ε_b) V, so that the energy and
mass balances use the same adsorbent mass.
"""
struct MaterialProperties
    # Activated carbon properties
    ρₛ::Float64        # Density of the solid phase / kg/m³
    cₛ::Float64        # Specific heat capacity of the solid phase / J/(kg·K)
    mₛ::Float64        # Mass of the solid phase = ρₛ (1 - ε_b) V / kg
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
    dr::Float64             # Radial step size (actual grid spacing r[2]-r[1]) / m
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
    U::Float64      # External heat transfer coefficient (wall-bath) / W/(m²·K)
    T_air::Float64  # Temperature of the external fluid (water bath) / K
    m_in            # Mass flow rate of hydrogen into the tank (negative for discharge) / kg/s
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

    # Universal gas constant / J/(mol·K)
    R = 8.314

    # Modified Dubinin-Astakov isotherm equation
    # P is clamped to a small positive value so that a rejected trial step of the DAE solver
    # (which may probe P < 0) does not throw a DomainError; it has no effect on valid states.
    nₐ = n₀ .* exp.(-(R .* T ./ (α .+ β .* T)) .^ m .* (log.(p₀ ./ max.(P, 1.0))) .^ m)

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

    # Universal gas constant / J/(mol·K)
    R = 8.314

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
    dH = α .* ((log.(n₀ ./ min.(adsorption_isotherm(params, P, T), 0.999 * n₀))) .^ (1 / m))

    return dH
end

"""
DAE residual of the 1-D radial model (method of lines), `out = F(du, u, p, t)`, with
`p::AdsorptionParameters`. This is the single definition of the model used by every script.

State vector `u` (length 3 n_r + 3):
- `u[1:n_r]`            bed temperature at each radial node, T / K (differential; node 1
                        is the symmetry condition, node n_r the bed-wall interface, both
                        algebraic constraints)
- `u[n_r+1:2n_r]`       absolute adsorption at each node, n_a / mol/kg (algebraic, Eq. 1)
- `u[2n_r+1]`           volume-averaged gas density ρ̄ / kg/m³ (differential, Eq. 3)
- `u[2n_r+2]`           pressure P / Pa (algebraic: P = ρ̄ R_H2 T_H, T_H = 1/<1/T>)
- `u[2n_r+3:3n_r+2]`    local gas density at each node (algebraic, local ideal gas)
- `u[3n_r+3]`           lumped wall temperature T_w / K (differential, Eq. 8)

Equations:
- energy balance (enthalpy form with compression work), see `heat_equation`
- symmetry at r = 0: d/dt(-3T_1 + 4T_2 - T_3) = 0
- bed-wall interface: k_eff ∂T/∂r|_{r_i} = k_w (T_w - T(r_i)) / (e/2), second-order
  one-sided difference, imposed algebraically
- volume-averaged mass balance, Eq. (3), with the trapezoidal cylindrical average of ∂n_a/∂t
- ideal-gas closure with the harmonic volume-mean temperature (exact for uniform P)
- wall energy balance, Eq. (8)
"""
function MDA_adsorption!(out, du, u, p::AdsorptionParameters, t)
    iso = p.isotherm
    mp = p.material
    gp = p.geometric
    op = p.operational

    ρₛ, ε_b, M_H2, R, k_eff = mp.ρₛ, mp.ε_b, mp.M_H2, mp.R, mp.k_eff
    n_r, dr, V, r_span, R_T, e = gp.n_r, gp.dr, gp.V, gp.r_span, gp.R_T, gp.e
    Ao, Ai, c_wall, m_wall, k_wall = gp.Ao, gp.Ai, gp.c_wall, gp.m_wall, gp.k_wall
    U, T_air, m_in = op.U, op.T_air, op.m_in

    # Unpack state variables
    T = u[1:n_r]
    nₐ = u[n_r+1:2*n_r]
    ρ_avg = u[2*n_r+1]
    P = u[2*n_r+2]
    ρ = u[2*n_r+3:3*n_r+2]
    T_wall = u[3*n_r+3]

    # Isosteric heat of adsorption
    dH = isosteric_heat_of_adsorption(iso, P, T)

    # Energy balance
    out[1:n_r] .= du[1:n_r] .- heat_equation(mp, gp, op, u, du, dH)

    # Symmetry condition at the axis (time derivative of the one-sided difference)
    out[1] = du[1] - (4 * du[2] - du[3]) / 3

    # Bed-wall interface: flux continuity with the lumped wall (algebraic)
    out[n_r] = (3 + 4 * k_wall * dr / k_eff / e) * T[n_r] - (4 * T[n_r-1] - T[n_r-2] + 4 * k_wall * dr / k_eff / e * T_wall)

    # Adsorption equilibrium at each node
    out[n_r+1:2*n_r] .= nₐ .- adsorption_isotherm(iso, P, T)

    # Volume-averaged mass balance, Eq. (3)
    dna_avg = cyl_avg(du[n_r+1:2*n_r], r_span, R_T)
    out[2*n_r+1] = du[2*n_r+1] - (m_in / (V * ε_b) - ρₛ * (1 - ε_b) * M_H2 / ε_b * dna_avg)

    # Ideal-gas closure with the harmonic volume-mean temperature
    T_H = 1 / cyl_avg(1 ./ T, r_span, R_T)
    out[2*n_r+2] = P - ρ_avg * R / M_H2 * T_H

    # Local gas density
    out[2*n_r+3:3*n_r+2] .= P .- ρ .* R .* T ./ M_H2

    # Tank wall energy balance, Eq. (8)
    Q_bed = -k_wall * Ai * (T_wall - T[n_r]) / (e / 2)          # Heat flow from the bed to the wall
    Q_wall = (e / 2 / k_wall + 1 / U)^(-1) * Ao * (T_wall - T_air) # Heat flow from the wall to the bath
    out[3*n_r+3] = du[3*n_r+3] - (Q_bed - Q_wall) / (c_wall * m_wall)
end

"""
Initial state, initial derivative guess and differential-variable flags for the DAE.

- `Pᵢ`, `T_0`: initial pressure and (uniform) temperature. For a stage restart, pass the
  final state of the previous stage through the keyword arguments `T_init` (vector),
  `ρ_avg_init`, `T_wall_init` and `du_init` (final derivative vector); the integrated ρ̄ is then
  carried over exactly, which keeps the global mass balance closed across stages.
"""
function dae_setup(isotermParams::IsothermParameters,
    materialProps::MaterialProperties,
    geometricParams::GeometricParameters,
    operationalParams::OperationalParameters,
    Pᵢ::Float64,
    T_0::Float64,
    R=8.314;
    T_init=nothing, ρ_avg_init=nothing, T_wall_init=nothing, du_init=nothing)

    # Unpack parameters
    M_H2, k_eff, ε_b, ρₛ = materialProps.M_H2, materialProps.k_eff, materialProps.ε_b, materialProps.ρₛ
    U, T_air, m_in = operationalParams.U, operationalParams.T_air, operationalParams.m_in
    n_r, dr, V, R_T, r_span = geometricParams.n_r, geometricParams.dr, geometricParams.V, geometricParams.R_T, geometricParams.r_span
    e = geometricParams.e
    k_wall = geometricParams.k_wall

    # Initial temperature vector
    Tᵢ = T_init === nothing ? ones(n_r) * T_0 : copy(T_init)

    # Initial adsorption vector
    nₐᵢ = adsorption_isotherm(isotermParams, Pᵢ, Tᵢ)

    # Initial local density vector
    ρᵢ = ideal_gas_equation(Tᵢ, R, M_H2, P=Pᵢ)

    # Initial volume-averaged density, consistent with P = ρ̄ R_H2 T_H
    ρ_avgᵢ = ρ_avg_init === nothing ? Pᵢ * M_H2 / R * cyl_avg(1 ./ Tᵢ, r_span, R_T) : ρ_avg_init

    # Initial wall temperature
    T_wallᵢ = T_wall_init === nothing ? Tᵢ[end] : T_wall_init

    # Initial state vector
    u₀ = vcat(Tᵢ, nₐᵢ, ρ_avgᵢ, Pᵢ, ρᵢ, T_wallᵢ)

    if T_init === nothing
        # Interface temperature consistent with the algebraic bed-wall condition
        u₀[n_r] = (4 * u₀[n_r-1] - u₀[n_r-2] + 4 * k_wall * dr / k_eff / e * T_wallᵢ) / (3 + 4 * k_wall * dr / k_eff / e)
        u₀[2*n_r] = adsorption_isotherm(isotermParams, Pᵢ, u₀[n_r])
        u₀[3*n_r+2] = ideal_gas_equation(u₀[n_r], R, M_H2, P=Pᵢ)
    end

    # Initial time derivative guess (IDA computes consistent values). At a stage restart the
    # final derivative of the previous stage must be supplied (du_init) so that the
    # consistent-initialisation Newton iteration of IDA starts close to the solution.
    du₀ = du_init === nothing ? ones(length(u₀)) * 1e-5 : copy(du_init)
    du_init === nothing && (du₀[n_r] = (4 * du₀[n_r-1] - du₀[n_r-2]) / (3 + 2 * U * dr / k_eff))
    du₀[2*n_r+1] = m_in / V / ε_b - ρₛ * (1 - ε_b) * M_H2 / ε_b * cyl_avg(du₀[n_r+1:2*n_r], r_span, R_T)

    # Differential / algebraic flags
    differential_vars = [true for _ in 1:(3*n_r+3)]
    differential_vars[n_r] = false                      # bed-wall interface
    differential_vars[n_r+1:2*n_r] .= false             # adsorption equilibrium
    differential_vars[2*n_r+2] = false                  # ideal-gas closure
    differential_vars[2*n_r+3:3*n_r+2] .= false         # local densities

    return u₀, du₀, differential_vars
end
