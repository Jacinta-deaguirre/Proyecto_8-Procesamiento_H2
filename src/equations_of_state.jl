"""
Non differentiated ideal gas equation.

# Inputs:
- `T`: Temperature in the tank / Vector of temperatures at each radial node / K
- `R`: Ideal gas constant / J/(mol·K)
- `M`: Molar mass of the gas / kg/mol
- `ρ`: Density of the gas / kg/m³ (optional, if provided, pressure will be calculated)
- `P`: Pressure of the gas / Pa (optional, if provided, density will be calculated)

# Outputs:
- `ρ`: Density of the gas / kg/m³ if `P` is provided
- `P`: Pressure of the gas / Pa if `ρ` is provided
"""
function ideal_gas_equation(T, R::Float64, M::Float64; ρ=nothing, P=nothing)
    if isnothing(ρ) && isnothing(P)
        error("Either density (ρ) or pressure (P) must be provided.")
    end

    if !isnothing(ρ) && !isnothing(P)
        error("Only one of density (ρ) or pressure (P) should be provided.")
    end
    
    if !isnothing(ρ)
        # Calculate pressure from density
        P = ρ .* R .* T ./ M
        return P
    else
        # Calculate density from pressure
        ρ = P .* M ./ (R .* T)
        return ρ
    end
end

"""
Differentiated ideal gas equation.

# Inputs:
- `T`: Temperature in the tank / Vector of temperatures at each radial node / K
- `dTdt`: Rate of change of temperature / Vector of rates of change of temperature at each radial node / K/s
- `R`: Ideal gas constant / J/(mol·K)
- `M`: Molar mass of the gas / kg/mol
- `R_T`: Radius of the tank / m
- `r_span`: Radial span of the tank / Vector of radial positions at each node / m
- `ρ`: Density of the gas / kg/m³
- `dρdt`: Rate of change of density / Vector of rates of change of density at each radial node / kg/m³/s

# Outputs:
- `dPdt`: Rate of change of pressure / Pa/s

Given that the temperature input is a vector, 
the function computes the average temperature and the average rate of change of temperature across the radial span of the tank. 
It then uses these averages to compute the rate of change of pressure using the ideal gas law.
"""

function cyl_avg(f, r_span, R_T)
    dr = r_span[2] - r_span[1]
    n_r=length(r_span)
    f_avg = (2/R_T^2) * sum( (f[2:n_r] .* r_span[2:end] + f[1:n_r-1] .* r_span[1:end-1]) / 2 * dr)
    return f_avg
end

function ideal_gas_equation(T::Vector{Float64}, dTdt::Vector{Float64}, R::Float64, M::Float64, R_T::Float64, r_span ,ρ::Float64, dρdt::Float64)
    # T_avg = (mean(T .* r_span) / R_T)
    T_avg = cyl_avg(T, r_span, R_T)
    dTdt_avg = cyl_avg(dTdt, r_span, R_T)

    dPdt = R / M * (dρdt * T_avg + ρ * dTdt_avg)
    return dPdt
end