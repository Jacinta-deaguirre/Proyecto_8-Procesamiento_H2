using Revise
using HydrogenAdsorption
using Sundials
using Statistics
using Symbolics
using Plots
using CSV
using DataFrames
using DataInterpolations
using DifferentialEquations

# Parameters
# Material properties
# Activated carbon properties
cₛ = 825            # Specific heat capacity of carbon / J/kg K
mₛ = 0.440          # Mass of activated carbon / kg
ρₛ = 1990           # Density of activated carbon / kg/m³
kₛ = 0.9            # Thermal conductivity of activated carbon / W/m K [Not Necessary for this case]
k_eff = 0.21       # Effective thermal conductivity / W/(m·K)
ε_b = 0.88         # Bed porosity     

# Hydrogen properties
cₚ = 14700.0        # Specific heat capacity of hydrogen / J/kg K
M_H2 = 2.0159e-3    # Molar mass of hydrogen / kg/mol
R = 8.314           # Ideal gas constant / J/mol K
k_g = 0.206         # Thermal conductivity of hydrogen / W/m K

# Generate struct for material properties
material_props = MaterialProperties(ρₛ, cₛ, mₛ, kₛ, ε_b, cₚ, M_H2, R, k_g, k_eff)

# Isotherm parameters
# Dubinin-Astakov parameters
P_lim = 77.75e6     # Pa
ψ = 7.3235          # mmol g-1
β = -0.0088         # mol kg-1 K-1
κ = 772.92          # J mol-1
γ = 18.828
m = 2.0            # Exponent in the isotherm equation
DA_params = DAParameters(P_lim, ψ, β, κ, γ, m)

# Tank parameters
L = 255e-3          # Length of the tank / m
R_T = 96e-3         # Interal radius of the tank / m
V = π * R_T^2 * L   # Volume of tank / m³
U = 36              # Heat transfer coefficient / W/(m²·K)
T₀ = 295.0          # Initial temperature of the tank / K
T_air = 295.0       # Ambient temperature / K
dr = 0.00025 / 2    # Radial step size / m
n_r, r_span, A, b = coefficient_matrix(R_T, dr, k_eff, U, T_air)
geometric_params = GeometricParameters(n_r, dr, V, A, b, r_span, R_T)

# Operational parameters
m_in = 9.5e-4 # Mass flow rate of hydrogen / m3/s
operational_params = OperationalParameters(U, T_air, m_in)

# Adsorption system parameters
par = AdsorptionParameters(DA_params, material_props, geometric_params, operational_params)

# Variable flow rate function
"""
Logistic fit for the “9.5×10^-4 m³/s” dataset:

    volumetric_flow_fit(t) = L1 / (1 + exp(k1 * (t - t01)))

"""
function volumetric_flow_fit(t)
    L1 = 0.000974    # plateau flow [m³/s]
    k1 = 0.0451      # 1/s
    t01 = 175.6       # s
    return L1 ./ (1 .+ exp.(k1 .* (t .- t01)))
end

flow_plt = plot(volumetric_flow_fit, 0, 500, label="Mass Flow Rate", xlabel="Time (s)", ylabel="Flow Rate (m³/s)", title="Variable Mass Flow Rate", legend=:topright)
display(flow_plt)

# Generate heat of adsorption function
@variables P T
symbolic_expresion = adsorption_isotherm(DA_params, P, T)
dnₐ_dT = Symbolics.derivative(symbolic_expresion, T)
dnₐ_dP = Symbolics.derivative(symbolic_expresion, P)
dH = -R * T^2 / P * dnₐ_dT / dnₐ_dP
ΔH = build_function(dH, T, P, expression=Val{false})

# Find initial conditions
Tᵢ = ones(n_r) * T₀ # Initial temperature / K
Pᵢ = 0.102564103e6 # Initial pressure / Pa
nₐᵢ = adsorption_isotherm(DA_params, Pᵢ, Tᵢ)
ρᵢ = ideal_gas_equation(Tᵢ, R, M_H2, P=Pᵢ)

# Evaluate dH at initial conditions
#dH_ᵢ = ΔH.(Tᵢ, Pᵢ)
#display(dH_ᵢ)

df_charge_pressure = CSV.read("hermosilla_pressure_exp.csv", DataFrame)
pressure_interpolation = CubicSpline(df_charge_pressure.P, df_charge_pressure.t)
pressure_plt = plot(pressure_interpolation)
display(pressure_plt)

# Find initial conditions with new function
u₀, du₀, differential_vars = dae_setup(DA_params, material_props, geometric_params, operational_params, Pᵢ, T₀)

function adsorption!(out, du, u, p, t)
    # Unpack parameters
    # Structs
    DA_params = p.isotherm
    material_props = p.material
    geometric_params = p.geometric
    operational_params = p.operational

    # Scalars from structs
    # Material properties
    ρₛ = material_props.ρₛ
    ε_b = material_props.ε_b
    M_H2 = material_props.M_H2
    R = material_props.R
    k_eff = material_props.k_eff

    # Geometric parameters
    n_r = geometric_params.n_r
    dr = geometric_params.dr
    V = geometric_params.V
    r_span = geometric_params.r_span
    R_T = geometric_params.R_T

    # Operational parameters
    U = operational_params.U
    m_in = operational_params.m_in

    # Unpack state variables
    T = u[1:n_r]
    nₐ = u[n_r+1:2*n_r]
    ρ_avg = u[2*n_r+1]
    P = u[2*n_r+2]
    ρ = u[2*n_r+3:end]

    P_charge = pressure_interpolation(t) * 1e5 # Pa
    ρ_charge = ideal_gas_equation(last(T), R, M_H2, P=P_charge)
    mass_flow = volumetric_flow_fit(t) .* ρ_charge # kg/s
    
    # Isosteric heat of adsorption
    dH = ΔH.(T, P)

    # Heat equation
    out[1:n_r] .= du[1:n_r] .- heat_equation(material_props, geometric_params, u, du, dH)

    # Neumann BC time derivative for the tank centre 
    out[1] = du[1] - (4 * du[2] - du[3]) / 3

    # Robin BC time derivative to apply method of lines
    out[n_r] = du[n_r] - (4 * du[n_r-1] - du[n_r-2]) / (3 + 2 * U * dr / k_eff)

    # Adsorption isotherm
    out[n_r+1:2*n_r] .= nₐ .- adsorption_isotherm(DA_params, P, T)

    # Macroscopic mass balance
    # mean(n_a .* r_span) / R computes the average adsorption of H2 
    out[2*n_r+1] = du[2*n_r+1] - (mass_flow / (V * ε_b) - ρₛ * (1 - ε_b) * M_H2 / ε_b * mean(du[n_r+1:2*n_r] .* r_span) / R_T)

    # Ideal gas equation
    out[2*n_r+2] = du[2*n_r+2] - ideal_gas_equation(T, du[1:n_r], R, M_H2, R_T, r_span, ρ_avg, du[2*n_r+1])

    # Density of the gas
    out[2*n_r+3:end] = P .- ideal_gas_equation(T, R, M_H2, ρ=ρ)
end

t₀ = 0.0 # Initial time // s
t_f = 500 # Final time // s
tspan = (t₀, t_f) # Time span for the simulation

#println("Testing isosteric heat of adsorption...")
#dH_test = isosteric_heat_of_adsorption(DA_params, Pᵢ, Tᵢ[1])

# Create the DAE problem
println("Simulating...")
prob = DAEProblem(adsorption!, du₀, u₀, tspan, p=par, differential_vars=differential_vars);
prob = remake(prob, p=par)
display(DifferentialEquations.EnsembleThreads())
sol = DifferentialEquations.solve(prob, IDA(linear_solver=:LapackDense), DifferentialEquations.EnsembleThreads())
println("End of simulation...")

# Extract the solution
t = sol.t
T = [sol.u[i][1:n_r] for i in 1:length(sol.u)]
nₐ = [sol.u[i][n_r+1:2*n_r] for i in 1:length(sol.u)]
ρ_avg = [sol.u[i][2*n_r+1] for i in 1:length(sol.u)]
P = [sol.u[i][2*n_r+2] for i in 1:length(sol.u)]
ρ = [sol.u[i][2*n_r+3:end] for i in 1:length(sol.u)]

# Extract the first element (center temperature) of the Temperature vector for each time instance
T_center = [T[i][1] for i in eachindex(T)]

# Experimental data
t_exp = [2.591936954156509, 26.90575469756226, 57.64671979298379, 81.34301761402062,
    97.34337048254767, 121.02814126503364, 149.8078630870115, 179.20592819125463,
    210.49113417825734, 244.2999382480077, 301.6922398329755, 362.8991678183902,
    431.1186520422266, 506.36057282324214]

T_exp = [293.7161290322581, 299.90967741935486, 310.43870967741935, 318.56774193548387,
    323.05806451612904, 330.10322580645163, 336.2193548387097, 340.47741935483873,
    342.18064516129033, 341.1741935483871, 337.69032258064516, 332.89032258064515,
    327.4709677419355, 322.36129032258066]

plt_T_center = plot(t, T_center, xlabel="Time (s)", ylabel="Center Temperature (K)", title="Temperature at Tank Center vs Time")
scatter!(plt_T_center, t_exp, T_exp, label="Experimental Temperature vs Time", legend=false, markersize=4, color=:red)
display(plt_T_center)

### Plotting ###
r_span = range(0, stop=R_T, length=n_r) # Generates radial nodes // m
generate_profiles_plot(t, r_span, T, nₐ, P, ρ, 50, :tab20b) # Generates the profiles plot for time_step = 200s
#generate_profiles_plot(t, r_span, T, nₐ, P, ρ, 100, :tab20b) # Generates the profiles plot for time_step = 100s
#generate_profiles_plot(t, r_span, T, nₐ, P, ρ, 50, :tab20b) # Generates the profiles plot for time_step = 50s

dH_vals = ΔH.(T_center, P)
plot(t, dH_vals, xlabel="Time (s)", ylabel="Isosteric Heat of Adsorption (kJ/mol)", label="Isosteric Heat of Adsorption vs Time", legend=false, size=(1280, 720), margin=Plots.cm)

m_in_vals = volumetric_flow_fit(t) .* ideal_gas_equation([last(T[i]) for i in eachindex(T)], R, M_H2, P=pressure_interpolation(t) * 1e5)
plot(t, m_in_vals, xlabel="Time (s)", ylabel="Mass Flow Rate (kg/s)", label="Mass Flow Rate vs Time", legend=false, size=(1280, 720), margin=Plots.cm)