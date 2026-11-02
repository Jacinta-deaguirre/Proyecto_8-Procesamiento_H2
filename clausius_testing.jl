using Revise
using HydrogenAdsorption
using Sundials
using Statistics
using Plots
using DataFrames
using CSV

# parameters
# Tank parameters
V = 2.4946e-3 # Volume of tank / m³
L = 0.4 # Length of the tank / m

# Inputs needed for the coefficient_matrix function
R_T = sqrt(V / (pi * L)) # Radius of the tank / m
dr = 0.0025 / 2 # Radial step size / m
k_eff = 0.4304 # Effective thermal conductivity / W/(m·K)
U = 36 # Heat transfer coefficient / W/(m²·K)
T₀ = 281.0 # Initial temperature of the tank / K
T_air = 281

n_r, r_span, A, b = coefficient_matrix(R_T, dr, k_eff, U, T_air)

# Isotherm parameters
# Modified Dubinin-Astakov Isotherm parameters
α = 3080.0 # Enthalpic Factor / J/mol
β = 18.9 # Entropic Factor / J/mol K
m = 2.0 # Exponential factor
p₀ = 1470e6 # Saturation pressure / Pa
n₀ = 71.6 # Limit adsoption / mol/kg

MDA_params = MDAParameters(n₀, p₀, α, β, m)

# Dubinin-Astakov parameters
P_lim = 77.75e6     # Pa
ψ = 7.3235          # mmol g-1
β = -0.0088         # mol kg-1 K-1
κ = 772.92          # J mol-1
γ = 18.828
m = 2.0            # Exponent in the isotherm equation
DA_params = DAParameters(P_lim, ψ, β, κ, γ, m)

# Material properties
# Activated carbon properties
ρₛ = 517.6 # Density of activated carbon / kg/m³
cₛ = 825 # Specific heat capacity of carbon / J/kg K
mₛ = 0.671 # Mass of activated carbon / kg
kₛ = 0.646 # Thermal conductivity of activated carbon / W/m K
ε_b = 0.49 # Bed porosity

# Hydrogen properties
cₚ = 14700.0 # Specific heat capacity of hydrogen / J/kg K
M_H2 = 2.0159e-3 # Molar mass of hydrogen / kg/mol
R = 8.314 # Ideal gas constant / J/mol K
k_g = 0.206 # Thermal conductivity of hydrogen / W/m K

k_eff = (kₛ * (1 - ε_b) + k_g * ε_b) # Effective thermal conductivity / W/(m·K)
material_props = MaterialProperties(ρₛ, cₛ, mₛ, kₛ, ε_b, cₚ, M_H2, R, k_g, k_eff)

# Geometric parameters
geometric_params = GeometricParameters(n_r, dr, V, A, b, r_span, R_T)

# Operational Parameters
U = 36.0 # Heat transfer coefficient / W/(m²·K)
T_air = 281.0 # Ambient temperature / K
m_in = 2.023e-5 # Mass flow rate of hydrogen / kg / s
operational_params = OperationalParameters(U, T_air, m_in)

# Adsorption system parameters
par = AdsorptionParameters(MDA_params, material_props, geometric_params, operational_params)

# Find initial conditions
Tᵢ = ones(n_r) * T₀ # Initial temperature / K
Pᵢ = 0.033e6 # Initial pressure / Pa
nₐᵢ = adsorption_isotherm(MDA_params, Pᵢ, Tᵢ)
ρᵢ = ideal_gas_equation(Tᵢ, R, M_H2, P=Pᵢ)

# Find intiial conditions with new function
u₀, du₀, differential_vars = dae_setup(MDA_params, material_props, geometric_params, operational_params, Pᵢ, T₀)

function adsorption!(out, du, u, p, t)
    # Unpack parameters
    # Structs
    MDA_params = p.isotherm
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

    # Isosteric heat of adsorption
    dH = isosteric_heat_of_adsorption(P, du[2*n_r+2], T, du[1:n_r], R)

    # Heat equation
    out[1:n_r] .= du[1:n_r] .- heat_equation(material_props, geometric_params, u, du, dH)

    # Neumann BC time derivative for the tank centre 
    out[1] = du[1] - (4 * du[2] - du[3]) / 3

    # Robin BC time derivative to apply method of lines
    out[n_r] = du[n_r] - (4 * du[n_r-1] - du[n_r-2]) / (3 + 2 * U * dr / k_eff)

    # Adsorption isotherm
    out[n_r+1:2*n_r] .= nₐ .- adsorption_isotherm(MDA_params, P, T)

    # Macroscopic mass balance
    # mean(n_a .* r_span) / R computes the average adsorption of H2 
    out[2*n_r+1] = du[2*n_r+1] - (m_in / (V * ε_b) - ρₛ * (1 - ε_b) * M_H2 / ε_b * mean(du[n_r+1:2*n_r] .* r_span) / R_T)

    # Ideal gas equation
    out[2*n_r+2] = du[2*n_r+2] - ideal_gas_equation(T, du[1:n_r], R, M_H2, R_T, r_span, ρ_avg, du[2*n_r+1])

    # Density of the gas
    out[2*n_r+3:end] = P .- ideal_gas_equation(T, R, M_H2, ρ=ρ)
end

t₀ = 0.0 # Initial time // s
t_f = 1042 # Final time // s
tspan = (t₀, t_f) # Time span for the simulation

# Create the DAE problem
prob = DAEProblem(adsorption!, du₀, u₀, tspan, p=par, differential_vars=differential_vars);
prob = remake(prob, p=par);
sol = solve(prob, IDA())

# Extract the solution
t = sol.t
T = [sol.u[i][1:n_r] for i in 1:length(sol.u)]
nₐ = [sol.u[i][n_r+1:2*n_r] for i in 1:length(sol.u)]
ρ_avg = [sol.u[i][2*n_r+1] for i in 1:length(sol.u)]
P = [sol.u[i][2*n_r+2] for i in 1:length(sol.u)]
ρ = [sol.u[i][2*n_r+3:end] for i in 1:length(sol.u)]

### Plotting ###
r_span = range(0, stop=R_T, length=n_r) # Generates radial nodes // m
#generate_profiles_plot(t, r_span, T, nₐ, P, ρ, 200, :tab20b) # Generates the profiles plot for time_step = 200s
#generate_profiles_plot(t, r_span, T, nₐ, P, ρ, 100, :tab20b) # Generates the profiles plot for time_step = 100s
generate_profiles_plot(t, r_span, T, nₐ, P, ρ, 50, :tab20b) # Generates the profiles plot for time_step = 50s

### Average properties ###
T_avg = [mean(T[i]) for i in eachindex(T)]
n_avg = [mean(nₐ[i]) for i in eachindex(nₐ)]
ρ_avg_nodes = [mean(ρ[i]) for i in eachindex(ρ)]

# Experimental data from Xiao et al. 2013
# Temperature data
t_exp = [18.633540372670836, 204.96894409937886, 409.9378881987577, 605.5900621118014, 805.9006211180124, 1006.2111801242236]
T_exp = [281.0010449320794, 285.6530825496343, 288.43260188087777, 290.15882967607104, 291.1243469174504, 291.6802507836991]
# Pressure data
t_exp2 = [0, 189.8370086, 402.6845638, 598.274209, 799.6164909, 1000.958773]
p_exp = [0.102564103, 1.487179487, 3.115384615, 4.858974359, 6.666666667, 8.58974359]
# Adsorption data
t_expna = [189.364710393251476, 398.85705983848175, 594.8740683936833, 802.068106113255, 1000.2396414963216]
na_exp = [0.002467612451772145, 0.004852345640921178, 0.007182067626830262, 0.00942904215293922, 0.011543818447602383] / M_H2 / mₛ
# Total mass data
t_expmt = [190.04768865776782, 397.97862397852805, 596.90503007500, 804.8179922835443, 1003.7803446044716]
mt_exp = [0.004142958134630587, 0.008197542476455221, 0.012164016631119842, 0.016174512928658725, 0.020229163171894844]
m_gas = mt_exp .- [0.002467612451772145, 0.004852345640921178, 0.007182067626830262, 0.00942904215293922, 0.011543818447602383]
ρ_exp = m_gas ./ (V * ε_b)

plt = plot(layout=(2, 2)) # Displat in a 2x2 grid

# Temperature
# Line plot for the average temperature and scatter plot for the experimental data
plot!(plt[1], t, T_avg, xlabel="Time (s)", ylabel="Average Temperature (K)", label="Average Temperature vs Time", legend=false)
scatter!(plt[1], t_exp, T_exp, label="Experimental Temperature vs Time", legend=false, markersize=4, color=:red)

# Adsorption
# Line plot for the average adsorption and scatter plot for the experimental data
plot!(plt[2], t, n_avg, xlabel="Time (s)", ylabel="q_a (mol/kg_ads)", label="Average Adsorption vs Time", legend=false)
scatter!(plt[2], t_expna, na_exp, label="Experimental Adsorption vs Time", legend=false, markersize=4, color=:pink)

# Pressure
# Line plot for the pressure and scatter plot for the experimental data
plot!(plt[3], t, P .* 1e-6, xlabel="Time (s)", ylabel="Pressure (MPa)", label="Pressure vs Time", legend=false)
scatter!(plt[3], t_exp2, p_exp, label="Experimental Pressure vs Time", legend=false, markersize=4, color=:orange)

# Gas density
# Line plot for the average gas density and scatter plot for the experimental data
plot!(plt[4], t, ρ_avg_nodes, xlabel="Time (s)", ylabel="Average Gas Density (kg/m³)", label="Average Gas Density vs Time", legend=false)
scatter!(plt[4], t_expmt, ρ_exp, label="Experimental Gas Density vs Time", legend=false, markersize=4, color=:green)

plot!(size=(1280, 720))
plot!(margin=Plots.cm)
display(plt)

∇T = [sol.du[i][1:n_r] for i in 1:length(sol.du)]
dT_avg = [mean(∇T[i]) for i in eachindex(∇T)]
∇P = [sol.du[i][2*n_r+2] for i in 1:length(sol.du)]
∇nₐ= nₐ = [sol.du[i][n_r+1:2*n_r] for i in 1:length(sol.du)]
dH_vals = isosteric_heat_of_adsorption(P, ∇P, T_avg, dT_avg, R) # Evaluate the isosteric heat of adsorption with Clausius
plot(t, dH_vals, xlabel="Time (s)", ylabel="Isosteric Heat of Adsorption (J/mol)", label="Isosteric Heat of Adsorption vs Time", legend=false, size=(1280, 720), margin=Plots.cm)
df = DataFrame(t=t, T_avg=T_avg, n_avg=n_avg, ρ_avg_nodes=ρ_avg_nodes, P=P .* 1e-6, dH=dH_vals)
CSV.write("xiao_clausius_simulation.csv", df) # Save the results to a CSV file


