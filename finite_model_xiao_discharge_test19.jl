using Revise
using HydrogenAdsorption
using Sundials
using Statistics
using Plots
using CSV
using DataFrames
VERBOSE = false # Set to true to print more information during the simulation

# parameters
# Tank parameters
V = 2.4946e-3 # Volume of tank / m³
L = 0.4 # Length of the tank / m
e = 0.0508 - 0.0469 # Wall thickness / m
c_wall = 468 # Specific heat capacity of the tank wall / J/kg K
m_wall = 3.714 # Mass of the tank wall / kg
k_wall = 13 # Thermal conductivity of the tank wall / W/m K

# Inputs needed for the coefficient_matrix function
R_T = sqrt(V / (pi * L)) # Radius of the tank / m
dr = (1 / 2)^5 * 0.0025 # Radial step size / m
U = 36 # Heat transfer coefficient / W/(m²·K)

# Tank surface areas
Ai = 2 * pi * R_T * L # Inner surface area of the tank / m²
Ao = 2 * pi * (R_T + e) * L # Outer surface area of the tank / m²

T₀ = 281.0 # Initial temperature of the tank / K
#T_air = 281.0 # Ambient temperature / K
T_air = 302.5 # Ambient temperature / K

#T_H2 = 281.6 # Temperature of the incoming hydrogen gas / K
T_H2 = 277.7 # Data from finite element model / K


# Isotherm parameters
# Modified Dubinin-Astakov Isotherm parameters
α = 3080.0 # Enthalpic Factor / J/mol
β = 18.9 # Entropic Factor / J/mol K
m = 2.0 # Exponential factor
p₀ = 1470e6 # Saturation pressure / Pa
n₀ = 71.6 # Limit adsoption / mol/kg

MDA_params = MDAParameters(n₀, p₀, α, β, m)

# Material properties
# Activated carbon properties
ρₛ = 517.6 # Density of activated carbon / kg/m³
cₛ = 825 # Specific heat capacity of carbon / J/kg K
mₛ = 0.671 # Mass of activated carbon / kg
kₛ = 0.646 # Thermal conductivity of activated carbon / W/m K
ε_b = 0.49 # Bed porosity

# Hydrogen properties
cₚ = 14700.0 # Specific heat capacity of hydrogen / J/kg K
cᵥ = 10134.0 # Specific heat capacity at constant volume / J/kg K
M_H2 = 2.0159e-3 # Molar mass of hydrogen / kg/mol
R = 8.314 # Ideal gas constant / J/mol K
k_g = 0.206 # Thermal conductivity of hydrogen / W/m K

k_eff = (kₛ * (1 - ε_b) + k_g * ε_b) # Effective thermal conductivity / W/(m·K)
n_r, r_span, A, b = coefficient_matrix(R_T, dr, k_eff, U, T_air)
material_props = MaterialProperties(ρₛ, cₛ, mₛ, kₛ, ε_b, cₚ, cᵥ, M_H2, R, k_g, k_eff)

# Geometric parameters
geometric_params = GeometricParameters(n_r, dr, V, L, A, b, r_span, R_T, e, Ao, Ai, c_wall, m_wall, k_wall)

# Operational Parameters
U = 1 * 36.0 # Heat transfer coefficient / W/(m²·K)
m_in = -2.186e-5 # Mass flow rate of hydrogen / kg / s
operational_params = OperationalParameters(U, T_air, m_in, T_H2)

# Adsorption system parameters
par = AdsorptionParameters(MDA_params, material_props, geometric_params, operational_params)

charge_temperature = CSV.read("outputs/xiao19_dormant_temperature_profiles.csv", DataFrame)
charge_final_temp = charge_temperature[!, Symbol("t_3822.0s")]
charge_timeseries_data = CSV.read("outputs/xiao19_dormant_timeseries_data.csv", DataFrame)
charge_final_pressure = charge_timeseries_data[end, :Pressure_Pa]
charge_final_wall_temp = charge_timeseries_data[end, :Wall_Temperature_K]

display("Charge final pressure: $charge_final_pressure Pa")
display("Charge final wall temperature: $charge_final_wall_temp K")

# Find initial conditions
Tᵢ = charge_final_temp # Initial temperature / 
Pᵢ = charge_final_pressure # Initial pressure / Pa
nₐᵢ = adsorption_isotherm(MDA_params, Pᵢ, Tᵢ)
ρᵢ = ideal_gas_equation(Tᵢ, R, M_H2, P=Pᵢ)

# Find intiial conditions with new function
u₀, du₀, differential_vars = dae_setup(MDA_params, material_props, geometric_params, operational_params, Pᵢ, T₀)

u₀[1:n_r] = Tᵢ
u₀[n_r+1:2*n_r] = nₐᵢ
u₀[2*n_r+1] = mean(ρᵢ)
u₀[2*n_r+2] = Pᵢ
u₀[2*n_r+3:3*n_r+2] = ρᵢ
u₀[3*n_r+3] = charge_final_wall_temp # Tank wall temperature

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
    e = geometric_params.e
    Ao = geometric_params.Ao
    Ai = geometric_params.Ai
    c_wall = geometric_params.c_wall
    m_wall = geometric_params.m_wall
    k_wall = geometric_params.k_wall

    # Operational parameters
    U = operational_params.U
    m_in = operational_params.m_in

    # Unpack state variables
    T = u[1:n_r]
    nₐ = u[n_r+1:2*n_r]
    ρ_avg = u[2*n_r+1]
    P = u[2*n_r+2]
    ρ = u[2*n_r+3:3*n_r+2]
    T_wall = u[3*n_r+3] # Tank wall temperature

    # Isosteric heat of adsorption
    dH = isosteric_heat_of_adsorption(MDA_params, P, T)
    #dH = isosteric_heat_of_adsorption(P, du[2*n_r+2], T, du[1:n_r], R)

    # Heat equation
    out[1:n_r] .= du[1:n_r] .- heat_equation(material_props, geometric_params, operational_params, u, du, dH)

    # Neumann BC time derivative for the tank centre 
    out[1] = du[1] - (4 * du[2] - du[3]) / 3

    # Robin BC time derivative to apply method of lines
    #out[n_r] = du[n_r] - (4 * du[n_r-1] - du[n_r-2]) / (3 + 4 * k_wall * dr / k_eff / e)
    out[n_r] = (3 + 4 * k_wall * dr / k_eff / e) * T[n_r] - (4 * T[n_r-1] - T[n_r-2] + 4 * k_wall * dr / k_eff / e * T_wall) # DAE Approach

    # Adsorption isotherm
    out[n_r+1:2*n_r] .= nₐ .- adsorption_isotherm(MDA_params, P, T)

    # Macroscopic mass balance
    # mean(n_a .* r_span) / R computes the average adsorption of H2 
    dna_avg = sum((du[n_r+2:2*n_r] .* r_span[2:end] + du[n_r+1:2*n_r-1] .* r_span[1:end-1]) / 2 * (2 / R_T^2) * dr)
    # dna_avg_simple = mean(du[n_r+1:2*n_r])
    # println("dna_avg: $dna_avg, dna_avg_simple: $dna_avg_simple")
    out[2*n_r+1] = du[2*n_r+1] - (m_in / (V * ε_b) - ρₛ * (1 - ε_b) * M_H2 / ε_b * dna_avg)
    #out[2*n_r+1] = du[2*n_r+1] - (m_in / (V * ε_b) - ρₛ  * M_H2 / ε_b * dna_avg)

    # Ideal gas equation
    out[2*n_r+2] = du[2*n_r+2] - ideal_gas_equation(T, du[1:n_r], R, M_H2, R_T, r_span, ρ_avg, du[2*n_r+1])

    # Density of the gas
    out[2*n_r+3:3*n_r+2] = P .- ideal_gas_equation(T, R, M_H2, ρ=ρ)


    # Tank wall energy balance
    Q_bed = -k_wall * Ai * (T_wall - T[n_r]) / (e / 2) # Heat flow from the bed to the wall
    Q_wall = (e / 2 / k_wall + 1 / U)^(-1) * Ao * (T_wall - T_air) # Heat flow from the wall to the ambient
    out[3*n_r+3] = du[3*n_r+3] - (Q_bed - Q_wall) / (c_wall * m_wall)
end

t₀ = 3822 # Initial time // s
t_f = 4694 # Final time // s
tspan = (t₀, t_f) # Time span for the simulation

# Create the DAE problem
prob = DAEProblem(adsorption!, du₀, u₀, tspan, p=par, differential_vars=differential_vars);
prob = remake(prob, p=par);
# @btime sol = solve(prob, IDA(linear_solver=:LapackDense), progress=true) # 838.288 s for dr = 0.000025
sol = solve(prob, IDA(linear_solver=:LapackDense), progress=true)
# Extract the solution
t = sol.t
T = [sol.u[i][1:n_r] for i in 1:length(sol.u)]
T_wall = [sol.u[i][3*n_r+3] for i in 1:length(sol.u)]
nₐ = [sol.u[i][n_r+1:2*n_r] for i in 1:length(sol.u)]
ρ_avg = [sol.u[i][2*n_r+1] for i in 1:length(sol.u)]
P = [sol.u[i][2*n_r+2] for i in 1:length(sol.u)]
ρ = [sol.u[i][2*n_r+3:3*n_r+2] for i in 1:length(sol.u)]

### Plotting ###
plot_palette = palette(:viridis, 13, rev=true)
r_span = range(0, stop=R_T, length=n_r) # Generates radial nodes // m
#generate_profiles_plot(t, r_span, T, nₐ, P, ρ, 200, :tab20b) # Generates the profiles plot for time_step = 200s
#generate_profiles_plot(t, r_span, T, nₐ, P, ρ, 100, :tab20b) # Generates the profiles plot for time_step = 100s
generate_profiles_plot(t, r_span, T, nₐ, P, ρ, 50, plot_palette) # Generates the profiles plot for time_step = 50s

middle_index = div(n_r, 2) # Index of the middle node
println("Middle index: $middle_index, final index: $(n_r)")


T_center = [T[i][1] for i in eachindex(T)] # Temperature at the tank center
T_middle = [T[i][middle_index] for i in eachindex(T)] # Temperature at the middle of the tank

# Reconstruct final free and adsorber hydrogen mass
#  in the tank 
T_avg = [mean(T[i]) for i in eachindex(T)]
n_avg = [mean(nₐ[i]) for i in eachindex(nₐ)]
ρ_avg_nodes = [mean(ρ[i]) for i in eachindex(ρ)]

m_H2_gas = ρ_avg .* V * ε_b # Mass of hydrogen in the gas phase / kg
m_H2_ads = n_avg .* mₛ * M_H2 # Mass of hydrogen in the adsorbed phase / kg
m_H2_total = m_H2_gas .+ m_H2_ads # Total mass of hydrogen in the tank / kg 
mass_df = DataFrame(Time_s=t[2:end], Gas_Mass_kg=m_H2_gas[2:end], Adsorbed_Mass_kg=m_H2_ads[2:end], Total_Mass_kg=m_H2_total[2:end])
CSV.write("outputs/xiao19_discharge_mass_profiles.csv", mass_df)

if VERBOSE
    println("Final mass of hydrogen in the gas phase: $(m_H2_gas[end] * 1000) g")
    println("Final mass of hydrogen in the adsorbed phase: $(m_H2_ads[end] * 1000) g")
    println("Final total mass of hydrogen in the tank: $(m_H2_total[end] * 1000) g")
    println("Final average temperature in the tank: $(T_avg[end]) K")
    println("Final pressure in the tank: $(P[end] / 1e6) MPa")
    println("Final average gas density in the tank: $(ρ_avg_nodes[end]) kg/m³")
    println("Final average adsorption in the tank: $(n_avg[end]) mol/kg_ads")
end

# Compare with m_in * t_f
println("Total mass of hydrogen injected: $(m_in * t_f * 1000) g")

# Experimental data
# Temperature profiles
df_avg_temp = CSV.read("inputs/xiao19_avg_temperature_discharge.csv", DataFrame)


# --- 1. Global Plot Styling ---
# Apply a theme with larger fonts, similar to the Python 'seaborn-talk' style
theme(:default,
    fontfamily="sans-serif",
    titlefontsize=22,
    guidefontsize=20,     # X/Y axis labels
    tickfontsize=16,
    legendfontsize=10,
    margin=1.5Plots.cm,   # Add generous margins
    grid=true,
    gridstyle=:dash,
    gridalpha=0.5,
    framestyle=:box,      # Adds a full box, looks professional
    palette=:seaborn_pastel   # Use a modern, colorblind-friendly palette
)

# Define common line/marker styles for clarity
# (Simulated = thick solid, Experimental = thinner dash + markers)
sim_lw = 2.5
exp_lw = 1.5
exp_marker = (:circle, 6, 0.7) # (shape, size, alpha)

# --- 2. Create Individual Plots ---
# We create three separate plot objects first

# Plot A: Temperature
p_temp = plot(
    title="Temperature Profiles",
    ylabel="Temperature / K",
    xlabel="Time / s",
    legend=:topright, # Move legend inside
    legend_columns=1,     # Arrange in 3 columns (2 rows)
    size=(1200, 900)
)

# Avg temperature
plot!(p_temp, df_avg_temp.Time, df_avg_temp.Temperature,
    label="Exp. Average", color=2, marker=exp_marker, linestyle=:dash, lw=exp_lw)
plot!(p_temp, t, T_avg,
    label="Sim. Average", color=2, lw=sim_lw)

# Plot B: Pressure
p_pressure = plot(
    title="Pressure Profile",
    ylabel="Pressure / MPa",
    xlabel="Time / s",
    legend=:best, # Move legend inside
    legend_columns=1,     # Arrange in 2 columns (1 row)
    size=(1200, 900)
)


plot!(p_pressure, t, P ./ 1e6,
    label="Simulated", color=1, lw=sim_lw)


# Plot C: Mass Profiles
p_mass = plot(
    title="Mass Profiles",
    ylabel="Mass / g",
    xlabel="Time / s",
    legend=:best, # Move legend inside
    legend_columns=1,     # Arrange in 3 columns (2 rows)
    size=(1200, 900)
)

# Gas mass

plot!(p_mass, t, m_H2_gas .* 1000,
    label="Sim. Gas", color=1, lw=sim_lw)

# Adsorbed mass

plot!(p_mass, t, m_H2_ads .* 1000,
    label="Sim. Adsorbed", color=2, lw=sim_lw)

# Total mass

plot!(p_mass, t, m_H2_total .* 1000,
    label="Sim. Total", color=3, lw=sim_lw)


# --- 3. Combine Plots into a Single Figure ---
# Use a 2x2 layout, placing Temp (a) and Pressure (b) side-by-side,
# and Mass (c) spanning the full width below them.
l = @layout [a b; c{0.5h}] # Bottom plot (c) takes 50% of the total height
final_plot = plot(p_temp, p_pressure, p_mass,
    layout=l,
    size=(1800, 1400) # A large, wide figure suitable for presentations
)

# Display the final combined plot
display(final_plot)


# --- 4. Save the Figure ---
# Save as SVG (vector format) for maximum quality, just like the Python script
output_filename = "outputs/xiao19_discharge_simulation_profiles.svg"
savefig(final_plot, output_filename)
println("\nPlot saved successfully as '$output_filename'")

# Save individuals
savefig(p_temp, "outputs/xiao19_discharge_julia_temperature_profile.svg")
savefig(p_pressure, "outputs/xiao19_discharge_julia_pressure_profile.svg")
savefig(p_mass, "outputs/xiao19_discharge_julia_mass_profiles.svg")


# --- Export Results to CSV Files ---
println("Exportando resultados a archivos CSV (Formato 'Ancho')...")

# --- 1. Preparar Datos de Perfil (n_r x n_t) ---

# Convertir el Vector de Vectores de la solución a una Matriz
# El resultado es una matriz de [n_r filas x n_t columnas]
T_matrix = reduce(hcat, T)
n_a_matrix = reduce(hcat, nₐ)

# Crear encabezados de tiempo (ej: "t_1.23s")
# Usamos el vector de tiempo 't' de la simulación
time_headers = [Symbol("t_$(round(time, digits=2))s") for time in t]

# --- 2. Perfiles de Temperatura (n_r filas x (n_t + 1) columnas) ---

# Iniciar el DataFrame con la primera columna: Radio
df_temp = DataFrame()
df_temp.Radius_m = r_span

# Añadir la matriz de temperatura (columna por columna) con los encabezados de tiempo
for (i, header) in enumerate(time_headers)
    df_temp[!, header] = T_matrix[:, i]
end
CSV.write("outputs/xiao19_discharge_temperature_profiles.csv", df_temp)

# --- 3. Perfiles de Adsorción (n_r filas x (n_t + 1) columnas) ---

# Iniciar el DataFrame con la primera columna: Radio
df_adsorption = DataFrame()
df_adsorption.Radius_m = r_span

# Añadir la matriz de adsorción (columna por columna)
for (i, header) in enumerate(time_headers)
    df_adsorption[!, header] = n_a_matrix[:, i]
end
CSV.write("outputs/xiao19_discharge_adsorption_profiles.csv", df_adsorption)

# --- 4. Datos de Series de Tiempo (P, ρ_avg, T_wall) ---
# Aquí es donde T_wall lógicamente pertenece, junto a otras series de tiempo.

df_timeseries = DataFrame(
    Time_s=t,
    Pressure_Pa=P,
    Avg_Density_kg_m3=ρ_avg,
    Wall_Temperature_K=T_wall  # <-- T_wall AÑADIDO AQUÍ
)
CSV.write("outputs/xiao19_discharge_timeseries_data.csv", df_timeseries)

println("Resultados exportados exitosamente.")
df_simulation_avg_temp = DataFrame(Time=t, Temperature=T_avg)
CSV.write("outputs/xiao19_discharge_avg_temperature.csv", df_simulation_avg_temp)