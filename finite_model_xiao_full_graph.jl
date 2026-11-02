using Revise
using HydrogenAdsorption
using Sundials
using Statistics
using Plots
using CSV
using DataFrames
using Interpolations
using SearchSortedNearest

# --- Load Experimental Data for Full Cycle Plotting ---
pressure_exp_charge_data = CSV.read("inputs/xiao_pressure_finite_element.csv", DataFrame)
pressure_exp_dormant_data = CSV.read("inputs/xiao_pressure_finite_element_dormant.csv", DataFrame)
pressure_exp_discharge_data = CSV.read("inputs/xiao_pressure_finite_element_discharge.csv", DataFrame)

# --- Combine Experimental Data ---
pressure_full_cycle_exp = vcat(pressure_exp_charge_data, pressure_exp_dormant_data, pressure_exp_discharge_data)

# --- Load Time Series Data for Full Cycle Plotting ---
time_series_charge = CSV.read("outputs/charge_timeseries_data.csv", DataFrame)
time_series_dormant = CSV.read("outputs/dormant_timeseries_data.csv", DataFrame)
time_series_discharge = CSV.read("outputs/discharge_timeseries_data.csv", DataFrame)

# --- Combine Time Series Data ---
time_series_full_cycle = vcat(time_series_charge, time_series_dormant, time_series_discharge)
pressure_full_cycle = time_series_full_cycle[!, :Pressure_Pa]
time_full_cycle = time_series_full_cycle[!, :Time_s]

# -----------------------------------------------------
# --- Load Temperature Data for Full Cycle Plotting ---
# -----------------------------------------------------

temperature_charge_data = CSV.read("outputs/charge_temperature_profiles.csv", DataFrame)
temperature_dormant_data = CSV.read("outputs/dormant_temperature_profiles.csv", DataFrame)
temperature_discharge_data = CSV.read("outputs/discharge_temperature_profiles.csv", DataFrame)

# --- Combine Temperature Data ---
temperature_full_cycle = hcat(temperature_charge_data, temperature_dormant_data[:, 3:end], temperature_discharge_data[:, 3:end])

# --- Extract Centre Temperature Data ---
temperature_centre_exp_charge = CSV.read("inputs/xiao_center_temp_finite_element.csv", DataFrame)
temperature_centre_exp_dormant = CSV.read("inputs/xiao_center_temp_finite_element_dormant.csv", DataFrame)
temperature_centre_exp_discharge = CSV.read("inputs/xiao_center_temp_finite_element_discharge.csv", DataFrame)
temperature_centre_exp = vcat(temperature_centre_exp_charge, temperature_centre_exp_dormant, temperature_centre_exp_discharge)

temperature_centre = temperature_full_cycle[1, 2:end]
T_time_vec = []

for str in names(temperature_centre)
    time = strip(str, ['t', 's', '_'])
    push!(T_time_vec, parse(Float64, time))
end

temperature_centre = DataFrame(Time_s=T_time_vec, Temperature_K=collect(temperature_centre))

# --- Extract Middle Node Temperature Data ---
temperature_middle_exp_charge = CSV.read("inputs/xiao_middle_temp_finite_element.csv", DataFrame)
temperature_middle_exp_dormant = CSV.read("inputs/xiao_middle_temp_finite_element_dormant.csv", DataFrame)
temperature_middle_exp_discharge = CSV.read("inputs/xiao_middle_temp_finite_element_discharge.csv", DataFrame)
temperature_middle_exp = vcat(temperature_middle_exp_charge, temperature_middle_exp_dormant, temperature_middle_exp_discharge)

middle_index = div(length(temperature_full_cycle[!, 1]), 2) # Index of the middle node
temperature_middle = temperature_full_cycle[middle_index, 2:end]
temperature_middle = DataFrame(Time_s=T_time_vec, Temperature_K=collect(temperature_middle))

# --- Extract Wall Temperature Data ---
temperature_wall_exp_charge = CSV.read("inputs/xiao_wall_temp_finite_element.csv", DataFrame)
temperature_wall_exp_dormant = CSV.read("inputs/xiao_wall_temp_finite_element_dormant.csv", DataFrame)
temperature_wall_exp_discharge = CSV.read("inputs/xiao_wall_temp_finite_element_discharge.csv", DataFrame)
temperature_wall_exp = vcat(temperature_wall_exp_charge, temperature_wall_exp_dormant, temperature_wall_exp_discharge)

temperature_wall = CSV.read("inputs/xiao_wall_temp_finite_element.csv", DataFrame)
temperature_wall = DataFrame(Time_s=time_series_full_cycle.Time_s, Temperature_K=time_series_full_cycle.Wall_Temperature_K)


# -----------------------------------------------------
# --- Load Mass Data for Full Cycle Plotting ---
# -----------------------------------------------------
mass_charge_data = CSV.read("outputs/charge_mass_profiles.csv", DataFrame)
mass_dormant_data = CSV.read("outputs/dormant_mass_profiles.csv", DataFrame)
mass_discharge_data = CSV.read("outputs/discharge_mass_profiles.csv", DataFrame)
# --- Combine Mass Data ---
mass_full_cycle = vcat(mass_charge_data, mass_dormant_data, mass_discharge_data)

# Load Experimental Mass Data
mass_charge_data = CSV.read("inputs/xiao_total_mass_finite_element.csv", DataFrame)
mass_dormant_data = CSV.read("inputs/xiao_total_mass_finite_element_dormant.csv", DataFrame)
mass_discharge_data = CSV.read("inputs/xiao_total_mass_finite_element_discharge.csv", DataFrame)
mass_exp_full_cycle = vcat(mass_charge_data, mass_dormant_data, mass_discharge_data)

ads_charge_data = CSV.read("inputs/xiao_adsorbed_mass_finite_element.csv", DataFrame)
ads_dormant_data = CSV.read("inputs/xiao_adsorbed_mass_finite_element_dormant.csv", DataFrame)
ads_discharge_data = CSV.read("inputs/xiao_adsorbed_mass_finite_element_discharge.csv", DataFrame)
ads_exp_full_cycle = vcat(ads_charge_data, ads_dormant_data, ads_discharge_data)

gass_charge_data = CSV.read("inputs/xiao_gas_mass_finite_element.csv", DataFrame)
gass_dormant_data = CSV.read("inputs/xiao_gas_mass_finite_element_dormant.csv", DataFrame)
gass_discharge_data = CSV.read("inputs/xiao_gas_mass_finite_element_discharge.csv", DataFrame)
gas_exp_full_cycle = vcat(gass_charge_data, gass_dormant_data, gass_discharge_data)


# --- 1. Global Plot Styling ---
# Apply a theme with larger fonts, similar to the Python 'seaborn-talk' style
theme(:default,
    fontfamily="sans-serif",
    titlefontsize=22,
    guidefontsize=20,     # X/Y axis labels
    tickfontsize=16,
    legendfontsize=10,
    margin=0.75Plots.cm,   # Add generous margins
    grid=true,
    gridstyle=:dash,
    gridalpha=0.5,
    framestyle=:box,      # Adds a full box, looks professional
    palette=:seaborn_pastel,   # Use a modern, colorblind-friendly palette
)

# Define common line/marker styles for clarity
# (Simulated = thick solid, Experimental = thinner dash + markers)
sim_lw = 2.5
exp_lw = 0
exp_marker = (:circle, 6, 0.7) # (shape, size, alpha)

# Plot A: Temperature
p_temp = plot(
    ylabel="Temperature / K",
    xlabel="Time / s",
    legend=:bottomleft, # Move legend inside
    legend_columns=1,     # Arrange in 3 columns (2 rows)
    size=(1200, 900)
)

# Center temperature
plot!(p_temp, temperature_centre.Time_s, temperature_centre.Temperature_K,
    label="Sim. Center", color="#D62728", lw=sim_lw)
plot!(p_temp, temperature_centre_exp.Time, temperature_centre_exp.Temperature,
    label="Exp. Center", color="#D62728", marker=exp_marker, linestyle=:dash, lw=exp_lw)

# Middle node temperature
plot!(p_temp, temperature_middle.Time_s, temperature_middle.Temperature_K,
    label="Sim. Middle", color="#9467BD", lw=sim_lw)
plot!(p_temp, temperature_middle_exp.Time, temperature_middle_exp.Temperature,
    label="Exp. Middle", color="#9467BD", marker=exp_marker, linestyle=:dash, lw=exp_lw)

plot!(p_temp, temperature_wall.Time_s, temperature_wall.Temperature_K,
    label="Sim. Wall", color="#808080", lw=sim_lw)
plot!(p_temp, temperature_wall_exp.Time, temperature_wall_exp.Temperature,
    label="Exp. Wall", color="#808080", marker=exp_marker, linestyle=:dash, lw=exp_lw)

#display(p_temp)

# Plot B: Pressure
p_pressure = plot(
    ylabel="Pressure / MPa",
    xlabel="Time / s",
    legend=:bottomright, # Move legend inside
    legend_columns=1,     # Arrange in 2 columns (1 row)
    size=(1200, 900)
)

plot!(p_pressure, time_series_full_cycle[!, :Time_s], pressure_full_cycle ./ 1e6,
    label="Simulated", color=1, lw=sim_lw)
plot!(p_pressure, pressure_full_cycle_exp.Time, pressure_full_cycle_exp.Pressure,
    label="Experimental", color=1, marker=exp_marker, linestyle=:dash, lw=exp_lw)

display(p_pressure)

# Plot C: Mass Profiles
p_mass = plot(
    ylabel="Mass / g",
    xlabel="Time / s",
    legend=:best, # Move legend inside
    legend_columns=1,     # Arrange in 3 columns (2 rows)
    size=(1200, 900)
)

# Gas mass
plot!(p_mass, gas_exp_full_cycle.Time, gas_exp_full_cycle.Mass .* 1000,
    label="Exp. Gas", color=1, marker=exp_marker, linestyle=:dash, lw=exp_lw)
plot!(p_mass, mass_full_cycle.Time_s, mass_full_cycle.Gas_Mass_kg .* 1000,
    label="Sim. Gas", color=1, lw=sim_lw)

# Adsorbed mass
plot!(p_mass, ads_exp_full_cycle.Time, ads_exp_full_cycle.Mass .* 1000,
    label="Exp. Adsorbed", color=2, marker=exp_marker, linestyle=:dash, lw=exp_lw)
plot!(p_mass, mass_full_cycle.Time_s, mass_full_cycle.Adsorbed_Mass_kg .* 1000,
    label="Sim. Adsorbed", color=2, lw=sim_lw)

# Total mass
plot!(p_mass, mass_exp_full_cycle.Time, mass_exp_full_cycle.Mass .* 1000,
    label="Exp. Total", color=3, marker=exp_marker, linestyle=:dash, lw=exp_lw)
plot!(p_mass, mass_full_cycle.Time_s, mass_full_cycle.Total_Mass_kg .* 1000,
    label="Sim. Total", color=3, lw=sim_lw)

#display(p_mass)

# --- 3. Combine Plots into a Single Figure ---
# Use a 2x2 layout, placing Temp (a) and Pressure (b) side-by-side,
# and Mass (c) spanning the full width below them.
l = @layout [a b; c{0.5h}] # Bottom plot (c) takes 50% of the total height
final_plot = plot(p_temp, p_pressure, p_mass,
    layout=l,
    size=(1800, 1400) # A large, wide figure suitable for presentations
)

# --- 4. Save the Figure ---
# Save as SVG (vector format) for maximum quality
output_filename = "outputs/full_simulation_profiles.svg"
savefig(final_plot, output_filename)
println("\nPlot saved successfully as '$output_filename'")

# Save individuals
savefig(p_temp, "outputs/full_julia_temperature_profile.svg")
savefig(p_pressure, "outputs/full_julia_pressure_profile.svg")
savefig(p_mass, "outputs/full_julia_mass_profiles.svg")

# Save individuals
savefig(p_temp, "outputs/full_julia_temperature_profile.png")
savefig(p_pressure, "outputs/full_julia_pressure_profile.png")
savefig(p_mass, "outputs/full_julia_mass_profiles.png")

# Calculate performance metrics
# average percentage absolute deviation (AAD)
display(time_series_full_cycle)
display(pressure_full_cycle_exp)



