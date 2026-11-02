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
pressure_full_cycle_exp = CSV.read("inputs/xiao19_pressure.csv", DataFrame)
# --- Load Time Series Data for Full Cycle Plotting ---
time_series_charge = CSV.read("outputs/xiao19_charge_timeseries_data.csv", DataFrame)
time_series_dormant = CSV.read("outputs/xiao19_dormant_timeseries_data.csv", DataFrame)
time_series_discharge = CSV.read("outputs/xiao19_discharge_timeseries_data.csv", DataFrame)

# --- Combine Time Series Data ---
time_series_full_cycle = vcat(time_series_charge, time_series_dormant, time_series_discharge)
pressure_full_cycle = time_series_full_cycle[!, :Pressure_Pa]
time_full_cycle = time_series_full_cycle[!, :Time_s]

# -----------------------------------------------------
# --- Load Temperature Data for Full Cycle Plotting ---
# -----------------------------------------------------

temperature_charge_data = CSV.read("outputs/xiao19_charge_avg_temperature.csv", DataFrame)
temperature_dormant_data = CSV.read("outputs/xiao19_dormant_avg_temperature.csv", DataFrame)
temperature_dormant_data = temperature_dormant_data[2:end, :] # Initial condition
temperature_discharge_data = CSV.read("outputs/xiao19_discharge_avg_temperature.csv", DataFrame)
temperature_discharge_data = temperature_discharge_data[2:end, :] # Initial condition

# --- Combine Temperature Data ---
temperature_full_cycle = vcat(temperature_charge_data, temperature_dormant_data, temperature_discharge_data)

# --- Load Experimental Temperature Data ---
temperature_charge = CSV.read("inputs/xiao19_avg_temperature.csv", DataFrame)
temperature_dormant = CSV.read("inputs/xiao19_avg_temperature_dormant.csv", DataFrame)
temperature_discharge = CSV.read("inputs/xiao19_avg_temperature_discharge.csv", DataFrame)

# --- Combine Experimental Temperature Data ---
temperature_exp_full_cycle = vcat(temperature_charge, temperature_dormant, temperature_discharge)

# -----------------------------------------------------
# --- Load Mass Data for Full Cycle Plotting ---
# -----------------------------------------------------
mass_charge_data = CSV.read("outputs/xiao19_charge_mass_profiles.csv", DataFrame)
mass_dormant_data = CSV.read("outputs/xiao19_dormant_mass_profiles.csv", DataFrame)
mass_discharge_data = CSV.read("outputs/xiao19_discharge_mass_profiles.csv", DataFrame)
# --- Combine Mass Data ---
mass_full_cycle = vcat(mass_charge_data, mass_dormant_data, mass_discharge_data)

# --- Load Experimental Mass Data ---
mass_exp = CSV.read("inputs/xiao19_total_mass.csv", DataFrame)
mass_gas_exp = CSV.read("inputs/xiao19_gas_mass.csv", DataFrame)
mass_adsorbed_exp = CSV.read("inputs/xiao19_ads_mass.csv", DataFrame)

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
    legend=:bottomright, # Move legend inside
    legend_columns=1,     # Arrange in 3 columns (2 rows)
    size=(1200, 900)
)

plot!(p_temp, temperature_exp_full_cycle.Time, temperature_exp_full_cycle.Temperature,
    label="Exp. Average", color=2, marker=exp_marker, linestyle=:dash, lw=exp_lw)
plot!(p_temp, temperature_full_cycle.Time, temperature_full_cycle.Temperature,
    label="Sim. Average", color=2, lw=sim_lw)


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


# Plot C: Mass Profiles
p_mass = plot(
    ylabel="Mass / g",
    xlabel="Time / s",
    legend=:best, # Move legend inside
    legend_columns=1,     # Arrange in 3 columns (2 rows)
    size=(1200, 900)
)

# Gas mass
plot!(p_mass, mass_gas_exp.Time, mass_gas_exp.Mass .* 1000,
    label="Exp. Gas", color=1, marker=exp_marker, linestyle=:dash, lw=exp_lw)
plot!(p_mass, mass_full_cycle.Time_s, mass_full_cycle.Gas_Mass_kg .* 1000,
    label="Sim. Gas", color=1, lw=sim_lw)

# Adsorbed mass
plot!(p_mass, mass_adsorbed_exp.Time, mass_adsorbed_exp.Mass .* 1000,
    label="Exp. Adsorbed", color=2, marker=exp_marker, linestyle=:dash, lw=exp_lw)
plot!(p_mass, mass_full_cycle.Time_s, mass_full_cycle.Adsorbed_Mass_kg .* 1000,
    label="Sim. Adsorbed", color=2, lw=sim_lw)

# Total mass
plot!(p_mass, mass_exp.Time, mass_exp.Mass .* 1000,
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

# Display the final combined plot
display(final_plot)

# --- 4. Save the Figure ---
# Save as SVG (vector format) for maximum quality
output_filename = "outputs/xiao19_full_simulation_profiles.svg"
savefig(final_plot, output_filename)
println("\nPlot saved successfully as '$output_filename'")

# Save individuals
savefig(p_temp, "outputs/xiao19_full_julia_temperature_profile.svg")
savefig(p_pressure, "outputs/xiao19_full_julia_pressure_profile.svg")
savefig(p_mass, "outputs/xiao19_full_julia_mass_profiles.svg")

# Save individuals
savefig(p_temp, "outputs/xiao19_full_julia_temperature_profile.png")
savefig(p_pressure, "outputs/xiao19_full_julia_pressure_profile.png")
savefig(p_mass, "outputs/xiao19_full_julia_mass_profiles.png")

# Calculate performance metrics
# average percentage absolute deviation (AAD)
display(time_series_full_cycle)
display(pressure_full_cycle_exp)




