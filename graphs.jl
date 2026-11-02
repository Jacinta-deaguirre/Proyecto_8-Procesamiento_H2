using CSV
using DataFrames
using Plots

df_mda = CSV.read("xiao_mda_simulation.csv", DataFrame)
df_clausius = CSV.read("xiao_clausius_simulation.csv", DataFrame)
df_heat_exp = CSV.read("xiao_heat_exp.csv", DataFrame)

plt = plot(layout=(3, 1), size=(900, 900), margin=Plots.cm)
plot!(plt[1], df_mda.t, df_mda.dH, label="MDA Isosteric Heat", xlabel="Time (s)", ylabel="Isosteric Heat (J/mol)", color=:blue)
plot!(plt[1], df_clausius.t, df_clausius.dH, label="Clausius Isosteric Heat", color=:red)

plot!(plt[2], df_mda.t, df_mda.T_avg, label="MDA Average Temperature", xlabel="Time (s)", ylabel="Average Temperature (K)", color=:blue)
plot!(plt[2], df_clausius.t, df_clausius.T_avg, label="Clausius Average Temperature", color=:red)

plot!(plt[3], df_mda.t, df_mda.P, label="MDA Pressure", xlabel="Time (s)", ylabel="Pressure (MPa)", color=:blue)
plot!(plt[3], df_clausius.t, df_clausius.P, label="Clausius Pressure", color=:red)

scatter!(plt[1], df_heat_exp[:, 1], df_heat_exp[:, 2]*1000, label="Experimental Data", color=:green, markersize=4)

# Temperature data
t_exp = [18.633540372670836, 204.96894409937886, 409.9378881987577, 605.5900621118014, 805.9006211180124, 1006.2111801242236]
T_exp = [281.0010449320794, 285.6530825496343, 288.43260188087777, 290.15882967607104, 291.1243469174504, 291.6802507836991]
# Pressure data
t_exp2 = [0, 189.8370086, 402.6845638, 598.274209, 799.6164909, 1000.958773]
p_exp = [0.102564103, 1.487179487, 3.115384615, 4.858974359, 6.666666667, 8.58974359]

scatter!(plt[2], t_exp, T_exp, label="Experimental Temperature vs Time", legend=false, markersize=4, color=:green)
scatter!(plt[3], t_exp2, p_exp, label="Experimental Pressure vs Time", legend=false, markersize=4, color=:green)