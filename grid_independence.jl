using DataFrames
using CSV
using Plots
using Interpolations
using LinearAlgebra

function number_of_nodes(dr)
    n_r = floor(Int32, R_T / dr) + 1
    return n_r
end

function radial_range(n_r, R_T)
    return range(0, stop=R_T, length=n_r) # Generates radial nodes // m
end

df_T = DataFrame(CSV.File("outputs/grid_independence/temperature_profile_normal_mesh.csv"))
df_T_fine = DataFrame(CSV.File("outputs/grid_independence/temperature_profile_fine_mesh.csv"))
df_T_coarse = DataFrame(CSV.File("outputs/grid_independence/temperature_profile_coarse_mesh.csv"))

dr_normal = (1 / 2)^5 * 0.0025 # Radial step size for normal mesh / m
dr_fine = (1 / 2)^6 * 0.0025 # Radial step size for fine mesh / m
dr_coarse = (1 / 2)^4 * 0.0025 # Rad


# Number of nodes
n_r_normal = number_of_nodes(dr_normal)
n_r_fine = number_of_nodes(dr_fine)
n_r_coarse = number_of_nodes(dr_coarse)
println("Number of nodes in normal mesh: ", n_r_normal)
println("Number of nodes in fine mesh: ", n_r_fine)
println("Number of nodes in coarse mesh: ", n_r_coarse)


# Radial ranges
r_span_normal = radial_range(n_r_normal, R_T)
r_span_fine = radial_range(n_r_fine, R_T)
r_span_coarse = radial_range(n_r_coarse, R_T)

#show(names(df_T))
#show(names(df_T_fine))
#show(names(df_T_coarse))


plot(r_span_normal, df_T[!, end], label="Normal mesh", xlabel="Radius (m)", ylabel="Temperature (K)")
plot!(r_span_fine, df_T_fine[!, end], label="Fine mesh", xlabel="Radius (m)", ylabel="Temperature (K)")
plot!(r_span_coarse, df_T_coarse[!, end], label="Coarse mesh", xlabel="Radius (m)", ylabel="Temperature (K)")

# Create the interpolation
# Normal mesh
itp_T = interpolate(df_T[!, end], BSpline(Cubic(Line(OnGrid()))))
sitp_T = scale(itp_T, r_span_normal)

# Fine mesh 
itp_T_fine = interpolate(df_T_fine[!, end], BSpline(Cubic(Line(OnGrid()))))
sitp_T_fine = scale(itp_T_fine, r_span_fine)

# Coarse mesh
itp_T_coarse = interpolate(df_T_coarse[!, end], BSpline(Cubic(Line(OnGrid()))))
sitp_T_coarse = scale(itp_T_coarse, r_span_coarse)

# Interpolate coarse mesh temperature profile onto the normal mesh radial points
T_coarse_on_normal = sitp_T_coarse.(r_span_normal)
T_fine_on_normal = sitp_T_fine.(r_span_normal)


# Examples evaluating the spline and comparing with the original values
plt = plot()
plot!(r_span_normal, df_T[!, end], label="Normal mesh", xlabel="Radius (m)", ylabel="Temperature (K)")
plot!(r_span_normal, sitp_T(r_span_normal), label="Interpolation", xlabel="Radius (m)", ylabel="Temperature (K)")
plot!(r_span_fine, df_T_fine[!, end], label="Fine mesh", xlabel="Radius (m)", ylabel="Temperature (K)")
plot!(r_span_normal, sitp_T_fine.(r_span_normal), label="Interpolation", xlabel="Radius (m)", ylabel="Temperature (K)")

#plot!(r_span_coarse, df_T_coarse[!, end], label="Coarse mesh", xlabel="Radius (m)", ylabel="Temperature (K)")
#plot!(r_span_normal, sitp_T_coarse.(r_span_normal), label="Interpolation", xlabel="Radius (m)", ylabel="Temperature (K)")
display(plt)

err = norm(sitp_T_fine.(r_span_normal) .- df_T[!, end]) / n_r_normal
println("Average Porcentual Error between normal and fine mesh: ", err * 100, "%")
err_coarse = norm(sitp_T_coarse.(r_span_normal) .- df_T[!, end]) / n_r_normal
println("Average Porcentual Error between normal and coarse mesh: ", err_coarse * 100, "%")