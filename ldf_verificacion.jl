# ldf_verification.jl - Proyecto 8
# Que hace: simula la etapa de carga del test 13 con los mismos
#parámetros de finite_model_xiao.jl pero para tres valores de k.
# k = inf : debe dar el modelo base
# k = 300 1/s : debe dar muy similar al modelo base
# k = 3 1/s : LDF de un pellet de 2mm.

ENV["GKSwstype"] = "100" # headless plotting
using HydrogenAdsorption
using Sundials
using Statistics
using Plots
using CSV
using DataFrames
using Printf

# parameters
# Tank parameters
V = 2.4946e-3 # Bed volume / m³ (measured, Richard et al. 2009)
R_T = 0.0469 # Bed radius / m (physical internal radius, Xiao et al. 2012)
L = V / (pi * R_T^2) # Equivalent cylindrical length of the bed / m (= 0.361 m; the internal height of 0.45 m includes the non-cylindrical ends)
e = 0.0508 - 0.0469 # Wall thickness / m
c_wall = 468 # Specific heat capacity of the tank wall / J/kg K
m_wall = 3.714 # Mass of the tank wall / kg
k_wall = 13 # Thermal conductivity of the tank wall / W/m K

# Inputs needed for the coefficient_matrix function
dr = (1 / 2)^2 * 0.0025 # Radial step size / m
U = 36 * 1 # Heat transfer coefficient / W/(m²·K)

# Tank surface areas
Ai = 2 * pi * R_T * L # Inner surface area of the tank / m²
Ao = 2 * pi * (R_T + e) * L # Outer surface area of the tank / m²

T₀ = 281.0 # Initial temperature of the tank / K
#T_air = 281.0 # Ambient temperature / K
T_air = 282.5 # Ambient temperature / K

#T_H2 = 281.6 # Temperature of the incoming hydrogen gas / K
T_H2 = 297.6 # Data from finite element model / K


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
# mₛ (mass of activated carbon) is set below as ρₛ (1 - ε_b) V, consistent with Eqs. (2)-(3)
kₛ = 0.646 # Thermal conductivity of activated carbon / W/m K
ε_b = 0.49 # Bed porosity
mₛ = ρₛ * (1 - ε_b) * V # Mass of activated carbon implied by the bed / kg

# Hydrogen properties
cₚ = 14700.0 # Specific heat capacity of hydrogen / J/kg K
cᵥ = 10134.0 # Specific heat capacity at constant volume / J/kg K
M_H2 = 2.0159e-3 # Molar mass of hydrogen / kg/mol
R = 8.314 # Ideal gas constant / J/mol K
k_g = 0.206 # Thermal conductivity of hydrogen / W/m K

k_eff = (kₛ * (1 - ε_b) + k_g * ε_b) # Effective thermal conductivity / W/(m·K)
n_r, r_span, A, b = coefficient_matrix(R_T, dr, k_eff, U, T_air)
r_span = collect(r_span)
dr = r_span[2] - r_span[1] # Actual grid spacing (grid covers R_T exactly)
material_props = MaterialProperties(ρₛ, cₛ, mₛ, kₛ, ε_b, cₚ, cᵥ, M_H2, R, k_g, k_eff)

# Geometric parameters
geometric_params = GeometricParameters(n_r, dr, V, L, A, b, r_span, R_T, e, Ao, Ai, c_wall, m_wall, k_wall)

# Operational Parameters
m_in = 1 * 2.023e-5 # Mass flow rate of hydrogen / kg / s

Pᵢ = 0.033e6 # Initial pressure / Pa
t_f = 1042.0 # Final time // s

function correr_carga(k; reltol=1e-3, abstol=1e-6)
    #fija la operación con 5 valores
    op = OperationalParameters(U, T_air, m_in, T_H2, k) 
    #Ficha del sistema con isoterma, material, geometría y operacion
    par = AdsorptionParameters(MDA_params, material_props, geometric_params, op) 
    #Condicion inicial
    u₀, du₀, dv = dae_setup(MDA_params, material_props, geometric_params, op, Pᵢ, T₀) 
    #Problema que resuelve
    prob = DAEProblem(MDA_adsorption!, du₀, u₀, (0.0,t_f), par, differential_vars=dv)
    #Resolver con IDA
    t_inicio = time()
    sol = solve(prob, IDA(linear_solver=:LapackDense), reltol=reltol, abstol=abstol, progress=true)
    tiempo = time() - t_inicio
    return sol, tiempo
end

function balance_masa(sol)
    #Sacar p y n_a de la solución
    ρ_avg = [u[2*n_r+1] for u in sol.u]
    n_avg = [cyl_avg(u[n_r+1:2*n_r], r_span, R_T) for u in sol.u]

    #masas de H2 
    m_H2_gas = ρ_avg .* V * ε_b # Mass of hydrogen in the gas phase / kg
    m_H2_ads = n_avg .* mₛ * M_H2 # Mass of hydrogen in the adsorbed phase / kg (mₛ = ρₛ (1 - ε_b) V)
    m_H2_total = m_H2_gas .+ m_H2_ads # Total mass of hydrogen in the tank / kg 

    #masa esperada y error
    m_esp = m_H2_total[1] .+ m_in .* sol.t 
    err_final = 100 * (m_H2_total[end] - m_esp[end])/m_esp[end]
    err_max = 100 * maximum(abs.(m_H2_total .- m_esp)) / m_esp[end]

    return m_H2_gas, m_H2_ads, m_H2_total, err_final, err_max
end

u_base = CSV.read(joinpath(@__DIR__, "corrida_prueba_base2", "charge_final_state.csv"), DataFrame).u 
sol_def, tiempo_def = correr_carga(Inf)
dif_identidad = maximum(abs.(sol_def.u[end] .- u_base) ./ abs.(u_base))
println("Prueba 1 (k = Inf v/s corrida base guardada) : mayor diferencia relativa =", dif_identidad)

k_valores = [Inf, 300.0, 3.0]
i_medio = div(n_r, 2)
tabla = DataFrame()
soluciones = Dict()

for k in k_valores
    sol_k, tiempo_k = correr_carga(k, reltol = 1e-8, abstol =1e-10)
    m_gas_k, m_ads_k, m_tot_k, err_final_k, err_max_k = balance_masa(sol_k)
    u_f = sol_k.u[end]
    push!(tabla, (k_1_s=k, tiempo_s = tiempo_k, P_MPa = u_f[2*n_r+2]/1e6, T_centro_K = u_f[1], T_medio_K = u_f[i_medio], 
                    m_ads_g = m_ads_k[end]*1000, m_total_g = m_tot_k[end]*1000, err_masa_final_pct = err_final_k, err_masa_max_pct = err_max_k))
    soluciones[k] = sol_k
end

tabla.dif_t_centro_K = tabla.T_centro_K .- tabla.T_centro_K[1]
tabla.dif_P_pct = 100 .* (tabla.P_MPa .- tabla.P_MPa[1]) ./ tabla.P_MPa[1]

show(tabla, allcols=true)
println()
mkpath(joinpath(@__DIR__, "outputs"))
CSV.write(joinpath(@__DIR__, "outputs", "ldf_verificacion.csv"), tabla)

estilos = Dict(Inf => :solid, 300.0 => :dash, 3.0 => :dot)
p1 = plot(xlabel="Tiempo / s", ylabel="Presión / Mpa", legend=:bottomright)
p2 = plot(xlabel="Tiempo / s", ylabel="T centro / K", legend=false)
p3 = plot(xlabel="Tiempo / s", ylabel="Masa adsorbida / g", legend=false)

for k in k_valores
    s = soluciones[k]
    m_ads_k = balance_masa(s)[2]
    etiqueta = isinf(k) ? "k = inf (base)" : "k=$(k) 1/s"
    plot!(p1, s.t, [u[2*n_r+2] / 1e6 for u in s.u], label=etiqueta, linestyle=estilos[k], linewidth=2)
    plot!(p2, s.t, [u[1] for u in s.u], label=etiqueta, linestyle=estilos[k], linewidth=2)
    plot!(p3, s.t, m_ads_k .* 1000, label=etiqueta, linestyle=estilos[k], linewidth=2)
end

fig = plot(p1, p2, p3, layout=(1, 3), size=(1500, 450), margin=5Plots.mm)
savefig(fig, joinpath(@__DIR__, "outputs", "ldf_verificacion.svg"))
display(fig)