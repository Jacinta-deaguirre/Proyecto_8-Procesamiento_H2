# ldf_sensibilidad_kT.jl - Proyecto 8
# Sensibilidad del llenado al coeficiente LDF k(T) de Arrhenius
# k(T) = k_ref*exp[-(E_a/R)*(1/T - 1/T_ref)], con T_ref = T0 = 281k
# Mismo caso que ldf_sensibilidad_carga.jl (Test 13)

ENV["GKSwstype"] = "100"
include("ldf_cycle_runner.jl")
using Plots

#Parametros del Test 13
prm = build_parameters(13; dr_factor=8.0)
n_r = prm.n_r
T_H2 = 297.6
m_total = 2.023e-5 * 1042.0

#Funcion que simula la carga con coeficiente k y tiempo de llenado t_carga
function correr_carga(k, t_carga; E_a = 0.0)
    m_in = m_total / t_carga
    op = OperationalParameters(prm.U, prm.T_bath, m_in, T_H2, k, E_a, prm.T₀)
    par = AdsorptionParameters(prm.MDA_params, prm.material, prm.geometric, op)
    u₀, du₀, dv = dae_setup(prm.MDA_params, prm.material, prm.geometric, op, prm.Pᵢ, prm.T₀)
    prob = DAEProblem(MDA_adsorption!, du₀, u₀, (0.0, t_carga), par, differential_vars=dv)
    sol=solve(prob, IDA(linear_solver=:LapackDense), reltol=1e-8, abstol=1e-10)
    return sol, m_in
end

#Funcion para ver H2 cargado hasta la presion limite. Busco el instante en que la presion llega al limite y calculo la masa inyectada
function h2_hasta_presion(sol, m_in, P_lim)
    t_carga = sol.t[end]
    t_malla = range(0, t_carga, length=2001)
    P = [sol(t)[2*n_r+2] for t in t_malla]
    i = findfirst(P .>= P_lim)
    if i === nothing
        t_lim = t_carga
    else
        t_lim = t_malla[i-1] + (P_lim - P[i-1])*(t_malla[i] - t_malla[i-1])/(P[i] - P[i-1])
    end
    return t_lim, m_in*t_lim
end

# Ciclo con k(T), para cada tiempo de llenado se corre el equilibrio y E_a = 0 es k constante
t_cargas = [540.0, 180.0] # 9 y 3 minutos
k_valores = [0.03, 4.8e-3] #pellets de 2 y 5cm
E_a_valores = [0.0, 1200.0, 10000.0] #J/mol
tabla = DataFrame()

for t_c in t_cargas
    sol_eq, m_eq = correr_carga(Inf, t_c)
    P_lim_c = sol_eq.u[end][2*n_r+2]
    T_max_eq = maximum(u[1] for u in sol_eq.u)
    for k in k_valores, E_a in E_a_valores
        sol_c, m_c = correr_carga(k, t_c; E_a=E_a)
        t_lim_c, m_carg_c = h2_hasta_presion(sol_c, m_c, P_lim_c)
        T_max = maximum(u[1] for u in sol_c.u)
        f_max= exp(-(E_a/prm.R)*(1/T_max - 1/prm.T₀))
        push!(tabla, (t_carga_s= t_c, k_1_s = k, E_a_kJ_mol= E_a/1000, N = k*t_c,
            P_lim_MPa = P_lim_c / 1e6,
            P_fin_MPa = sol_c.u[end][2*n_r+2]/1e6,
            t_lim_s = t_lim_c,
            H2_cargado_g = m_carg_c * 1000,
            H2_cargado_pct = 100*t_lim_c/t_c,
            T_max_eq_K = T_max_eq,
            T_max_K = T_max,
            dT_max_K = T_max - T_max_eq,
            kTmax_sobre_k = f_max))
        println("t = $t_c s, k = $k, E_a = $(E_a/1000) kJ/mol: H2 = ", round(100*t_lim_c/t_c, digits=2), " %")
    end
end

show(tabla, allrows=true, allcols=true)
println()
CSV.write(joinpath(@__DIR__, "resultados_F4", "ldf_sensibilidad_kT.csv"), tabla)