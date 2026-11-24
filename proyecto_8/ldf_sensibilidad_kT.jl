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

#Funcion para ver H2 cargado hasta la presion limite. Busco el el nstante en que la presion llega al limite y calculo la masa inyectada
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

#Verificación con E_a=0 se deben tener los mismos resultados de resultados_F4
ref = CSV.read(joinpath(@__DIR__, "resultados_F4", "ldf_sensibilidad_carga.csv"), DataFrame)
for (k, t_c) in [(0.03, 540.0), (4.8e-3, 180.0)]
    sol_eq, _=correr_carga(Inf, t_c)
    P_lim = sol_eq.u[end][2*n_r+2]
    sol, m_in = correr_carga(k, t_c; E_a=.0)
    t_lim, _=h2_hasta_presion(sol, m_in, P_lim)
    fila = ref[(ref.t_carga_s .==t_c) .& (ref.k_1_s .== k), :]
    println("k=$k, t = $t_c s")
    println(" P_lim: ahora", P_lim/1e6, "MPa | antes", fila.P_lim_MPa[1])
    println(" H2 cargado ahora:", 100*t_lim/t_c, " %| antes", fila.H2_cargado_pct[1])
end