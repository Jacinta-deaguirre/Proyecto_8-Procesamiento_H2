# ldf_sensibilidad_carga.jl - Proyecto 8 (IIQ3843, modelo LDF)
# Que hace: sensibilidad del llenado al coeficiente LDF k y al tiempo de llenado t_carga.
# Mismo estanque y masa que el Test 13 (21,08g), con cuatro tiempos de llenado:
# 1042s (Test 13), 540s (mencionado CR1), 300s y 180s, y 5 valores de k (pellets de 2mm a 10cm): 20 corridas más 4 en equilibrio
# Para cada caso calcula el H2 cargado hasta la presion limite (la final del caso en equilibrio).
# Como correr: en PowerShell poner julia --project=. proyecto_8/ldf_sensibilidad_carga.jl 
# Salida: en proyecto_8/resultados_F4/ldf_sensibilidad_carga.csv (despues la leen fig2_h2_vs_N.jl, tabla1_h2_cargado.jl y fig_ldf_kT.jl)

include("ldf_cycle_runner.jl")
carpeta = joinpath(@__DIR__, "resultados_F4")

#Parametros del Test 13
prm = build_parameters(13; dr_factor=8.0)
n_r = prm.n_r
T_H2 = 297.6
m_total = 2.023e-5 * 1042.0

#Funcion que simula la carga con coeficiente k y tiempo de llenado t_carga
function correr_carga(k, t_carga)
    m_in = m_total / t_carga
    op = OperationalParameters(prm.U, prm.T_bath, m_in, T_H2, k)
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

# Ciclo, para cada tiempo de llenado se corre el equilibrio
t_cargas = [1042.0, 540.0, 300.0, 180.0]
k_valores = [3.0, 0.12, 0.03, 4.8e-3, 1.2e-3]
tabla = DataFrame()

for t_c in t_cargas
    sol_eq, m_eq = correr_carga(Inf, t_c)
    P_lim_c = sol_eq.u[end][2*n_r+2]
    T_max_eq = maximum(u[1] for u in sol_eq.u)
    for k in k_valores
        sol_c, m_c = correr_carga(k, t_c)
        t_lim_c, m_carg_c = h2_hasta_presion(sol_c, m_c, P_lim_c)
        T_max = maximum(u[1] for u in sol_c.u)
        push!(tabla, (t_carga_s= t_c, k_1_s = k, N = k*t_c,
            P_lim_MPa = P_lim_c / 1e6,
            P_fin_MPa = sol_c.u[end][2*n_r+2]/1e6,
            t_lim_s = t_lim_c,
            H2_cargado_g = m_carg_c * 1000,
            H2_cargado_pct = 100*t_lim_c/t_c,
            T_max_eq_K = T_max_eq,
            T_max_K = T_max,
            dT_max_K = T_max - T_max_eq))
    end
end

show(tabla, allrows=true, allcols=true)
println()
mkpath(carpeta)
CSV.write(joinpath(carpeta, "ldf_sensibilidad_carga.csv"), tabla)

mkpath(carpeta)
CSV.write(joinpath(carpeta, "ldf_sensibilidad_carga.csv"), tabla)