#ldf_curva_diseno.jl - Proyecto 8, Fase 4
#Curva de diseno - pellet maximo para cada tiempo de llenado con N=k*t ≥ N_crit
# k = 15*D_eff/R_p^2 -> d_max = 2*raiz(15*D_eff*t/N_crit)

ENV["GKSwstype"] = "100"
include("ldf_cycle_runner.jl")

#Parametros del Test 13
prm = build_parameters(13; dr_factor=8.0)
n_r = prm.n_r
T_H2 = 297.6
m_total = 2.023e-5 * 1042.0
D_eff = 2e-7

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

# Ciclo con k(T), para cada tiempo de llenado se corre el equilibrio y E_a = 0 es k constante
t_cargas = [1042.0, 540.0, 300.0, 180.0] # 17, 9, 5 y 3 minutos
N_valores = [50.0, 100.0] #criterio a verificar
tabla = DataFrame()

for t_c in t_cargas
    sol_eq, m_eq = correr_carga(Inf, t_c)
    P_lim_c = sol_eq.u[end][2*n_r+2]
    for N in N_valores
        k = N/t_c
        d_max_mm = 2*sqrt(15*D_eff / k)*1000
        sol_c, m_c = correr_carga(k, t_c)
        t_lim_c, m_carg_c = h2_hasta_presion(sol_c, m_c, P_lim_c)
        push!(tabla, (t_carga_s= t_c, N = N, k_1_s = k, d_max_mm = d_max_mm,
            P_lim_MPa = P_lim_c / 1e6,
            H2_cargado_g = m_carg_c * 1000,
            H2_cargado_pct = 100*t_lim_c/t_c))
        println("t = $t_c s, N = $N: k = ", round(k, digits=4), " 1/s, d_max = ",
                round(d_max_mm, digits=1), " mm, H2 = ", round(100*t_lim_c/t_c, digits=2), " %")
    end
end

show(tabla, allrows=true, allcols=true)
println()
CSV.write(joinpath(@__DIR__, "resultados_F4", "ldf_curva_diseno.csv"), tabla)