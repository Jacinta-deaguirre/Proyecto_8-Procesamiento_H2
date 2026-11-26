#ldf_fig1_carga.jl - Proyecto 8, Fase 4
#Carga completa con equilibrio y los 5 valores de k, para llenados de 17 min (test 13) y 9min (CR1)

ENV["GKSwstype"] = "100"
include("ldf_cycle_runner.jl")

#Parametros del Test 13
prm = build_parameters(13; dr_factor=8.0)
n_r = prm.n_r
T_H2 = 297.6
m_total = 2.023e-5 * 1042.0

#Funcion que simula la carga con coeficiente k y tiempo de llenado t_carga
function correr_carga(k, t_carga, E_a = 0.0)
    m_in = m_total / t_carga
    op = OperationalParameters(prm.U, prm.T_bath, m_in, T_H2, k, E_a, prm.T₀)
    par = AdsorptionParameters(prm.MDA_params, prm.material, prm.geometric, op)
    u₀, du₀, dv = dae_setup(prm.MDA_params, prm.material, prm.geometric, op, prm.Pᵢ, prm.T₀)
    prob = DAEProblem(MDA_adsorption!, du₀, u₀, (0.0, t_carga), par, differential_vars=dv)
    sol=solve(prob, IDA(linear_solver=:LapackDense), reltol=1e-8, abstol=1e-10)
    return sol, m_in
end

#Casos
t_cargas = [1042.0, 540.0] #17 min (Test 13) y 9 min (CR1)
casos = [(Inf, "Equilibrio"), (3.0, "2 mm"), (0.12, "1 cm"), (0.03, "2 cm" ), (4.8e-3, "5 cm"), (1.2e-3, "10 cm")]
datos=DataFrame()

for t_c in t_cargas, (k, caso) in casos
    sol, m_in = correr_carga(k, t_c)
    for t in range(0, t_c, length=301)
        u=sol(t)
        n_prom= cyl_avg(u[n_r+1:2*n_r], prm.r_span, prm.R_T)
        push!(datos, (t_carga_s = t_c, k_1_s = k, caso = caso, t_s = t,
            P_MPa = u[2*n_r+2]/1e6,
            T_centro_k = u[1],
            m_ads_g = n_prom * prm.material.mₛ * prm.M_H2*1000))
    end
    println("Listo: t = $t_c s, $caso")
end

CSV.write(joinpath(@__DIR__, "resultados_F4", "ldf_fig1_carga.csv"), datos)
