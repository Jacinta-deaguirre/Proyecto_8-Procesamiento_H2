# ldf_grid_study.jl - Proyecto 8 (IIQ3843, modelo LDF)
#Que hace: convergencia de malla del modelo LDF en el ciclo completo del Test 13 (carga, almacenamiento y descarga), 
#   con tolerancias finas del solver para que las diferencias entre mallas midan solo el error espacial.
# Basado en grid_study.jl
# Mallas: dr_factor 8, 4 y 2 (76, 151 y 301 nodos), para k = 3 y 0,03 1/s
# Como correr: en PowerShell poner julia --project=. proyecto_8/ldf_grid_study.jl
# Salida: en proyecto_8/verificacion_F3/ -> ldf_grid_study_k3.csv y ldf_grid_study_k0.03.csv, y en pantalla una tabla de diferencias respecto de la malla más fina.

include("ldf_cycle_runner.jl")
using Printf
using LinearAlgebra

factors = isempty(ARGS) ? [8.0, 4.0, 2.0] : parse.(Float64, ARGS)
reltol = 1e-8; abstol = 1e-10
carpeta = joinpath(@__DIR__, "verificacion_F3")
mkpath(carpeta)

# Repite el estudio para k=3 y 0,03
for k in (3, 0.03)
    rows = DataFrame(dr_factor=Float64[], n_r=Int[], dr_m=Float64[], time_s=Float64[],
        T_c_max=Float64[], T_c_1042=Float64[], T_mid_1042=Float64[], T_w_1042=Float64[], P_1042=Float64[], m_ads_1042_g=Float64[],
        T_c_4131=Float64[], P_4131=Float64[])
    profiles = Dict{Float64,Tuple{Vector{Float64},Vector{Float64}}}()
    for f in factors
        res = run_cycle(13; k=k, dr_factor=f, reltol=reltol, abstol=abstol)
        s = res["summary"]; prm = res["params"]
        ch = s[s.stage.=="charge", :]; di = s[s.stage.=="discharge", :]
        Ms = prm.ρₛ * (1 - prm.ε_b) * prm.V
        push!(rows, (f, prm.n_r, prm.r_span[2] - prm.r_span[1], sum(res["timing"]),
            maximum(ch.T_center), ch.T_center[end], ch.T_mid[end], ch.T_wall[end], ch.P[end], Ms * prm.M_H2 * ch.na_avg_cyl[end] * 1e3,
            di.T_center[end], di.P[end]))
        sol = res["sols"]["charge"]
        profiles[f] = (prm.r_span, sol.u[end][1:prm.n_r])
        @printf("dr factor %.3g: n_r = %d, T_c,max = %.5f K, P(1042) = %.6f MPa, %.1f s\n", f, prm.n_r, maximum(ch.T_center), ch.P[end] / 1e6, sum(res["timing"]))
    end
    CSV.write(joinpath(carpeta, "ldf_grid_study_k$(k).csv"), rows)

    # Differences with respect to the finest grid
    fmin = minimum(factors); ref = rows[rows.dr_factor.==fmin, :][1, :]
    rref, Tref = profiles[fmin]
    println("\nDifferences with respect to the finest grid (k = $(k) 1/s, dr factor $(fmin), n_r = $(ref.n_r)):")
    println("dr_factor  n_r   dT_c,max      dT_c(1042)    dP(1042)/Pa   dm_ads/g      dT_w(1042)    dT_c(4131)    profile RMS (vol.-weighted) K")
    for r in eachrow(rows)
        r.dr_factor == fmin && continue
        rr, TT = profiles[r.dr_factor]
        Ti = [Tref[argmin(abs.(rref .- x))] for x in rr]     # nearest-node reference (grids are nested)
        w = copy(rr); w[1] *= 0.5; w[end] *= 0.5
        rms = sqrt(sum(w .* (TT .- Ti) .^ 2) / sum(w))
        @printf("%-9.3g %5d  %+.3e  %+.3e  %+.3e  %+.3e  %+.3e  %+.3e  %.3e\n", r.dr_factor, r.n_r,
            r.T_c_max - ref.T_c_max, r.T_c_1042 - ref.T_c_1042, r.P_1042 - ref.P_1042, r.m_ads_1042_g - ref.m_ads_1042_g,
            r.T_w_1042 - ref.T_w_1042, r.T_c_4131 - ref.T_c_4131, rms)
    end
end
