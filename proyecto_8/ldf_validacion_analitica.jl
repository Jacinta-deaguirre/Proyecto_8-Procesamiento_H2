#ldf_validacion_analitica.jl-Proyrcto 8
# Validacion del modelo LDF contra solucion analitica (T, P y k constantes):
# estanque presurizado a 9MPa, sin flujo, T inicial (281K),
# poco carbón (p_s 1e-8) e inicialmente vacío
ENV["GKSwstype"] = "100"
include("ldf_cycle_runner.jl")
using Plots

#Parametros del Test 13
prm = build_parameters(13; dr_factor=8.0, T_air=281.0)
n_r = prm.n_r
T₀ = 281.0
P₀ = 9.0e6

#Poco carbón: lo que se absorbe no cambia P ni libera calor apreciable
f_ads = 1e-8
mat = prm.material
mat_chico = MaterialProperties(mat.ρₛ * f_ads, mat.cₛ, mat.mₛ * f_ads, mat.kₛ, mat.ε_b, mat.cₚ, mat.cᵥ,
    mat.M_H2, mat.R, mat.k_g, mat.k_eff)

n_eq = adsorption_isotherm(prm.MDA_params, P₀, T₀)

function correr_relajacion(k)
    op = OperationalParameters(prm.U, T₀, 0.0, T₀, k)
    par = AdsorptionParameters(prm.MDA_params, mat_chico, prm.geometric, op)
    u₀, du₀, dv = dae_setup(prm.MDA_params, mat_chico, prm.geometric, op, P₀, T₀; 
        T_init = fill(T₀, n_r), nₐ_init = zeros(n_r))
    prob = DAEProblem(MDA_adsorption!, du₀, u₀, (0.0, 10 / k), par, differential_vars=dv)
    return solve(prob, IDA(linear_solver=:LapackDense), reltol=1e-8, abstol=1e-10)
end

# Comparacion con la solucion analitica 
k_valores = [3.0, 0.12, 0.03, 4.8e-3, 1.2e-3]
marcadores = [:circle, :square, :diamond, :utriangle, :star5]
tabla = DataFrame()

x = range(0, 10, length=200)
fig = plot(x, 1 .- exp.(-x), label="Solucion analitica: 1 - exp(-k t)", xlabel="k*t (adimensional)", ylabel="n_a / n_eq (promedio en el lecho)", legend=:bottomright, size=(800,550))

for (i,k) in enumerate(k_valores)
    sol_k = correr_relajacion(k)
    n_mod = [cyl_avg(u[n_r+1:2*n_r], prm.r_span, prm.R_T) for u in sol_k.u]
    n_ana = n_eq .* (1 .- exp.(-k.* sol_k.t))
    T_c = [u[1] for u in sol_k.u]
    P = [u[2*n_r+2] for u in sol_k.u]
    push!(tabla, (k_1_s = k,
        error_max_pct = 100* maximum(abs.(n_mod .- n_ana)) / n_eq,
        cambio_max_T_K =maximum(abs.(T_c .- T₀)),
        cambio_max_P_pct= 100* maximum(abs.(P .- P₀))/P₀))
    x_m = (0:0.5:9.5) .+ .1 *(i-1)
    n_m = [cyl_avg(sol_k(xx / k)[n_r+1:2*n_r], prm.r_span, prm.R_T) for xx in x_m]
    scatter!(fig, x_m, n_m ./n_eq, label="Modelo, k = $(k) 1/s", marker=marcadores[i], markersize=5)
end

show(tabla, allrows=true, allcols=true)
println()
mkpath(joinpath(@__DIR__,"outputs"))
CSV.write(joinpath(@__DIR__, "outputs", "ldf_validacion_analitica.csv"), tabla)
savefig(fig, joinpath(@__DIR__, "outputs", "ldf_validacion_analitica.svg"))
display(fig)