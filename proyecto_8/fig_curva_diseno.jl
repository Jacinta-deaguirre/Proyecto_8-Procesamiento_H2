# fig_curva_diseno.jl - Proyecto 8 (Fase 4)
#Figura de diametro maximo del pellet v/s tiempo de llenado (para N = 50 y N=100)

ENV["GKSwstype"] = "100"
using CSV, DataFrames, Plots

diseno = CSV.read(joinpath(joinpath(@__DIR__, "resultados_F4"), "ldf_curva_diseno.csv"), DataFrame)
D_eff = 2e-7 # difusividad efectiva / m^2/s

fig = plot(xlabel="Tiempo de llenado (min)", ylabel="Diámetro máximo del pellet (mm)",
    legend=:topleft, size=(800,550), xlims=(0, 18), ylims=(0, 20), dpi=300, framestyle=:box)
vspan!(fig, [3,5], color=:lightgray, alpha=0.5, label="Meta DOE de llenado (3 a 5 min)")
hline!(fig, [2.0], color=:gray, linestyle=:dash, linewidth=1, label="Pellet del paper base (2 mm)")

t_min = range(0, 18, length=200)
curvas = [(50.0, "N = 50 (99,0% del H2)", :blue),
            (100.0, "N = 100 (99,5% del H2)", :red)]
for (N, etiqueta, color) in curvas
    d_mm = 2 .* sqrt.(15 .* D_eff .* (t_min .* 60) ./ N) .* 1000
    plot!(fig, t_min, d_mm, color=color, linewidth=2, label=etiqueta)
end

scatter!(fig, diseno.t_carga_s ./ 60, diseno.d_max_mm, color=:black, marker=:circle, markersize=5,
    markerstrokecolor=:white, markerstrokewidth=1.5, label="Corridas del modelo")

annotate!(fig, 9.5, 4.5, text("Pellets bajo la curva: el equilibrio sirve", 8, :left, :black))

savefig(fig, joinpath(joinpath(@__DIR__, "resultados_F4"), "ldf_curva_diseno.svg"))
savefig(fig, joinpath(joinpath(@__DIR__, "resultados_F4"), "ldf_curva_diseno.png"))
display(fig)
