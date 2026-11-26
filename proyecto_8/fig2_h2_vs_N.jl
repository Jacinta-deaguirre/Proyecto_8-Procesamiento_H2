#fig2_h2_vs_n.jl - proyecto_8 (Fase 4)
#H2 cargado hasta P límite v/s N 

ENV["GKSwstype"] = "100"
using CSV, DataFrames, Plots

carpeta = joinpath(@__DIR__, "resultados_F4")
f4 = CSV.read(joinpath(carpeta, "ldf_sensibilidad_carga.csv"), DataFrame)

fig = plot(xscale=:log10, xlabel="N = k*t_llenado", ylabel="H2 cargado hasta P límite (%)",
    legend=:bottomright, size=(800,550), ylims=(45,101), dpi=300, framestyle=:box,
    xticks=([0.1,1,10,100,1000], ["0,1", "1", "10", "100", "1000"]))
hline!(fig, [99,95], color=:lightgray, linestyle=:dash, linewidth=1, label="99% y 95%")
vline!(fig, [50], color=:gray, linestyle=:dot, linewidth=1.5, label="N = 50")

#Curva de la f4
f4_orden = sort(f4, :N)
plot!(fig, f4_orden.N, f4_orden.H2_cargado_pct, color=:gray, linewidth=1, label = "")

#Puntos
llenados =[(1042.0, "17 min (Test 13)", :black, :circle),
            (540.0, "9 min (CR1)", :blue, :diamond),
            (300.0, "5 min", :green, :square),
            (180.0, "3 min", :red, :utriangle)]

for (t_c, etiqueta, color, forma) in llenados
    filas = f4[f4.t_carga_s .== t_c, :]
    scatter!(fig, filas.N, filas.H2_cargado_pct, color=color, marker=forma, markersize=6,
    markerstrokecolor=:white, markerstrokewidth=1, label="Llenado de "*etiqueta)
end

savefig(fig, joinpath(carpeta, "fig2_h2_vs_N.svg"))
savefig(fig, joinpath(carpeta, "fig2_h2_vs_N.png"))
display(fig)