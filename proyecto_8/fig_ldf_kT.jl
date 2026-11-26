#fig_ldf_kT.jl - Proyecto 8, Fase 4
# Figura de H2 cargado v/s N con curva de k constante y puntos con k(T)

ENV["GKSwstype"] = "100"
using CSV, DataFrames, Plots

f4 = CSV.read(joinpath(joinpath(@__DIR__, "resultados_F4"), "ldf_sensibilidad_carga.csv"), DataFrame)
kT = CSV.read(joinpath(joinpath(@__DIR__, "resultados_F4"), "ldf_sensibilidad_kT.csv"), DataFrame)

fig = plot(xscale=:log10, xlabel="N = k_ref*t_llenado", ylabel="H2 cargado hasta P límite (%)",
    legend=:bottomright, size=(800,550), ylims=(45,101), dpi=300, framestyle=:box,
    xticks=([0.1,1,10,100,1000], ["0,1", "1", "10", "100", "1000"]))
hline!(fig, [99,95], color=:lightgray, linestyle=:dash, linewidth=1, label="99% y 95%")

#Curva de la f4
f4_orden = sort(f4, :N)
plot!(fig, f4_orden.N, f4_orden.H2_cargado_pct, color=:gray, linewidth=2,
    label = "Curva de la F4 (k constante, 20 corridas)")

#Puntos con k(T)
estilos =[(10.0, "E_a = 10 kJ/mol", :red, :diamond, 7),
            (1.2, "E_a = 1,2 kJ/mol", :blue, :square, 7),
            (0.0, "E_a = 0 (k constante)", :black, :circle, 4)]
for (Ea, etiqueta, color, forma, tamano) in estilos
    filas = kT[kT.E_a_kJ_mol .== Ea, :]
    scatter!(fig, filas.N, filas.H2_cargado_pct, color=color, marker=forma, markersize=tamano,
    markerstrokecolor=:white, markerstrokewidth=1.5, label=etiqueta)
end

#etiqueta de cada caso
pellet = Dict(0.03 => "2 cm", 0.0048 => "5 cm")
llenado = Dict( 540.0 => "9 min", 180.0 => "3 min")
for fila in eachrow(kT[kT.E_a_kJ_mol .== 0.0, :])
    annotate!(fig, fila.N*1.08, fila.H2_cargado_pct - 3.5,
        text(pellet[fila.k_1_s]*" · " * llenado[fila.t_carga_s], 8, :left, :black))
end

savefig(fig, joinpath(joinpath(@__DIR__, "resultados_F4"), "fig4_sensibilidad_kT.svg"))
savefig(fig, joinpath(joinpath(@__DIR__, "resultados_F4"), "fig4_sensibilidad_kT.png"))
display(fig)