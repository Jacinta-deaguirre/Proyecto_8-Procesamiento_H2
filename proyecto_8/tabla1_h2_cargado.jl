#tabla1_h2_cargado.jl - proyecto8 (Fase 4)
# H2 cargado hasta P límite por tamaño de pellet y tiempo de llenado

using CSV, DataFrames

carpeta = joinpath(@__DIR__, "resultados_F4")
f4 = CSV.read(joinpath(carpeta, "ldf_sensibilidad_carga.csv"), DataFrame)

#H2 cargado (%)
function h2_pct(k, t_c)
    fila = f4[(f4.k_1_s .== k) .& (f4.t_carga_s .== t_c), :]
    return round(fila.H2_cargado_pct[1], digits=1)
end

pellets = [(3.0, "2 mm"), (0.12, "1 cm"), (0.03, "2 cm"), (4.8e-3, "5 cm"), (1.2e-3, "10 cm")]
tabla1 = DataFrame()

for (k, nombre) in pellets
    push!(tabla1, (Pellet = nombre, k_1_s = k, tau_s = round(1/k, digits=1),
        min17 = h2_pct(k, 1042.0), min9 = h2_pct(k, 540.0),
        min5 = h2_pct(k, 300.0), min3 = h2_pct(k, 180.0)))
end

show(tabla1, allrows=true, allcols=true)
println()
CSV.write(joinpath(carpeta, "tabla1_h2_cargado.csv"), tabla1)