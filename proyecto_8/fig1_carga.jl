#fig1_carga.jl - Proyecto 8 (F4)
#Carga con equilibrio y 5 tamaños de pellet (presión, T del centro, masa adsorbida vs tiempo)
# Lee ldf_fig1_carga.csv 

ENV["GKSwstype"]="100"
using CSV, DataFrames, Plots

carpeta = joinpath(@__DIR__, "resultados_F4")
datos = CSV.read(joinpath(carpeta, "ldf_fig1_carga.csv"), DataFrame)

estilos=[("2 mm", :blue, :solid), ("1 cm", :green, :solid), ("2 cm", :orange, :solid),
        ("5 cm", :red, :solid), ("10 cm", :purple, :solid), ("Equilibrio", :black, :dash)]
filas_fig = [(1042.0, "Llenado de 17 min (Test 13)"), (540.0, "Llenado de 9 min (CR1)")]

paneles=[]
for (j, (t_c, titulo)) in enumerate(filas_fig)
    pP = plot(xlabel="Tiempo (min)", ylabel="Presión (MPa)", title=titulo, titlefontsize=10,
        ylims= (0, 17), legend = j == 1 ? :topleft : false)
    pT = plot(xlabel="Tiempo (min)", ylabel="T en el centro (K)", ylims= (280, 320), legend = false)
    pM = plot(xlabel="Tiempo (min)", ylabel="Masa adsorbida (g)", legend = false)
    for (caso, color, linea) in estilos
        d= datos[(datos.t_carga_s .== t_c) .& (datos.caso .== caso), :]
        plot!(pP, d.t_s ./ 60, d.P_MPa, color=color, linestyle=linea, linewidth=2, label=caso)
        plot!(pT, d.t_s ./ 60, d.T_centro_k, color=color, linestyle=linea, linewidth=2, label=caso)
        plot!(pM, d.t_s ./ 60, d.m_ads_g, color=color, linestyle=linea, linewidth=2, label=caso)
    end
    push!(paneles, pP, pT, pM)
end

fig = plot(paneles..., layout=(2, 3), size=(1500, 900), margin=5Plots.mm, dpi=300, framestyle=:box)
savefig(fig, joinpath(carpeta, "fig1_carga.svg"))
savefig(fig, joinpath(carpeta, "fig1_carga.png"))
display(fig)