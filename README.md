# HydrogenAdsorption

This project is used to model the adsorption of hydrogen inside a cylindrical tank filled with activated carbon. In particular, it consists of a reduced-order model, that couples 1-D radial heat transfer, with mass conservation and the ideal gas equation. A DAE approach has been used to solve numerically the problem.

The program outputs temperature, pressure and mass profiles.
## Version 1.1 (September 2026)

The changes made for the second revision of the article (compression work in the energy balance, algebraic ideal-gas closure, physical bed radius, full-state stage restarts, grid and sensitivity studies, figure pipeline) are documented in [CHANGES_R2.md](CHANGES_R2.md).

## Version 2 (Modelo LDF - Noviembre, 2026)
Proyecto 8 - Resistencia a la transferencia de masa (modelo LDF) en HydrogenAdsorption.jl

Proyecto del ramo IIQ3843 - Procesamiento de Hidrógeno para Energías Sostenibles.

El modelo base de HydrogenAdsorption.jl supone equilibrio instantáneo de adsorción (n_a = n_eq(P, T)). En este proyecto se reemplaza por el modelo LDF (∂nₐ/∂t = k(T) * [n_eq(P, T) - n_a]) y se estudia cuánto hidrógeno se alcanza a cargar hasta la presión límite según el tamaño del pellet (k) y el tiempo de llenado. Con k = infinito, da el modelo base.

Estructura
- src\adsorption_systems.jl -> se encuentra el modelo LDF
- proyecto_8\verificacion_F3 -> Resultados de la verificación (Fase 3)
- proyecto_8\resultados_F4 -> Datos, figuras y tabla del informe (Fase 4)

### Para generar las figuras del informe
Cada script se corre en PowerShell con: julia --project=. proyecto_8/nombre.jl
Primero corre el script de datos (que guarda un csv en proyecto_8/resultados_F4) y despues el de la figura (lee ese csv). De todas formas, los csv ya se generaron y están en el repositorio así que las figuras se pueden hacer directamente.

Resultados:
- Fig 1: script figura = fig1_carga.jl y script datos = ldf_fig1_carga.jl
- Fig 2: script figura = fig2_h2_vs_N.jl y script datos = ldf_sensibilidad_carga.jl
- Tabla 1: script tabla = tabla1_h2_cargado.jl y script datos = ldf_sensibilidad_carga.jl
- Fig 3: script figura = fig_curva_diseno.jl y script datos = ldf_curva_diseno.jl
- Fig 4: script figura = fig_ldf_kT.jl y script datos = ldf_sensibilidad_kT.jl y ldf_sensibilidad_carga.jl

### Verificación (Fase 3)
- ldf_verificacion.jl (k = Ing y balance de masa)
- ldf_convergencia.jl (paso de tiemp)
- ldf_grid_study-jl (malla)
- ldf_validacion_analitica.jl (solución analítica)
Guardan los resultados en proyecto_8/verificacion_F3