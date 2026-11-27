# HydrogenAdsorption

This project is used to model the adsorption of hydrogen inside a cylindrical tank filled with activated carbon. In particular, it consists of a reduced-order model, that couples 1-D radial heat transfer, with mass conservation and the ideal gas equation. A DAE approach has been used to solve numerically the problem.

The program outputs temperature, pressure and mass profiles.
## Version 1.1 (September 2026)

The changes made for the second revision of the article (compression work in the energy balance, algebraic ideal-gas closure, physical bed radius, full-state stage restarts, grid and sensitivity studies, figure pipeline) are documented in [CHANGES_R2.md](CHANGES_R2.md).

## Version 2 (Modelo LDF - Noviembre, 2026)
Proyecto 8 - Resistencia a la transferencia de masa (modelo LDF) en HydrogenAdsorption.jl

Proyecto del ramo IIQ3843 - Procesamiento de Hidrógeno para Energías Sostenibles.

El modelo base de HydrogenAdsorption.jl supone equilibrio instántaneo de adsorción (n_a = n_eq(P, T)). En este proyecto se reemplaza por el modelo LDF (∂nₐ/∂t = k(T) * [n_eq(P, T) - n_a]) y se estudia cuánto hidrógeno se alcanza a cargar hasta la presión límite según el tamaño del pellet (k) y el tiempo de llenado. Con k = infinito, da el modelo base.

#Estructura
-src\adsorption_systems.jl -> se encuentra el modelo LDF
-proyecto_8\verificacion_F3 -> Resultados de la verificación (Fase 3)
-proyecto_8\resultados_F4 -> Datos, figuras y tabla del informe (Fase 4)

Cada código se corre en PowerShell con julia --project=. proyecto_8/<script>.jl