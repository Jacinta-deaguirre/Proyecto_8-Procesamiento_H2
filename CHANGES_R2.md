# Changes made for the second revision (September 2026)

All changes respond to Reviewer 3 (second round). The model is now defined in one place,
`src/adsorption_systems.jl::MDA_adsorption!` (DAE residual) and
`src/HydrogenAdsorption.jl::heat_equation` (energy balance); every script calls these
functions. The previous copies of the residual inside each `finite_model_xiao*.jl` script
have been removed.

## Physics
1. **Compression work added to the energy balance** (`heat_equation`). The balance is now
   written in enthalpy form: the gas accumulation term uses `c_p`, and the source
   `Q_compression = ε_b dP/dt` is included (counterpart of `Q_p` in the finite-element model
   of Xiao et al. 2012). `dP/dt` is evaluated from the time derivative of the ideal-gas
   closure (item 3) using the differential variables only.
2. **Outflow enthalpy during discharge.** The inlet term `max(ṁ,0)/V c_p (T_H2 − T)` acts only
   while hydrogen enters the tank; the outflow leaves at the local bed temperature. The
   previous scripts applied `ṁ c_p (T_H2 − T)` with `ṁ < 0` and `T_H2 = 297.6 K`, which was a
   spurious cooling term.
3. **Ideal-gas closure with the harmonic volume-mean temperature.** `P = ρ̄ R_H2 T_H`,
   `T_H = 1/<1/T>`, imposed as an algebraic equation (`differential_vars[2n_r+2] = false`)
   instead of the integrated `dP/dt` ODE. This is exact for uniform pressure and removes the
   drift of the integrated pressure with respect to the equation of state.
4. **Adsorbent mass** `mₛ = ρₛ (1 − ε_b) V = 0.6585 kg` is used consistently in the energy and
   mass balances and in the post-processing (previously 0.671 kg in some places).
5. **Volume averages** (`cyl_avg`, composite trapezoidal rule with weight `r`) are used for
   `n̄_a` in Eqs. (2)–(3) and for every reported average; the post-processing previously used
   the arithmetic node mean.

## Numerics
6. **Stage restarts carry the full state.** Each stage writes `outputs/<stage>_final_state.csv`
   (state `u` and derivative `du`); the next stage starts from it through the new keyword
   arguments of `dae_setup` (`T_init`, `ρ_avg_init`, `T_wall_init`, `du_init`). Previously
   `ρ̄` was re-initialised from the arithmetic mean of the nodal densities, which broke the
   global mass balance by 0.06 g at each restart.
7. **Grid covers the radius exactly.** `coefficient_matrix` now uses
   `n_r = round(R/Δr) + 1`, `r = range(0, R, length = n_r)` and the actual spacing `step(r)`
   in the finite-difference coefficients; the scripts take `dr = r_span[2] − r_span[1]`.
   Previously the nominal `Δr` was used with a grid whose spacing differed by 0.05 %.
8. Isotherm and isosteric-heat evaluations clamp their logarithm arguments so that a trial
   step of the solver with `P < 0` raises no `DomainError` (no effect on valid states).
9. New scripts: `cycle_runner.jl` (full-cycle runner used by the studies),
   `grid_study.jl` (five grids, tolerances 1e-8/1e-10, several outputs),
   `sensitivity_study.jl` (replaces the manually edited `finite_model_xiao_sensibility.jl`),
   `extra_runs.jl`/`extra_runs2.jl` (runs used in the response), `plot_radial_profiles.py`
   (Fig. 5, extracted from `graphs.ipynb`). `grid_independence.jl` was removed.

## Reproducing the manuscript figures
```
julia --project=. finite_model_xiao.jl            # Test 13 charge
julia --project=. finite_model_xiao_dormant.jl    # storage
julia --project=. finite_model_xiao_discharge.jl  # discharge
julia --project=. finite_model_xiao_full_graph.jl # Figs. 2-4
python plot_radial_profiles.py                    # Fig. 5
julia --project=. finite_model_xiao_test19.jl ... _dormant_test19.jl ... _discharge_test19.jl ... _full_graph_test19.jl   # Figs. 6-8
julia --project=. grid_study.jl 4 2 1 0.5 0.25    # Table SM2.1
julia --project=. sensitivity_study.jl            # Figs. SM2.1-SM2.5
```
Set `GKSwstype=100` in the environment for headless plotting.

## Manuscript figures
`plot_manuscript_figures.py` produces Figs. 2-4 and 6-8 (PDF + PNG in `outputs/manuscript_figs/`) from the dense
1-s output of `extra_runs2.jl` (`outputs/extra/t13_base`, `t20_base`), with the +-2.2 K thermocouple uncertainty
drawn on every experimental temperature point and the stage boundaries marked. The Julia `*_full_graph.jl` scripts
still produce the quick-look SVG/PNG versions without error bars.

## Geometry (15 September 2026)
The bed radius is now the physical internal radius of the tank, `R_T = 0.0469 m` (Xiao et al. 2012), and the
cylindrical length is the equivalent length that conserves the measured bed volume, `L = V/(pi R_T^2) = 0.361 m`
(the 0.45 m internal height includes the non-cylindrical end sections). Previously `L = 0.40 m` and the
volume-equivalent radius `R_T = sqrt(V/(pi L)) = 0.04455 m` were used. With the physical radius the thermocouples
Cr (0.02345 m) and Cw (0.0469 m) are exactly at r_i/2 and r_i, and the grid has N_r = 601 nodes for the same
dr = 7.82e-5 m. Effect on the results: < 0.5 K in every temperature metric, < 0.2 % in the pressure AAD.
`cycle_runner.jl` keeps the old geometry available through `build_parameters(test; r_i=nothing, L=0.4)`.
