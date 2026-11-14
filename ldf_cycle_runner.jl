# Copia de cycle_runner.jl con el modelo LDF (Proyecto 8)
#Full charge -> storage -> discharge cycle runner built on the package residual
# `MDA_adsorption!`. Used by grid_study.jl and sensitivity_study.jl; the stage scripts
# finite_model_xiao*.jl solve exactly the same DAE and additionally write the figures.
#
#   include("cycle_runner.jl")
#   res = run_cycle(13; dr_factor=1.0, reltol=1e-3, abstol=1e-6)
#
# `res` is a Dict with the per-stage solutions and a summary DataFrame (1 s output).
using HydrogenAdsorption
using Sundials
using Statistics
using CSV
using DataFrames

function build_parameters(test::Int; dr_factor=1.0, U=36.0, k_eff_factor=1.0, k_wall_factor=1.0,
    c_wall_factor=1.0, m_wall_factor=1.0, T_air=nothing, r_i=0.0469, L=0.4)
    # Tank parameters (Table 1 of the manuscript). Two geometries are possible for the same bed volume V:
    #   r_i = 0.0469 (default): physical bed radius [Xiao et al. 2012], with the equivalent cylindrical length L = V/(pi r_i^2)
    #   r_i = nothing          : volume-equivalent radius R_T = sqrt(V/(pi L)) for the given cylindrical length L (previous versions)
    V = 2.4946e-3; e = 0.0508 - 0.0469
    if r_i !== nothing; L = V / (pi * r_i^2); end
    c_wall = 468 * c_wall_factor; m_wall = 3.714 * m_wall_factor; k_wall = 13 * k_wall_factor
    R_T = r_i === nothing ? sqrt(V / (pi * L)) : r_i
    dr_nominal = (1 / 2)^5 * 0.0025 * dr_factor
    Ai = 2 * pi * R_T * L; Ao = 2 * pi * (R_T + e) * L
    # Isotherm
    MDA_params = MDAParameters(71.6, 1470e6, 3080.0, 18.9, 2.0)
    # Material properties
    ρₛ = 517.6; cₛ = 825; kₛ = 0.646; ε_b = 0.49
    mₛ = ρₛ * (1 - ε_b) * V
    cₚ = 14700.0; cᵥ = 10134.0; M_H2 = 2.0159e-3; R = 8.314; k_g = 0.206
    k_eff = (kₛ * (1 - ε_b) + k_g * ε_b) * k_eff_factor
    T_bath = T_air === nothing ? (test == 13 ? 282.5 : 302.5) : T_air
    n_r, r_span, A, b = coefficient_matrix(R_T, dr_nominal, k_eff, U, T_bath)
    r_span = collect(r_span); dr = r_span[2] - r_span[1]
    material = MaterialProperties(ρₛ, cₛ, mₛ, kₛ, ε_b, cₚ, cᵥ, M_H2, R, k_g, k_eff)
    geometric = GeometricParameters(n_r, dr, V, L, A, b, r_span, R_T, e, Ao, Ai, c_wall, m_wall, k_wall)
    if test == 13
        T₀ = 281.0; Pᵢ = 0.033e6
        stages = [("charge", 0.0, 1042.0, 2.023e-5, 297.6), ("dormant", 1042.0, 3189.0, 0.0, 297.6), ("discharge", 3189.0, 4131.0, -2.184e-5, 297.6)]
    elseif test == 20
        T₀ = 302.4; Pᵢ = 0.032e6
        stages = [("charge", 0.0, 953.0, 2.048e-5, 276.5), ("dormant", 953.0, 3822.0, 0.0, 297.6), ("discharge", 3822.0, 4694.0, -2.186e-5, 277.7)]
    else
        error("test must be 13 or 20")
    end
    return (; MDA_params, material, geometric, U, T_bath, T₀, Pᵢ, stages, n_r, r_span, R_T, V, ε_b, ρₛ, M_H2, R)
end

function run_cycle(test::Int; k=Inf, reltol=1e-3, abstol=1e-6, saveat=1.0, stages_to_run=("charge", "dormant", "discharge"), kwargs...)
    prm = build_parameters(test; kwargs...)
    (; MDA_params, material, geometric, U, T_bath, T₀, Pᵢ, n_r, r_span, R_T, V, ε_b, ρₛ, M_H2, R) = prm
    summary = DataFrame(stage=String[], t=Float64[], P=Float64[], rho_avg=Float64[], T_wall=Float64[],
        T_center=Float64[], T_mid=Float64[], T_edge=Float64[], T_avg_cyl=Float64[], T_avg_node=Float64[], T_harm_cyl=Float64[],
        na_avg_cyl=Float64[], na_avg_node=Float64[], dH_avg_cyl=Float64[])
    sols = Dict{String,Any}()
    u_prev = nothing; du_prev = nothing; timing = Float64[]
    for (sname, t0, tf, m_in, T_H2) in prm.stages
        sname in stages_to_run || continue
        op = OperationalParameters(U, T_bath, m_in, T_H2,k)
        par = AdsorptionParameters(MDA_params, material, geometric, op)
        if u_prev === nothing
            u₀, du₀, dv = dae_setup(MDA_params, material, geometric, op, Pᵢ, T₀)
        else
            u₀, du₀, dv = dae_setup(MDA_params, material, geometric, op, u_prev[2*n_r+2], T₀;
                T_init=u_prev[1:n_r], ρ_avg_init=u_prev[2*n_r+1], T_wall_init=u_prev[3*n_r+3], du_init=du_prev,
                nₐ_init = isinf(k) ? nothing : u_prev[n_r+1:2*n_r])
        end
        prob = DAEProblem(MDA_adsorption!, du₀, u₀, (t0, tf), par, differential_vars=dv)
        tt = @elapsed sol = solve(prob, IDA(linear_solver=:LapackDense), reltol=reltol, abstol=abstol)
        push!(timing, tt)
        sols[sname] = sol
        for ts in t0:saveat:tf   # dense output every `saveat` seconds
            us = sol(ts)
            T = us[1:n_r]; na = us[n_r+1:2*n_r]; P = us[2*n_r+2]
            dH = isinf(k) ? isosteric_heat_of_adsorption(MDA_params, P, T) : isosteric_heat_of_adsorption(MDA_params, na)
            push!(summary, (sname, ts, P, us[2*n_r+1], us[3*n_r+3],
                T[1], T[div(n_r, 2)], T[end], cyl_avg(T, r_span, R_T), mean(T), 1 / cyl_avg(1 ./ T, r_span, R_T),
                cyl_avg(na, r_span, R_T), mean(na), cyl_avg(dH, r_span, R_T)))
        end
        u_prev = sol.u[end]; du_prev = sol.du[end]
    end
    return Dict("summary" => summary, "sols" => sols, "params" => prm, "timing" => timing)
end
