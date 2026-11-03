# Grid-independence study (Test No. 13, complete cycle) with tight time-integration
# tolerances so that the differences between grids measure the spatial error only.
# Outputs: outputs/grid_study.csv (one row per grid) and a printed table of differences
# with respect to the finest grid for several outputs.
#
#   julia --project=. grid_study.jl            # factors 4, 2, 1, 0.5 (0.25 optional, slow)
#   julia --project=. grid_study.jl 4 2 1 0.5 0.25
include("cycle_runner.jl")
using Printf
using LinearAlgebra

factors = isempty(ARGS) ? [4.0, 2.0, 1.0, 0.5] : parse.(Float64, ARGS)
reltol = 1e-8; abstol = 1e-10
mkpath("outputs")

rows = DataFrame(dr_factor=Float64[], n_r=Int[], dr_m=Float64[], time_s=Float64[],
    T_c_max=Float64[], T_c_1042=Float64[], T_mid_1042=Float64[], T_w_1042=Float64[], P_1042=Float64[], m_ads_1042_g=Float64[],
    T_c_4131=Float64[], P_4131=Float64[])
profiles = Dict{Float64,Tuple{Vector{Float64},Vector{Float64}}}()
for f in factors
    res = run_cycle(13; dr_factor=f, reltol=reltol, abstol=abstol)
    s = res["summary"]; prm = res["params"]
    ch = s[s.stage.=="charge", :]; di = s[s.stage.=="discharge", :]
    Ms = prm.ρₛ * (1 - prm.ε_b) * prm.V
    push!(rows, (f, prm.n_r, prm.r_span[2] - prm.r_span[1], sum(res["timing"]),
        maximum(ch.T_center), ch.T_center[end], ch.T_mid[end], ch.T_wall[end], ch.P[end], Ms * prm.M_H2 * ch.na_avg_cyl[end] * 1e3,
        di.T_center[end], di.P[end]))
    sol = res["sols"]["charge"]
    profiles[f] = (prm.r_span, sol.u[end][1:prm.n_r])
    @printf("dr factor %.3g: n_r = %d, T_c,max = %.5f K, P(1042) = %.6f MPa, %.1f s\n", f, prm.n_r, maximum(ch.T_center), ch.P[end] / 1e6, sum(res["timing"]))
end
CSV.write("outputs/grid_study.csv", rows)

# Differences with respect to the finest grid
fmin = minimum(factors); ref = rows[rows.dr_factor.==fmin, :][1, :]
rref, Tref = profiles[fmin]
println("\nDifferences with respect to the finest grid (dr factor $(fmin), n_r = $(ref.n_r)):")
println("dr_factor  n_r   dT_c,max      dT_c(1042)    dP(1042)/Pa   dm_ads/g      dT_w(1042)    dT_c(4131)    profile RMS (vol.-weighted) K")
for r in eachrow(rows)
    r.dr_factor == fmin && continue
    rr, TT = profiles[r.dr_factor]
    Ti = [Tref[argmin(abs.(rref .- x))] for x in rr]     # nearest-node reference (grids are nested)
    w = copy(rr); w[1] *= 0.5; w[end] *= 0.5
    rms = sqrt(sum(w .* (TT .- Ti) .^ 2) / sum(w))
    @printf("%-9.3g %5d  %+.3e  %+.3e  %+.3e  %+.3e  %+.3e  %+.3e  %.3e\n", r.dr_factor, r.n_r,
        r.T_c_max - ref.T_c_max, r.T_c_1042 - ref.T_c_1042, r.P_1042 - ref.P_1042, r.m_ads_1042_g - ref.m_ads_1042_g,
        r.T_w_1042 - ref.T_w_1042, r.T_c_4131 - ref.T_c_4131, rms)
end
