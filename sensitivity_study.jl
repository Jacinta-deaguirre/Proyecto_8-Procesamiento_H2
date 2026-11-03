# Parametric sensitivity of the charging stage of Test No. 13 (Supplementary Material 2).
# Each of h, k_w, c_pw, m_w and k_eff is perturbed by -100 % (x1e-6), -50 %, -30 %, +30 %,
# +50 % and +100 % relative to the base case; the centre temperature history is written to
# outputs/sensitivity/<param>.csv and plotted in outputs/sensitivity/<param>.png/.svg.
ENV["GKSwstype"] = "100"
include("cycle_runner.jl")
using Plots
using Printf

mkpath("outputs/sensitivity")
factors = [1e-6, 0.5, 0.7, 1.0, 1.3, 1.5, 2.0]
labels = ["-100 %", "-50 %", "-30 %", "base", "+30 %", "+50 %", "+100 %"]
params = [("h", :U), ("k_w", :k_wall_factor), ("c_pw", :c_wall_factor), ("m_w", :m_wall_factor), ("k_eff", :k_eff_factor)]

summary = DataFrame(parameter=String[], factor=Float64[], T_c_max=Float64[], dT_peak_pct=Float64[])
for (pname, key) in params
    df = nothing
    plt = plot(xlabel="Time / s", ylabel="Centre temperature / K", title="Sensitivity to $(pname)", legend=:bottomright, size=(900, 600))
    Tbase = nothing
    for (f, lab) in zip(factors, labels)
        kw = key == :U ? Dict(:U => 36.0 * f) : Dict(key => f)
        res = run_cycle(13; stages_to_run=("charge",), kw...)
        s = res["summary"]
        if df === nothing; df = DataFrame(t=s.t); end
        df[!, Symbol(lab)] = s.T_center
        plot!(plt, s.t, s.T_center, label=lab, lw=2)
        f == 1.0 && (Tbase = maximum(s.T_center))
        push!(summary, (pname, f, maximum(s.T_center), NaN))
        @printf("%s x%.2g: T_c,max = %.2f K\n", pname, f, maximum(s.T_center))
    end
    base = summary[(summary.parameter.==pname).&(summary.factor.==1.0), :T_c_max][1]
    for r in eachrow(summary)
        r.parameter == pname && (r.dT_peak_pct = (r.T_c_max - base) / (base - 281.0) * 100)
    end
    CSV.write("outputs/sensitivity/$(pname).csv", df)
    savefig(plt, "outputs/sensitivity/$(pname).png"); savefig(plt, "outputs/sensitivity/$(pname).svg")
end
CSV.write("outputs/sensitivity/summary.csv", summary)
println(summary)
