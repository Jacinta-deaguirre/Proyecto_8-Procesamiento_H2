# Compare the volume-equivalent geometry (r_i = 0.04455 m, L = 0.40 m) with the physical radius (r_i = 0.0469 m, L = V/(pi r_i^2))
include("cycle_runner.jl")
mkpath("outputs/extra")
for (tag, kw) in [("t13_phys", (test=13, r_i=0.0469)), ("t20_phys", (test=20, r_i=0.0469))]
    args = Dict(pairs(kw)); test = pop!(args, :test)
    res = run_cycle(test; args...)
    mkpath("outputs/extra/$tag"); CSV.write("outputs/extra/$tag/summary.csv", res["summary"])
    prm = res["params"]
    open("outputs/extra/$tag/meta.txt", "w") do io; println(io, "tag=$tag test=$test n_r=$(prm.n_r) dr=$(prm.r_span[2]-prm.r_span[1]) R_T=$(prm.R_T) L=$(prm.geometric.L) Ai=$(prm.geometric.Ai) Ao=$(prm.geometric.Ao)\ntiming_s=", join(round.(res["timing"], digits=2), ",")); end
    println(tag, " done: R_T=", prm.R_T, " L=", prm.geometric.L, " n_r=", prm.n_r)
end
