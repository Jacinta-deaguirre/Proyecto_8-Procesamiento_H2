# Dense (1 s) output of the base cases with the manuscript tolerances, used for the error metrics.
include("cycle_runner.jl")
for (tag, test) in [("t13_base", 13), ("t20_base", 20)]
    res = run_cycle(test)
    mkpath("outputs/extra/$tag"); CSV.write("outputs/extra/$tag/summary.csv", res["summary"])
    open("outputs/extra/$tag/meta.txt", "w") do io; println(io, "tag=$tag test=$test n_r=$(res["params"].n_r) dr=$(res["params"].r_span[2]-res["params"].r_span[1]) base\ntiming_s=", join(round.(res["timing"], digits=2), ",")); end
    println(tag, " done: ", round.(res["timing"], digits=1), " s")
end
