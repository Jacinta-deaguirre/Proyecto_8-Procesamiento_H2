# Extra runs used in the response to reviewers: Test 20 with h = 60 and 100 W/m2K, and the
# tolerance study of Test 13 (rtol 1e-6). Writes outputs/extra/<tag>/summary.csv
include("cycle_runner.jl")
mkpath("outputs/extra")
for (tag, kw) in [("t20_h60", (test=20, U=60.0)), ("t20_h100", (test=20, U=100.0)), ("t13_tol6", (test=13, reltol=1e-6, abstol=1e-9)), ("t13_tol8", (test=13, reltol=1e-8, abstol=1e-10))]
    args = Dict(pairs(kw)); test = pop!(args, :test)
    res = run_cycle(test; args...)
    mkpath("outputs/extra/$tag"); CSV.write("outputs/extra/$tag/summary.csv", res["summary"])
    open("outputs/extra/$tag/meta.txt", "w") do io; println(io, "tag=$tag test=$test n_r=$(res["params"].n_r) dr=$(res["params"].r_span[2]-res["params"].r_span[1]) $(kw)\ntiming_s=", join(round.(res["timing"], digits=2), ",")); end
    println(tag, " done: ", round.(res["timing"], digits=1), " s")
end
