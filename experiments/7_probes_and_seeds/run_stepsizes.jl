# Appendix: the step sizes chosen by the adaptive MatfreeEK1 over time, for 10 seeds at the
# default m, on the 1D Brusselator at two stiffnesses α and two tolerances. Not plotted.
# The text only states how much the step sizes differ between seeds.
# Usage: julia --project experiments/7_probes_and_seeds/run_stepsizes.jl
include("../../src/setup.jl")

const ALPHAS = [1.0, 10.0]
const RTOLS = [1e-3, 1e-5]
const SEEDS = 1:10

rows = []
for alpha in ALPHAS, rtol in RTOLS
    prob = add_jvp_vjp(brusselator_1d(; alpha))
    for seed in SEEDS
        Random.seed!(seed)
        sol = solve(prob, method("MatfreeEK1"); abstol=rtol * 1e-3, reltol=rtol, dense=false)
        ts = sol.t
        append!(rows, [(; method="MatfreeEK1", alpha, rtol, seed, t=ts[i], dt=ts[i+1] - ts[i])
                       for i in 1:length(ts)-1])
        @printf("α=%-4g rtol=%.0e seed=%-2d %s, %d steps\n", alpha, rtol, seed, sol.retcode,
                length(ts) - 1)
        flush(stdout)
    end
end
mkpath(joinpath(@__DIR__, "data"))
CSV.write(joinpath(@__DIR__, "data", "stepsizes.csv"), rows)
