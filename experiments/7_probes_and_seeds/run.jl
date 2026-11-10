# Appendix: the MatfreeEK1 for different numbers of probes m and random seeds, on the 1D
# Brusselator (d = 80) at two stiffnesses α and two tolerances, with the EK1 and the EK1 with
# block-diagonal covariance as references.
# The runtime is the minimum over three solves (two for solves over a minute), after one
# warm-up solve per method. Run with one thread so the runtimes are comparable at this small d.
# Usage: julia -t 1 --project experiments/7_probes_and_seeds/run.jl 1|10
include("../../src/setup.jl")

const RTOLS = [1e-3, 1e-5]
const MS = [4, 6, 8, 16, 32, 64, 128, 256, 512]  # default 2(q+1) = 8, minimum q+1 = 4
const SEEDS = 1:10

alpha = parse(Float64, ARGS[1])
prob = add_jvp_vjp(brusselator_1d(; alpha))
u_ref = solve(prob, RadauIIA5(); abstol=1e-12, reltol=1e-12, save_everystep=false).u[end]

runs = Any[("EK1", method("EK1"), 0, 1), ("EK1, block-diagonal", method("EK1"), 0, 1)]
for m in MS, seed in SEEDS
    push!(runs, ("MatfreeEK1", method("MatfreeEK1"; nprobes=m), m, seed))
end

# Step the block-diagonal EK1 by hand to truncate after every step
function solve_run(name, alg, kw)
    if name == "EK1, block-diagonal"
        integ = init(prob, alg; save_everystep=false, kw...)
        while integ.t < prob.tspan[2] && SciMLBase.check_error(integ) == SciMLBase.ReturnCode.Success
            step!(integ)
            truncate_to_blockdiag!(integ.cache.x.Σ, length(prob.u0))
        end
        return SciMLBase.check_error(integ), integ.u, integ.stats
    end
    sol = solve(prob, alg; save_everystep=false, kw...)
    return sol.retcode, sol.u[end], sol.stats
end
for (name, alg, _, _) in unique(r -> r[1], runs)
    solve_run(name, alg, (; abstol=1e-5, reltol=1e-2))
end

path = joinpath(@__DIR__, "data", "alpha$(ARGS[1]).csv")
mkpath(dirname(path))
rm(path; force=true)
for rtol in RTOLS, (name, alg, m, seed) in runs
    Random.seed!(seed)
    kw = (; abstol=rtol * 1e-3, reltol=rtol, maxiters=10^5)
    time = @elapsed retcode, u, stats = solve_run(name, alg, kw)
    ok = SciMLBase.successful_retcode(retcode)
    for _ in 1:(ok ? nrepeats(time, 2) : 0)
        Random.seed!(seed)
        time = min(time, @elapsed solve_run(name, alg, kw))
    end
    row = (; method=name, rtol, m, seed, status=string(retcode),
           error=ok ? mean(abs, u .- u_ref) : NaN, time, stats.naccept, stats.nreject)
    @printf("α=%-4g rtol=%.0e %-19s m=%-3d seed=%-2d %-10s error=%.2e steps=%d\n", alpha, rtol,
            name, m, seed, row.status, row.error, row.naccept)
    flush(stdout)
    CSV.write(path, [row]; append=isfile(path))
end
