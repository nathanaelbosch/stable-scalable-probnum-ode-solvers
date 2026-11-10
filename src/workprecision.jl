# Repeat a solve `numruns` times for timing. Solves over a minute vary little, so they repeat
# only once. That is enough to exclude compilation from the minimum.
nrepeats(time, numruns) = time < 60 ? numruns : min(numruns, 1)

function solve_once(prob, alg; kwargs...)
    Random.seed!(1)  # the probes of the `MatfreeEK1` are random
    return solve(prob, alg; save_everystep=false, maxiters=10^7, kwargs...)
end

const ABSTOLS = 10.0 .^ -(2:12)
const RELTOLS = 10.0 .^ -(-1:9)

# Writes one row per method and tolerance to `dir/<method>.csv`, with the mean absolute error
# at the final time. Set `RTOLS` (e.g. "10,1e-3") to run only some of the tolerances.
function run_workprecision(dir, methods, make_problem, ref_alg, ref_tol; rtols=RELTOLS,
                           numruns=3)
    haskey(ENV, "RTOLS") && (rtols = intersect(rtols, parse.(Float64, split(ENV["RTOLS"], ","))))
    prob = make_problem(diag_jac=false)
    u_ref = solve(prob, ref_alg; abstol=ref_tol, reltol=ref_tol, save_everystep=false).u[end]
    prob = add_jvp_vjp(prob)
    prob_diag = make_problem(diag_jac=true)
    mkpath(dir)
    for name in methods
        path = joinpath(dir, "$name.csv")
        rm(path; force=true)
        p, alg = uses_diag_jac(name) ? prob_diag : prob, method(name)
        for (atol, rtol) in zip(ABSTOLS, RELTOLS)
            rtol in rtols || continue
            time = @elapsed sol = solve_once(p, alg; abstol=atol, reltol=rtol)
            ok = SciMLBase.successful_retcode(sol)
            for _ in 1:(ok ? nrepeats(time, numruns) : 0)
                time = min(time, @elapsed solve_once(p, alg; abstol=atol, reltol=rtol))
            end
            s = sol.stats
            # Krylov counts of the matrix-free solvers, zero for the others
            pn = sol isa ProbNumDiffEq.ProbODESolution ? sol.pnstats :
                 (; nsolve_mean=0, niter_mean=0, nsolve_covariance=0, niter_covariance=0,
                  nsolve_calibration=0, niter_calibration=0)
            row = (; method=name, rtol, status=string(sol.retcode),
                   error=ok ? mean(abs, sol.u[end] .- u_ref) : NaN, time, s.naccept, s.nreject,
                   s.nnonliniter, s.nnonlinconvfail, pn.nsolve_mean, pn.niter_mean,
                   pn.nsolve_covariance, pn.niter_covariance, pn.nsolve_calibration,
                   pn.niter_calibration)
            @printf("%-18s rtol=%.0e  %-10s err=%.2e  time=%.3fs  steps=%d\n", name, rtol,
                    row.status, row.error, time, s.naccept)
            flush(stdout)
            CSV.write(path, [row]; append=isfile(path))
        end
    end
end
