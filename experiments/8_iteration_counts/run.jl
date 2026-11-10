# Appendix: the Krylov and IEKF iteration counts on the sweeps of Figure 1. We solve once per
# successful point of Figure 1, with the same settings, so run 1_fhn first. The counts over the
# tolerance come from the work-precision runs of Figure 5.
# Usage: julia --project experiments/8_iteration_counts/run.jl domain|stiffness|refinement [method]
include("../../src/setup.jl")

const METHODS = ["MatfreeEK1", "MatfreeEK1+IEKF", "DiagonalEK1+IEKF"]

sweep = ARGS[1]
dir = joinpath(@__DIR__, "data", sweep)
mkpath(dir)
for name in (length(ARGS) > 1 ? ARGS[2:2] : METHODS)
    path = joinpath(dir, "$name.csv")
    rm(path; force=true)
    fig1 = CSV.File(joinpath(@__DIR__, "..", "1_fhn", "data", sweep, "$name.csv"))
    for r in fig1
        r.status == "Success" || continue
        (; N, L) = r
        prob = fitzhugh_nagumo_2d(; N, L, diag_jac=uses_diag_jac(name))
        prob = uses_diag_jac(name) ? prob : add_jvp_vjp(prob)
        sol = solve_once(prob, method(name); abstol=1e-3, reltol=1e-1)
        s, pn = sol.stats, sol.pnstats
        row = (; method=name, N, L, d=2N^2, dx=L / N, status=string(sol.retcode), s.naccept,
               s.nnonliniter, s.nnonlinconvfail, pn.nsolve_mean, pn.niter_mean,
               pn.nsolve_covariance, pn.niter_covariance, pn.nsolve_calibration,
               pn.niter_calibration)
        @printf("%-18s N=%-4d L=%-8.4g %-10s steps=%d\n", name, N, L, row.status, s.naccept)
        flush(stdout)
        CSV.write(path, [row]; append=isfile(path))
        GC.gc()
    end
end
