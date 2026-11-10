# Appendix: calibration of the final-time posterior on Lorenz96 (d = 80, T = 1) at fixed
# steps. For each solver we save the error and χ²/d of the final state under its calibrated
# posterior. The MatfreeEK1 runs with several numbers of probes m and 10 seeds each.
# Usage: julia --project experiments/6_calibration/run.jl
include("../../src/setup.jl")

const DTS = [1e-2, 5e-3, 2e-3]
const MS = [4, 6, 8, 16, 32, 64, 128, 256, 512]  # the default is 2(q+1) = 8
const SEEDS = 1:10

# Solve with fixed step `dt` and a global calibrated diffusion. We calibrate by hand. With a
# zero initial covariance and an initial diffusion of 1, the calibrated covariance is the
# global diffusion estimate times the uncalibrated one.
function final_posterior(prob, alg, dt; truncate=false)
    integ = init(prob, alg; dt, adaptive=false, save_everystep=false, dense=false)
    while integ.t < prob.tspan[2] - dt / 2
        t = integ.t
        step!(integ)
        integ.t > t || error("step failed at t = $t")
        truncate && truncate_to_blockdiag!(integ.cache.x.Σ, length(prob.u0))
    end
    pu = ProbNumDiffEq._gaussian_mul!(integ.cache.pu_tmp, integ.cache.SolProj, integ.cache.x)
    return copy(pu.μ), integ.cache.global_diffusion .* Matrix(pu.Σ)
end

prob = lorenz96(; N=80, tspan=(0.0, 1.0))
prob_diag = lorenz96(; N=80, tspan=(0.0, 1.0), diag_jac=true)
prob_mf = add_jvp_vjp(prob)
u_ref = solve(prob, Vern9(); abstol=1e-14, reltol=1e-14).u[end]

diffusion = FixedDiffusion()  # calibrated, global
runs = [("EK0", prob, method("EK0"; diffusionmodel=diffusion), false, 1),
        ("DiagonalEK1", prob_diag, method("DiagonalEK1"; diffusionmodel=diffusion), false, 1),
        ("EK1", prob, method("EK1"; diffusionmodel=diffusion), false, 1),
        ("EK1, block-diagonal", prob, method("EK1"; diffusionmodel=diffusion), true, 1),
        ("DiagonalEK1+IEKF", prob_diag, method("DiagonalEK1+IEKF"; diffusionmodel=diffusion), false, 1)]
for seed in SEEDS
    push!(runs, ("MatfreeEK1+IEKF, m=8", prob_mf, method("MatfreeEK1+IEKF"; diffusionmodel=diffusion),
                 false, seed))
end
for m in MS, seed in SEEDS
    push!(runs, ("MatfreeEK1, m=$m", prob_mf,
                 method("MatfreeEK1"; nprobes=m, diffusionmodel=diffusion), false, seed))
end

path = joinpath(@__DIR__, "data", "lorenz96.csv")
mkpath(dirname(path))
rm(path; force=true)
for dt in DTS, (name, p, alg, truncate, seed) in runs
    Random.seed!(seed)
    μ, Σ = final_posterior(p, alg, dt; truncate)
    diff = μ .- u_ref
    row = (; method=name, dt, seed, status="Success", error=mean(abs, diff),
           chi2_d=dot(diff, cholesky(Symmetric(Σ)) \ diff) / length(μ))
    @printf("dt=%.0e %-22s seed=%-2d error=%.1e chi2/d=%.3g\n", dt, name, seed, row.error,
            row.chi2_d)
    flush(stdout)
    CSV.write(path, [row]; append=isfile(path))
end
