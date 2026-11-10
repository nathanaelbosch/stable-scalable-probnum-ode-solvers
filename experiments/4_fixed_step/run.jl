# Fixed-step solves of Burgers' equation and the 1D Fisher-KPP equation, with an IWP(2) prior
# and a fixed diffusion (Figure 4).
# Usage: julia --project experiments/4_fixed_step/run.jl burgers|fisher_kpp [method]
include("../../src/setup.jl")

const PLOTTED = ["EK1", "ExpEK", "DiagonalEK1+IEKF", "MatfreeEK1", "MatfreeEK1+IEKF"]
const PROBLEMS = Dict(
    "burgers" => (; make=(; kw...) -> burgers(; N=200, kw...), dts=10.0 .^ -(0:0.5:3.5)),
    "fisher_kpp" => (; make=(; kw...) -> fisher_kpp_1d(; N=100, kw...), dts=10.0 .^ -(0:0.5:2.5)),
)

key = ARGS[1]
(; make, dts) = PROBLEMS[key]
dir = joinpath(@__DIR__, "data", key)
mkpath(dir)

# The reference solution, also saved at 200 time points for the plot of the solution
prob = make(diag_jac=false)
ts = range(prob.tspan..., length=200)
ref = solve(prob, RadauIIA5(); abstol=1e-20, reltol=1e-20, saveat=ts)
writedlm(joinpath(dir, "solution.txt"), reduce(hcat, ref.u), ',')

prob, prob_diag = add_jvp_vjp(prob), make(diag_jac=true)
for name in (length(ARGS) > 1 ? ARGS[2:2] : PLOTTED)
    p = uses_diag_jac(name) ? prob_diag : prob
    alg = method(name; order=2, diffusionmodel=FixedDiffusion())
    path = joinpath(dir, "$name.csv")
    rm(path; force=true)
    for dt in dts
        sol = solve_once(p, alg; adaptive=false, dt)
        ok = SciMLBase.successful_retcode(sol)
        row = (; method=name, dt, status=string(sol.retcode),
               error=ok ? mean(abs, sol.u[end] .- ref.u[end]) : NaN)
        @printf("%-18s dt=%.1e  %-10s err=%.2e\n", name, dt, row.status, row.error)
        flush(stdout)
        CSV.write(path, [row]; append=isfile(path))
    end
end
