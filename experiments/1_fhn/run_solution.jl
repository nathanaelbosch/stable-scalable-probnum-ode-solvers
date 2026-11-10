# The solution of the 2D FitzHugh-Nagumo model at t = 20 on a 250 × 250 grid (d = 125,000),
# computed with the MatfreeEK1 (Figure 1a).
# Usage: julia --project experiments/1_fhn/run_solution.jl
include("../../src/setup.jl")

N = 250
prob = add_jvp_vjp(fitzhugh_nagumo_2d(; N))
time = @elapsed sol = solve_once(prob, method("MatfreeEK1"); abstol=1e-3, reltol=1e-2)
println("$(sol.retcode) after $(sol.stats.naccept) steps and $(round(time)) s")
dir = joinpath(@__DIR__, "data")
mkpath(dir)
writedlm(joinpath(dir, "solution_u.txt"), reshape(sol.u[end][1:N^2], N, N), ',')
writedlm(joinpath(dir, "solution_v.txt"), reshape(sol.u[end][N^2+1:end], N, N), ',')
