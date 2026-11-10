# The runtime of solving the 2D FitzHugh-Nagumo model to t = 20 at rtol = 0.1, when the
# dimension grows with the domain (Figure 1b), when the problem gets stiffer (1c), and when
# the discretization gets finer (1d).
# Usage: julia -t 16 --project experiments/1_fhn/run.jl domain|stiffness|refinement [method]
include("../../src/setup.jl")

const METHODS = ["EK0", "DiagonalEK1", "EK1", "DiagonalEK1+IEKF", "MatfreeEK1", "MatfreeEK1+IEKF",
                 "KenCarp4+GMRES"]
# The grid sizes N and domain sizes L of each sweep. The dimension is d = 2N².
const SWEEPS = Dict(
    "domain" => [(N, N / 64) for N in [8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256]],
    "stiffness" => [(8, L) for L in [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125]],
    "refinement" => [(N, 1.0) for N in [8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256]],
)
const CUTOFF = 3600.0  # skip the rest of a sweep once a solve takes longer, in seconds

sweep = ARGS[1]
dir = joinpath(@__DIR__, "data", sweep)
mkpath(dir)
for name in (length(ARGS) > 1 ? ARGS[2:2] : METHODS)
    path = joinpath(dir, "$name.csv")
    rm(path; force=true)
    for (N, L) in SWEEPS[sweep]
        prob = fitzhugh_nagumo_2d(; N, L, diag_jac=uses_diag_jac(name))
        prob = uses_diag_jac(name) ? prob : add_jvp_vjp(prob)
        solve_fhn() = solve_once(prob, method(name); abstol=1e-3, reltol=1e-1)
        time = @elapsed sol = solve_fhn()  # includes compilation
        for _ in 1:(time < CUTOFF ? nrepeats(time, 3) : 0)
            time = min(time, @elapsed solve_fhn())
        end
        row = (; method=name, N, L, d=2N^2, dx=L / N, status=string(sol.retcode), time)
        @printf("%-18s N=%-4d L=%-8.4g %-10s time=%.3e s\n", name, N, L, row.status, time)
        flush(stdout)
        CSV.write(path, [row]; append=isfile(path))
        GC.gc()
        time < CUTOFF || break
    end
end
