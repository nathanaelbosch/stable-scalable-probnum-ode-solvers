# The runtime of a single step on Lorenz96, against the dimension d (Figure 2).
# Usage: julia -t 16 --project experiments/2_dim_scaling/run.jl [method]
include("../../src/setup.jl")

# Methods and their expected cost per step, as a power of d
const PLOTTED = ["EK0" => 1, "DiagonalEK1" => 1, "EK1" => 3, "ExpEK" => 3, "DiagonalEK1+IEKF" => 1,
                 "MatfreeEK1" => 1, "MatfreeEK1+IEKF" => 1]
const DIMS = 2 .^ (3:21)
const CUTOFF = 300.0  # skip a method's larger d once its next step is predicted to take longer

function time_step(prob, alg)
    integ = init(prob, alg; adaptive=false, dt=0.01, save_everystep=false)
    step!(integ)
    return minimum(1:5) do _
        reinit!(integ)
        @elapsed step!(integ)
    end
end

dir = joinpath(@__DIR__, "data")
mkpath(dir)
# A trajectory at d = 40 for the plot
sol = solve(lorenz96(), RadauIIA5(); abstol=1e-9, reltol=1e-9, saveat=0.05)
writedlm(joinpath(dir, "trajectory.txt"), reduce(hcat, sol.u), ',')

for (name, slope) in (length(ARGS) > 0 ? filter(p -> p[1] == ARGS[1], PLOTTED) : PLOTTED)
    path = joinpath(dir, "$name.csv")
    rm(path; force=true)
    for d in DIMS
        prob = uses_diag_jac(name) ? lorenz96(N=d, diag_jac=true) : add_jvp_vjp(lorenz96(N=d))
        time = time_step(prob, method(name))
        @printf("%-18s d=%-8d time=%.3e s\n", name, d, time)
        flush(stdout)
        CSV.write(path, [(; method=name, d, time)]; append=isfile(path))
        GC.gc()
        time * 2^slope < CUTOFF || break
    end
end
