# The 1D Brusselator at different stiffness coefficients α (Figure 3). The number of steps
# at rtol = 1e-3 for all α (a), and work-precision diagrams at α = 0.02 (b) and α = 10 (c).
# Usage: julia --project experiments/3_alpha_sweep/run.jl 0.02|0.1|0.3|1|3|10 [method]
include("../../src/setup.jl")

const PLOTTED = ["EK0", "DiagonalEK1", "EK1", "ExpEK", "DiagonalEK1+IEKF", "MatfreeEK1",
                 "MatfreeEK1+IEKF"]

# Figure 3 shows no runtimes, so each solve runs only once.
alpha = parse(Float64, ARGS[1])
rtols = alpha in (0.02, 10) ? filter(<=(0.1), RELTOLS) : [1e-3]
run_workprecision(joinpath(@__DIR__, "data", "alpha$(ARGS[1])"),
                  length(ARGS) > 1 ? ARGS[2:2] : PLOTTED,
                  (; kw...) -> brusselator_1d(; alpha, kw...), RadauIIA5(), 1e-12; rtols, numruns=0)
