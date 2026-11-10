# Work-precision diagrams on Burgers' equation and the 2D Fisher-KPP equation (Figure 5).
# Usage: julia -t 16 --project experiments/5_workprecision/run.jl burgers|fisher_kpp_2d [method]
# The EK1 on Fisher-KPP only at the loosest tolerance, for the runtime in the text:
# RTOLS=10 julia -t 16 --project experiments/5_workprecision/run.jl fisher_kpp_2d EK1
include("../../src/setup.jl")

const PLOTTED = ["DiagonalEK1", "EK1", "DiagonalEK1+IEKF", "MatfreeEK1", "MatfreeEK1+IEKF"]
# Classical solvers for comparison, in the appendix
const BASELINES = ["Rodas5P", "RadauIIA5", "KenCarp4+GMRES"]
const PROBLEMS = Dict(
    "burgers" => (; make=(; kw...) -> burgers(; N=500, kw...), ref_alg=RadauIIA5(), ref_tol=1e-12,
                  methods=[PLOTTED; BASELINES]),
    # The EK1 is too slow here. It takes about 26 minutes already at the loosest tolerance.
    "fisher_kpp_2d" => (; make=(; kw...) -> fisher_kpp_2d(; N=64, kw...), ref_alg=Rodas5P(),
                        ref_tol=1e-10, methods=[filter(!=("EK1"), PLOTTED); BASELINES]),
)

key = ARGS[1]
(; make, ref_alg, ref_tol, methods) = PROBLEMS[key]
run_workprecision(joinpath(@__DIR__, "data", key), length(ARGS) > 1 ? ARGS[2:2] : methods,
                  make, ref_alg, ref_tol)
