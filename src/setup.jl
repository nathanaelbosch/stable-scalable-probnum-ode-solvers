using ProbNumDiffEq, OrdinaryDiffEq, OrdinaryDiffEqFIRK, OrdinaryDiffEqSDIRK
using LinearSolve: KrylovJL_GMRES
using LinearAlgebra, Statistics, Printf, Random
import DifferentiationInterface as DI
using ADTypes: AutoForwardDiff, AutoEnzyme
import ForwardDiff, Enzyme
using CSV, DelimitedFiles
import SciMLBase

# Backport of JuliaLang/LinearAlgebra.jl#1814, not yet in a Julia release. Without it, dividing
# by a diagonal Cholesky factor runs a generic single-threaded solve. At d = 2048 the EK1
# initialization then took 2177 s instead of 6 s.
LinearAlgebra.rdiv!(A::AbstractMatrix, B::LinearAlgebra.UpperOrLowerTriangular) =
    LinearAlgebra._rdiv!(A, A, B)
LinearAlgebra.ldiv!(A::LinearAlgebra.UpperOrLowerTriangular, B::AbstractVecOrMat) =
    LinearAlgebra._ldiv!(B, A, B)

include("problems.jl")
include("methods.jl")
include("workprecision.jl")

# Add the Jacobian-vector and vector-Jacobian products that the `MatfreeEK1` needs
function add_jvp_vjp(prob)
    f, u, p, t = prob.f.f, prob.u0, prob.p, prob.tspan[1]
    du, v = similar(u), similar(u)
    jvp_backend, vjp_backend = AutoForwardDiff(), AutoEnzyme()
    jvp_prep = DI.prepare_pushforward(f, du, jvp_backend, u, (v,), DI.Constant(p), DI.Constant(t))
    vjp_prep = DI.prepare_pullback(f, du, vjp_backend, u, (v,), DI.Constant(p), DI.Constant(t))
    jvp(Jv, v, u, p, t) = DI.pushforward!(
        f, du, (Jv,), jvp_prep, jvp_backend, u, (v,), DI.Constant(p), DI.Constant(t))
    vjp(Jᵀv, v, u, p, t) = DI.pullback!(
        f, du, (Jᵀv,), vjp_prep, vjp_backend, u, (v,), DI.Constant(p), DI.Constant(t))
    return remake(prob; f=ODEFunction{true,SciMLBase.FullSpecialize}(f; jvp, vjp))
end
