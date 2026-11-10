# The solvers compared in the paper. The ODE filters use the ProbNumDiffEq defaults (order 3,
# IWP prior, dynamic diffusion, probabilistic error estimate) but don't smooth. The ExpEK
# uses an IOUP prior instead.
const FILTERS = Dict(
    "EK0" => EK0,
    "DiagonalEK1" => DiagonalEK1,
    "EK1" => EK1,
    "ExpEK" => RosenbrockExpEK,
    "DiagonalEK1+IEKF" => (; kw...) -> ODEFilter(; linearization=Iterated(DiagonalJacobian()), kw...),
    "MatfreeEK1" => MatfreeEK1,
    "MatfreeEK1+IEKF" => (; kw...) -> ODEFilter(; linearization=Iterated(MatrixFreeJacobian()), kw...),
)
const CLASSICAL = Dict(
    "Rodas5P" => Rodas5P(),
    "RadauIIA5" => RadauIIA5(),
    "KenCarp4+GMRES" => KenCarp4(linsolve=KrylovJL_GMRES()),
)
method(name; kwargs...) =
    haskey(FILTERS, name) ? FILTERS[name](; smooth=false, kwargs...) : CLASSICAL[name]

# The `DiagonalEK1` variants get the analytic diagonal Jacobian instead of a dense AD Jacobian
uses_diag_jac(name) = startswith(name, "DiagonalEK1")

# Truncate the EK1 covariance to block-diagonal. Done after every step, this is the limit of
# the MatfreeEK1 for many probes. The state is ordered by derivative, so the block of dimension
# i is i:d:D. The QR of R[:, block] gives a square factor of that block of Σ = RᵀR.
function truncate_to_blockdiag!(Σ, d)
    R = Σ.R
    blocks = [qr(R[:, i:d:end]).R for i in 1:d]
    fill!(R, 0)
    for i in 1:d
        R[i:d:end, i:d:end] .= blocks[i]
    end
    return Σ
end
