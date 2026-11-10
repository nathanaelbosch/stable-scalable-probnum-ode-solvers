# Test problems. `diag_jac=true` adds the analytic diagonal of the Jacobian as `f.jac`,
# which the `DiagonalEK1` uses instead of a dense AD Jacobian.

# Brusselator 1D [Hairer, Nørsett, Wanner, Solving ODEs II, p. 6]
# u[1:N] = species 1, u[N+1:2N] = species 2, Dirichlet BCs.
function brusselator_1d!(du, u, p, t)
    A, B, alpha_dx2 = p
    N = length(u) ÷ 2
    @inbounds for i in 2:(N-1)
        u1, u2 = u[i], u[N+i]
        du[i] = alpha_dx2 * (u[i-1] - 2u1 + u[i+1]) + A + u1^2 * u2 - (B + 1) * u1
        du[N+i] = alpha_dx2 * (u[N+i-1] - 2u2 + u[N+i+1]) + B * u1 - u1^2 * u2
    end
    for i in (1, N, N + 1, 2N)
        du[i] = zero(u[i])
    end
    return nothing
end

function brusselator_1d_diag_jac!(J, u, p, t)
    A, B, alpha_dx2 = p
    N = length(u) ÷ 2
    fill!(J.diag, 0)
    @inbounds for i in 2:(N-1)
        u1, u2 = u[i], u[N+i]
        J.diag[i] = -2alpha_dx2 + 2u1 * u2 - (B + 1)
        J.diag[N+i] = -2alpha_dx2 - u1^2
    end
    return J
end

function brusselator_1d(; N=40, A=1.0, B=3.0, alpha=0.02, tspan=(0.0, 10.0), diag_jac=false)
    x = range(0, 1, length=N)
    u0 = [1 .+ sin.(2π .* x); fill(3.0, N)]
    p = (A, B, alpha / (1 / N)^2)
    f = brusselator_1d!
    if diag_jac
        f = ODEFunction(f; jac=brusselator_1d_diag_jac!, jac_prototype=Diagonal(zero(u0)))
    end
    return ODEProblem(f, u0, tspan, p)
end

# Burgers 1D [Bosch et al., Probabilistic Exponential Integrators, NeurIPS 2023]
# ∂ₜu = -u ∂ₓu + D ∂²ₓu on [0, 1], zero-Dirichlet BCs.
function burgers!(du, u, p, t)
    D_dx2, inv_4dx = p
    N = length(u)
    @inbounds for i in 1:N
        ul = i > 1 ? u[i-1] : zero(u[i])
        ur = i < N ? u[i+1] : zero(u[i])
        du[i] = D_dx2 * (ul - 2u[i] + ur) - inv_4dx * (ur^2 - ul^2)
    end
    return nothing
end

function burgers_diag_jac!(J, u, p, t)
    fill!(J.diag, -2p[1])
    return J
end

function burgers(; N=500, D=0.075, tspan=(0.0, 1.0), diag_jac=false)
    dx = 1 / N
    x = range(0, 1, length=N)
    u0 = @. sin(3π * x)^3 * (1 - x)^(3 / 2)
    f = burgers!
    if diag_jac
        f = ODEFunction(f; jac=burgers_diag_jac!, jac_prototype=Diagonal(zero(u0)))
    end
    return ODEProblem(f, u0, tspan, (D / dx^2, 1 / (4dx)))
end

# Fisher-KPP 2D: ∂ₜu = D ∇²u + u(1 - u) on [0, 1]², zero-Neumann BCs.
function fisher_kpp_2d!(du, u, p, t)
    D_dx2, N = p
    U = reshape(u, N, N)
    dU = reshape(du, N, N)
    @inbounds for j in 1:N, i in 1:N
        uij = U[i, j]
        lap = (i > 1 ? U[i-1, j] : uij) + (i < N ? U[i+1, j] : uij) +
              (j > 1 ? U[i, j-1] : uij) + (j < N ? U[i, j+1] : uij) - 4uij
        dU[i, j] = D_dx2 * lap + uij * (1 - uij)
    end
    return nothing
end

function fisher_kpp_2d_diag_jac!(J, u, p, t)
    D_dx2, N = p
    D = reshape(J.diag, N, N)
    U = reshape(u, N, N)
    @inbounds for j in 1:N, i in 1:N
        nneighbors = (i > 1) + (i < N) + (j > 1) + (j < N)
        D[i, j] = -nneighbors * D_dx2 + (1 - 2U[i, j])
    end
    return J
end

function fisher_kpp_2d(; N=64, D=1.0, tspan=(0.0, 1.0), diag_jac=false)
    x = ((1:N) .- 0.5) ./ N
    u0 = vec(@. exp(-30 * (x^2 + x'^2)))
    f = fisher_kpp_2d!
    if diag_jac
        f = ODEFunction(f; jac=fisher_kpp_2d_diag_jac!, jac_prototype=Diagonal(zero(u0)))
    end
    return ODEProblem(f, u0, tspan, (D / (1 / N)^2, N))
end

# Fisher-KPP 1D [Bosch et al., Probabilistic Exponential Integrators, NeurIPS 2023]
# ∂ₜu = D ∂²ₓu + u(1 - u) on [0, 1], zero-Neumann BCs.
function fisher_kpp_1d!(du, u, p, t)
    D_dx2 = p[1]
    N = length(u)
    @inbounds for i in 1:N
        ul = i > 1 ? u[i-1] : u[i]
        ur = i < N ? u[i+1] : u[i]
        du[i] = D_dx2 * (ul + ur - 2u[i]) + u[i] * (1 - u[i])
    end
    return nothing
end

function fisher_kpp_1d_diag_jac!(J, u, p, t)
    D_dx2 = p[1]
    N = length(u)
    @inbounds for i in 1:N
        J.diag[i] = -((i > 1) + (i < N)) * D_dx2 + (1 - 2u[i])
    end
    return J
end

function fisher_kpp_1d(; N=100, D=0.25, tspan=(0.0, 2.0), diag_jac=false)
    x = range(0, 1, length=N)
    u0 = @. 1 - exp(30x - 10) / (1 + exp(30x - 10))
    f = fisher_kpp_1d!
    if diag_jac
        f = ODEFunction(f; jac=fisher_kpp_1d_diag_jac!, jac_prototype=Diagonal(zero(u0)))
    end
    return ODEProblem(f, u0, tspan, (D * N^2,))
end

# Lorenz96 with periodic BCs: du[i] = (u[i+1] - u[i-2]) u[i-1] - u[i] + F.
function lorenz96!(du, u, p, t)
    F, N = p
    @inbounds for i in 1:N
        du[i] = (u[mod1(i + 1, N)] - u[mod1(i - 2, N)]) * u[mod1(i - 1, N)] - u[i] + F
    end
    return nothing
end

lorenz96_diag_jac!(J, u, p, t) = (fill!(J.diag, -1.0); J)

function lorenz96(; N=40, F=8.0, tspan=(0.0, 30.0), diag_jac=false)
    u0 = fill(F, N)
    u0[1] += 0.01
    f = lorenz96!
    if diag_jac
        f = ODEFunction(f; jac=lorenz96_diag_jac!, jac_prototype=Diagonal(zero(u0)))
    end
    return ODEProblem(f, u0, tspan, (F, N))
end

# FitzHugh-Nagumo 2D [Krämer et al., Probabilistic ODE Solutions in Millions of Dimensions,
# ICML 2022]: ∂ₜu = a ∇²u + u - u³ - v + k, ∂ₜv = (b ∇²v + u - v) / τ on [0, L]², with
# zero-Neumann BCs and a random initial condition. u = state[1:N²], v = state[N²+1:2N²].
function fitzhugh_nagumo_2d!(du, u, p, t)
    a_dx2, b_dx2, k, tau, N = p
    N2 = N * N
    U, V = reshape(view(u, 1:N2), N, N), reshape(view(u, N2+1:2N2), N, N)
    dU, dV = reshape(view(du, 1:N2), N, N), reshape(view(du, N2+1:2N2), N, N)
    @inbounds for j in 1:N, i in 1:N
        uij, vij = U[i, j], V[i, j]
        lap_u = (i > 1 ? U[i-1, j] : uij) + (i < N ? U[i+1, j] : uij) +
                (j > 1 ? U[i, j-1] : uij) + (j < N ? U[i, j+1] : uij) - 4uij
        lap_v = (i > 1 ? V[i-1, j] : vij) + (i < N ? V[i+1, j] : vij) +
                (j > 1 ? V[i, j-1] : vij) + (j < N ? V[i, j+1] : vij) - 4vij
        dU[i, j] = a_dx2 * lap_u + uij - uij^3 - vij + k
        dV[i, j] = (b_dx2 * lap_v + uij - vij) / tau
    end
    return nothing
end

function fitzhugh_nagumo_2d_diag_jac!(J, u, p, t)
    a_dx2, b_dx2, k, tau, N = p
    N2 = N * N
    @inbounds for j in 1:N, i in 1:N
        l = (j - 1) * N + i
        nneighbors = (i > 1) + (i < N) + (j > 1) + (j < N)
        J.diag[l] = -nneighbors * a_dx2 + 1 - 3u[l]^2
        J.diag[N2+l] = (-nneighbors * b_dx2 - 1) / tau
    end
    return J
end

function fitzhugh_nagumo_2d(; N=64, L=1.0, a=2.8e-4, b=5e-3, k=-0.005, tau=0.1,
                            tspan=(0.0, 20.0), diag_jac=false)
    u0 = rand(Random.MersenneTwister(2), 2N^2)
    dx = L / N
    f = fitzhugh_nagumo_2d!
    if diag_jac
        f = ODEFunction(f; jac=fitzhugh_nagumo_2d_diag_jac!, jac_prototype=Diagonal(zero(u0)))
    end
    return ODEProblem(f, u0, tspan, (a / dx^2, b / dx^2, k, tau, N))
end
