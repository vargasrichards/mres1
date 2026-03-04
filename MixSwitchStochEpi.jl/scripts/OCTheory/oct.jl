"""
Activity-Structured SEIR Optimal Control
=========================================

Joint optimal control of:
  u_ij(t)  -- contact pattern modulation (NPI kernel)
  v_ij(t)  -- behavioural transition rate perturbations

State:  x = (s_i, e_i, ι_i, r_i) for i = 1..n   [4n-dimensional]
Derived: ρ_i = s_i + e_i + ι_i + r_i

Gradient computation via:
  - Forward pass:  DifferentialEquations.jl (Tsit5 / Rodas5)
  - Adjoint pass:  SciMLSensitivity (InterpolatingAdjoint)
  - Optimisation:  Optim.jl LBFGS
  - LQR subproblem: MatrixEquations.jl

References:
  Tkachenko et al. (2021) PNAS 118:e2015972118
  Novozhilov (2008) Math. Biosci. 215:177-185
"""

using DifferentialEquations
using SciMLSensitivity
using Zygote
using Optim
using LinearAlgebra
using Printf
using Plots

# 1.  PARAMETERS

"""
    SEIRParams

All fixed epidemiological and cost parameters.

Fields
------
- n        : number of activity classes
- a        : contact rates for each class [length n]
- β0       : global transmission scaling
- ϕ        : assortativity parameter ∈ [0,1]
              0 = proportionate mixing, 1 = fully assortative
- σ        : 1/latent period (days⁻¹)
- γ        : 1/infectious period (days⁻¹)
- q0       : baseline behavioural transition rate matrix [n×n]
              q0[i,j] = rate from class i to class j  (diagonal = 0)
- w        : class-specific burden weights [length n]
- α        : cost weight on contact control (NPI cost)
- μ        : cost weight on behavioural transition control
- ν        : terminal penalty on residual E+I
- κ        : terminal penalty on deviation from ρ_target
- ρ_target : desired terminal activity distribution [length n]
- u_min    : lower bound on u_ij controls (scalar, applied uniformly)
"""
struct SEIRParams
    n        :: Int
    a        :: Vector{Float64}
    β0       :: Float64
    ϕ        :: Float64
    σ        :: Float64
    γ        :: Float64
    q0       :: Matrix{Float64}
    w        :: Vector{Float64}
    α        :: Float64
    μ        :: Float64
    ν        :: Float64
    κ        :: Float64
    ρ_target :: Vector{Float64}
    u_min    :: Float64
end

"""
Default parameters: n=3 activity classes (low/medium/high)
"""
function default_params()
    n  = 3
    a  = [0.5, 1.0, 2.5]          # activity rates
    q0 = [0.0   0.02  0.005;      # baseline behavioural transitions
          0.02  0.0   0.02;
          0.005 0.02  0.0  ]
    SEIRParams(
        n, a,
        0.3,          # β0
        0.3,          # ϕ  (moderate assortativity)
        1/3,          # σ  (3-day latent period)
        1/5,          # γ  (5-day infectious period)
        q0,
        [1.0, 1.5, 2.0],   # w  (high-activity classes weigh more)
        5.0,          # α  (NPI cost)
        2.0,          # μ  (behavioural nudge cost)
        10.0,         # ν  (terminal infection penalty)
        1.0,          # κ  (terminal distribution penalty)
        [0.5, 0.35, 0.15],  # ρ_target (shift toward low activity)
        0.05          # u_min (at most 95% contact reduction)
    )
end

# 2.  MIXING KERNEL

"""
    mixing_kernel(p::SEIRParams, ρ::Vector)

Compute the n×n baseline transmission matrix β_ij using the
separable-plus-assortative parameterisation:

    β_ij = β0 [ (1-ϕ) a_i a_j / ā² + ϕ a_i δ_ij / ā ]

where ā = ∑_i a_i ρ_i is the population-weighted mean activity.
"""
function mixing_kernel(p::SEIRParams, ρ::AbstractVector)
    ā  = dot(p.a, ρ)
    β  = zeros(eltype(ρ), p.n, p.n)
    for i in 1:p.n, j in 1:p.n
        proportionate = p.a[i] * p.a[j] / ā^2
        assortative   = (i == j) ? p.a[i] / ā : 0.0
        β[i,j] = p.β0 * ((1 - p.ϕ) * proportionate + p.ϕ * assortative)
    end
    return β
end

# 3.  STATE EQUATIONS

"""
    unpack_state(x, n)

Extract (s, e, ι, r) each of length n from flat state vector x [4n].
"""
function unpack_state(x, n)
    s = x[1:n]
    e = x[n+1:2n]
    ι = x[2n+1:3n]
    r = x[3n+1:4n]
    return s, e, ι, r
end

"""
    pack_state(s, e, ι, r) -> Vector

Flatten compartment vectors into the 4n state vector.
"""
pack_state(s, e, ι, r) = vcat(s, e, ι, r)

"""
    unpack_controls(θ, n)

Extract u_ij (n×n symmetric) and v_ij (n×n, diagonal=0) from
flat parameter vector θ used in optimisation.

Layout of θ:
  [u_upper_triangle (n*(n+1)/2 entries), v_off_diagonal (n*(n-1) entries)]
"""
function unpack_controls(θ::AbstractVector, n::Int)
    nu  = n*(n+1)÷2
    nv  = n*(n-1)
    @assert length(θ) == nu + nv

    # Reconstruct symmetric u matrix from upper triangle
    u = zeros(eltype(θ), n, n)
    k = 1
    for i in 1:n, j in i:n
        u[i,j] = θ[k]
        u[j,i] = θ[k]
        k += 1
    end

    # Reconstruct v matrix (off-diagonal only)
    v = zeros(eltype(θ), n, n)
    k = nu + 1
    for i in 1:n, j in 1:n
        if i != j
            v[i,j] = θ[k]
            k += 1
        end
    end
    return u, v
end

"""
    force_of_infection(β_eff, ι, ρ) -> Vector[n]

λ_i = ∑_j β_eff_ij * ι_j / ρ_j
"""
function force_of_infection(β_eff::AbstractMatrix, ι::AbstractVector, ρ::AbstractVector)
    n = length(ι)
    λ = zeros(eltype(β_eff), n)
    for i in 1:n, j in 1:n
        λ[i] += β_eff[i,j] * ι[j] / (ρ[j] + 1e-12)
    end
    return λ
end

"""
    behavioural_flux(compartment, q) -> Vector[n]

Net influx into class i from behavioural transitions:
    ∑_{j≠i} [q_ji * comp_j - q_ij * comp_i]
"""
function behavioural_flux(comp::AbstractVector, q::AbstractMatrix)
    n = length(comp)
    flux = zeros(eltype(comp), n)
    for i in 1:n, j in 1:n
        if i != j
            flux[i] += q[j,i] * comp[j] - q[i,j] * comp[i]
        end
    end
    return flux
end

"""
    seir_rhs!(dx, x, (p, u_mat, v_mat), t)

In-place RHS for the activity-structured SEIR system with controls.

Arguments
---------
- x     : state vector [4n]
- p     : SEIRParams
- u_mat : n×n contact control matrix  (u_ij ∈ [0,1])
- v_mat : n×n behavioural transition perturbation (q = q0 + v)
"""
function seir_rhs!(dx, x, (p, u_mat, v_mat), t)
    n          = p.n
    s, e, ι, r = unpack_state(x, n)
    ρ          = s .+ e .+ ι .+ r

    # Effective transition rates
    q          = p.q0 .+ v_mat

    # Effective transmission kernel
    β          = mixing_kernel(p, ρ)
    β_eff      = u_mat .* β

    # Force of infection
    λ          = force_of_infection(β_eff, ι, ρ)

    # Behavioural fluxes
    flux_s = behavioural_flux(s, q)
    flux_e = behavioural_flux(e, q)
    flux_ι = behavioural_flux(ι, q)
    flux_r = behavioural_flux(r, q)

    # SEIR derivatives
    ds = -λ .* s .+ flux_s
    de =  λ .* s .- p.σ .* e .+ flux_e
    dι =  p.σ .* e .- p.γ .* ι .+ flux_ι
    dr =  p.γ .* ι .+ flux_r

    dx .= pack_state(ds, de, dι, dr)
    return nothing
end

# 4.  OBJECTIVE FUNCTION

"""
    running_cost(x, u_mat, v_mat, p) -> scalar

L(x, u, v) = ∑_i w_i ι_i + (α/2)∑_{i≤j}(1-u_ij)² + (μ/2)∑_{i≠j} v_ij²
"""
function running_cost(x::AbstractVector, u_mat::AbstractMatrix,
                      v_mat::AbstractMatrix, p::SEIRParams)
    n          = p.n
    _, _, ι, _ = unpack_state(x, n)

    burden_cost = dot(p.w, ι)

    npi_cost = 0.0
    for i in 1:n, j in i:n
        npi_cost += (1 - u_mat[i,j])^2
    end
    npi_cost *= p.α / 2

    beh_cost = 0.0
    for i in 1:n, j in 1:n
        if i != j
            beh_cost += v_mat[i,j]^2
        end
    end
    beh_cost *= p.μ / 2

    return burden_cost + npi_cost + beh_cost
end

"""
    terminal_cost(x_T, p) -> scalar

Φ(x(T)) = ν ∑_i (e_i+ι_i) + κ ∑_i (ρ_i - ρ*_i)²
"""
function terminal_cost(x_T::AbstractVector, p::SEIRParams)
    n          = p.n
    s, e, ι, r = unpack_state(x_T, n)
    ρ          = s .+ e .+ ι .+ r

    inf_penalty  = p.ν * sum(e .+ ι)
    dist_penalty = p.κ * sum((ρ .- p.ρ_target).^2)

    return inf_penalty + dist_penalty
end

# 5.  FORWARD SOLVE AND LOSS (Zygote-differentiable)

"""
    solve_forward(x0, θ, p, tspan; dt=1.0)

Solve the SEIR ODE forward in time with controls encoded in θ.
Controls are CONSTANT over the horizon (bang-bang relaxation baseline);
for time-varying controls replace θ with a neural network / spline.

Returns: ODE solution object.
"""
function solve_forward(x0, θ, p::SEIRParams, tspan; dt=1.0)
    u_mat, v_mat = unpack_controls(θ, p.n)

    # Project controls to admissible set
    u_mat = clamp.(u_mat, p.u_min, 1.0)
    v_mat = max.(v_mat, -p.q0)   # ensure q = q0 + v >= 0

    prob = ODEProblem(seir_rhs!, x0, tspan, (p, u_mat, v_mat))
    sol  = solve(prob, Tsit5();
                 saveat  = tspan[1]:dt:tspan[2],
                 abstol  = 1e-8,
                 reltol  = 1e-6,
                 sensealg = InterpolatingAdjoint(autojacvec=ZygoteVJP()))
    return sol
end

"""
    compute_loss(x0, θ, p, tspan; dt=1.0)

Full objective J = ∫L dt + Φ(x(T)), computed via quadrature over
the saved time points of the forward solve.

This function is differentiable with respect to θ via Zygote.
"""
function compute_loss(x0, θ, p::SEIRParams, tspan; dt=1.0)
    u_mat, v_mat = unpack_controls(θ, p.n)
    u_mat_proj   = clamp.(u_mat, p.u_min, 1.0)
    v_mat_proj   = max.(v_mat, -p.q0)

    sol = solve_forward(x0, θ, p, tspan; dt=dt)

    # Trapezoidal integration of running cost
    ts     = sol.t
    xs     = sol.u
    costs  = [running_cost(x, u_mat_proj, v_mat_proj, p) for x in xs]
    J_run  = sum(0.5 * (costs[k] + costs[k+1]) * (ts[k+1] - ts[k])
                 for k in 1:(length(ts)-1))

    J_term = terminal_cost(last(xs), p)
    return J_run + J_term
end

# 6.  GRADIENT VIA ZYGOTE AND OPTIMISATION

"""
    optimise_controls(x0, p, tspan;
                      dt       = 1.0,
                      n_iter   = 200,
                      verbose  = true)

Solve the optimal control problem using LBFGS with Zygote gradients.

Returns: (θ_opt, loss_history, sol_opt)
"""
function optimise_controls(x0, p::SEIRParams, tspan;
                           dt      = 1.0,
                           n_iter  = 200,
                           verbose = true)
    n  = p.n
    nu = n*(n+1)÷2
    nv = n*(n-1)

    # Initial guess: u = 1 (no NPI), v = 0 (no behavioural nudge)
    θ0 = vcat(ones(nu), zeros(nv))

    loss_history = Float64[]

    # Objective and gradient for Optim
    function fg!(F, G, θ)
        val, grads = Zygote.withgradient(θ_ -> compute_loss(x0, θ_, p, tspan; dt=dt), θ)
        if G !== nothing
            G .= grads[1]
        end
        if F !== nothing
            push!(loss_history, val)
            if verbose && length(loss_history) % 10 == 0
                @printf("  iter %4d  |  loss = %.6f\n",
                        length(loss_history), val)
            end
            return val
        end
    end

    result = optimize(Optim.only_fg!(fg!), θ0,
                      LBFGS(),
                      Optim.Options(iterations=n_iter,
                                    show_trace=false,
                                    store_trace=true))

    θ_opt   = Optim.minimizer(result)
    sol_opt = solve_forward(x0, θ_opt, p, tspan; dt=dt)

    return θ_opt, loss_history, sol_opt
end

# 7.  COSTATE / ADJOINT EQUATIONS (explicit, for verification)

"""
    costate_rhs!(dψ, ψ, (x_interp, u_mat, v_mat, p), t)

Explicit adjoint equations (run BACKWARDS from T to 0).
ψ = (ψ_s, ψ_e, ψ_ι, ψ_r) each of length n.

Useful for:
 (a) Verifying SciML adjoint gradients
 (b) Computing the optimal control analytically given ψ
"""
function costate_rhs!(dψ, ψ, (x_interp, u_mat, v_mat, p), t)
    n  = p.n
    x  = x_interp(t)
    s, e, ι, r = unpack_state(x, n)
    ρ  = s .+ e .+ ι .+ r

    ψ_s, ψ_e, ψ_ι, ψ_r = unpack_state(ψ, n)

    q     = p.q0 .+ v_mat
    β     = mixing_kernel(p, ρ)
    β_eff = u_mat .* β
    λ     = force_of_infection(β_eff, ι, ρ)

    Q_out = [sum(q[i,j] for j in 1:n if j != i) for i in 1:n]

    dψ_s = zeros(n)
    dψ_e = zeros(n)
    dψ_ι = zeros(n)
    dψ_r = zeros(n)

    for i in 1:n
        # ψ_s equation: ∂H/∂s_i
        dψ_s[i] = ( λ[i] * (ψ_s[i] - ψ_e[i])
                  + ψ_s[i] * Q_out[i]
                  - sum(q[j,i] * ψ_s[j] for j in 1:n if j != i) )

        # ψ_e equation: ∂H/∂e_i
        dψ_e[i] = ( p.σ * (ψ_e[i] - ψ_ι[i])
                  + ψ_e[i] * Q_out[i]
                  - sum(q[j,i] * ψ_e[j] for j in 1:n if j != i) )

        # ψ_ι equation: ∂H/∂ι_i  (includes coupling through λ_k)
        coupling = sum( (ψ_s[k] - ψ_e[k]) * s[k] *
                        β_eff[k,i] / (ρ[i] + 1e-12)
                        for k in 1:n )
        dψ_ι[i] = ( -p.w[i]
                  + p.γ * (ψ_ι[i] - ψ_r[i])
                  + ψ_ι[i] * Q_out[i]
                  - sum(q[j,i] * ψ_ι[j] for j in 1:n if j != i)
                  + coupling )

        # ψ_r equation: ∂H/∂r_i
        dψ_r[i] = ( ψ_r[i] * Q_out[i]
                  - sum(q[j,i] * ψ_r[j] for j in 1:n if j != i) )
    end

    dψ .= pack_state(dψ_s, dψ_e, dψ_ι, dψ_r)
    return nothing
end

"""
    solve_costate(sol_fwd, u_mat, v_mat, p)

Solve the costate equations backwards from T.
Terminal conditions:
  ψ_s(T) = 0,  ψ_e(T) = ν,  ψ_ι(T) = ν,  ψ_r(T) = 0.
"""
function solve_costate(sol_fwd, u_mat, v_mat, p::SEIRParams)
    n     = p.n
    tspan = (sol_fwd.t[end], sol_fwd.t[1])

    x_interp = t -> sol_fwd(t)

    ψ0 = pack_state(zeros(n), fill(p.ν, n), fill(p.ν, n), zeros(n))

    prob = ODEProblem(costate_rhs!, ψ0, tspan,
                      (x_interp, u_mat, v_mat, p))
    sol  = solve(prob, Tsit5(); saveat=0.5, abstol=1e-8, reltol=1e-6)
    return sol
end

# 8.  OPTIMAL CONTROL FROM COSTATES (analytical, for verification)

"""
    optimal_u_from_costate(ψ_s, ψ_e, ι, ρ, β, p) -> Matrix

Compute u*_ij = max(u_min, min(1, ũ_ij)) from costate variables,
using the first-order condition:

  ũ_ij = 1 + (β_ij/α)[ (ψ_e_i - ψ_s_i) s_i ι_j/ρ_j
                       + (ψ_e_j - ψ_s_j) s_j ι_i/ρ_i ]
"""
function optimal_u_from_costate(ψ_s, ψ_e, s, ι, ρ, β, p::SEIRParams)
    n = p.n
    u = ones(n, n)
    for i in 1:n, j in 1:n
        gap_i = (ψ_e[i] - ψ_s[i]) * s[i] * ι[j] / (ρ[j] + 1e-12)
        gap_j = (ψ_e[j] - ψ_s[j]) * s[j] * ι[i] / (ρ[i] + 1e-12)
        ũ = 1.0 + (β[i,j] / p.α) * (gap_i + gap_j)
        u[i,j] = clamp(ũ, p.u_min, 1.0)
    end
    return u
end

"""
    optimal_v_from_costate(ψ_s, ψ_e, ψ_ι, ψ_r, s, e, ι, r, p) -> Matrix

Compute v*_ij from first-order condition:

  v*_ij = -(1/μ)[ ψ_s_i s_i + ψ_e_i e_i + ψ_ι_i ι_i + ψ_r_i r_i
                  - ψ_s_j s_j - ψ_e_j e_j - ψ_ι_j ι_j - ψ_r_j r_j ]

projected to ensure q0 + v >= 0.
"""
function optimal_v_from_costate(ψ_s, ψ_e, ψ_ι, ψ_r, s, e, ι, r, p::SEIRParams)
    n       = p.n
    shadow  = [ψ_s[i]*s[i] + ψ_e[i]*e[i] + ψ_ι[i]*ι[i] + ψ_r[i]*r[i]
               for i in 1:n]
    v = zeros(n, n)
    for i in 1:n, j in 1:n
        if i != j
            v_raw    = -(shadow[i] - shadow[j]) / p.μ
            v[i,j]   = max(v_raw, -p.q0[i,j])   # admissibility
        end
    end
    return v
end

# 9.  STEADY-STATE ANALYSIS

"""
    R0(p::SEIRParams) -> Float64

Compute the basic reproduction number as the spectral radius of
the next-generation matrix K_ij = β_ij * s_j(0) / (γ ρ_j(0)).
"""
function compute_R0(x0, p::SEIRParams)
    s, e, ι, r = unpack_state(x0, p.n)
    ρ  = s .+ e .+ ι .+ r
    β  = mixing_kernel(p, ρ)
    K  = [β[i,j] * s[j] / (p.γ * (ρ[j] + 1e-12)) for i in 1:p.n, j in 1:p.n]
    return maximum(abs.(eigvals(K)))
end

# 10. LQR FLUCTUATION CONTROL (requires MatrixEquations.jl)

"""
    lqr_gains(x_star, u_star, v_star, p)

Compute the optimal LQR feedback gains around endemic equilibrium x*.
Returns: (P, K_u, K_v) where P solves the algebraic Riccati equation
and the optimal policy is:
  δu = -R_u⁻¹ B_u' P δx
  δv = -R_v⁻¹ B_v' P δx

Uses ForwardDiff to compute Jacobians A, B_u, B_v automatically.
"""
function lqr_gains(x_star, u_star_vec, v_star_vec, p::SEIRParams)
    # This function requires:  using MatrixEquations, ForwardDiff
    # Uncomment and add to Project.toml as needed
    #
    # n  = p.n
    # N  = 4n
    # nu = n*(n+1)÷2
    # nv = n*(n-1)
    #
    # Pack controls into flat vector for Jacobian
    # θ_star = vcat(u_star_vec, v_star_vec)
    #
    # # Jacobian A = ∂f/∂x at (x*, u*, v*)
    # f_x(x_) = begin
    #     dx = similar(x_)
    #     seir_rhs!(dx, x_, (p, reshape_u(u_star_vec,n), reshape_v(v_star_vec,n)), 0.0)
    #     dx
    # end
    # A = ForwardDiff.jacobian(f_x, x_star)
    #
    # # Jacobians B_u, B_v via ForwardDiff through θ
    # f_θ(θ_) = begin
    #     dx = similar(x_star)
    #     u_, v_ = unpack_controls(θ_, n)
    #     seir_rhs!(dx, x_star, (p, u_, v_), 0.0)
    #     dx
    # end
    # B = ForwardDiff.jacobian(f_θ, θ_star)
    # B_u = B[:, 1:nu]
    # B_v = B[:, nu+1:end]
    #
    # # LQR cost matrices
    # Q = Diagonal(vcat(zeros(n), zeros(n), p.w, zeros(n)))
    # R = Diagonal(vcat(fill(p.α, nu), fill(p.μ, nv)))
    #
    # # Solve algebraic Riccati:  A'P + PA - PBR⁻¹B'P + Q = 0
    # P = ared(A, B, R, Q).X   # MatrixEquations.ared
    # K = R \ (B' * P)
    # K_u = K[1:nu, :]
    # K_v = K[nu+1:end, :]
    # return P, K_u, K_v

    @warn "lqr_gains requires MatrixEquations.jl and ForwardDiff.jl. Stub returned."
    return nothing, nothing, nothing
end

# 11.  INITIAL CONDITIONS AND EXAMPLE RUN

"""
    initial_conditions(p::SEIRParams; ι0_frac=1e-4)

Construct x0 with a small seed of infectious individuals in the
highest-activity class and near-fully-susceptible population.
"""
function initial_conditions(p::SEIRParams; ι0_frac=1e-4)
    n  = p.n
    ρ0 = [0.5, 0.35, 0.15]       # initial activity distribution
    @assert length(ρ0) == n

    # Seed infections in highest-activity class only
    ι0 = zeros(n); ι0[end] = ι0_frac * ρ0[end]
    s0 = ρ0 .- ι0
    e0 = zeros(n)
    r0 = zeros(n)

    return pack_state(s0, e0, ι0, r0)
end

# 12.  PLOTTING UTILITIES

"""
    plot_epidemic(sol, p; title="SEIR Dynamics")

Plot total prevalence ∑_i ι_i(t) and class-stratified susceptibles
over the simulation horizon.
"""
function plot_epidemic(sol, p::SEIRParams; title="SEIR Dynamics")
    n  = p.n
    ts = sol.t
    xs = hcat(sol.u...)'    # rows = timepoints, cols = state dims

    total_I = [sum(xs[k, 2n+1:3n]) for k in 1:length(ts)]
    total_S = [sum(xs[k, 1:n])     for k in 1:length(ts)]
    total_R = [sum(xs[k, 3n+1:4n]) for k in 1:length(ts)]

    p1 = plot(ts, total_I .* 100, label="I (total %)", lw=2, colour=:red)
    plot!(p1, ts, total_S .* 100, label="S (total %)", lw=2, colour=:blue)
    plot!(p1, ts, total_R .* 100, label="R (total %)", lw=2, colour=:green)
    xlabel!(p1, "Time (days)")
    ylabel!(p1, "Population fraction (%)")
    title!(p1, title)

    # Class-stratified infectious
    p2 = plot()
    labels = ["Low activity", "Medium activity", "High activity"]
    colours = [:skyblue, :orange, :crimson]
    for i in 1:n
        plot!(p2, ts, xs[:, 2n+i] .* 100, label=labels[i],
              lw=2, colour=colours[i])
    end
    xlabel!(p2, "Time (days)")
    ylabel!(p2, "Infectious (%)")
    title!(p2, "Class-stratified prevalence")

    return plot(p1, p2, layout=(2,1), size=(800, 600))
end

"""
    plot_controls(θ_opt, p)

Bar chart of optimal u_ij and v_ij values.
"""
function plot_controls(θ_opt, p::SEIRParams)
    n            = p.n
    u_mat, v_mat = unpack_controls(θ_opt, n)
    u_mat        = clamp.(u_mat, p.u_min, 1.0)
    v_mat        = max.(v_mat, -p.q0)

    labels_u = ["u$i$j" for i in 1:n for j in i:n]
    vals_u   = [u_mat[i,j] for i in 1:n for j in i:n]

    labels_v = ["v$i→$j" for i in 1:n for j in 1:n if i != j]
    vals_v   = [v_mat[i,j] for i in 1:n for j in 1:n if i != j]

    pu = bar(labels_u, vals_u, title="Contact control u_ij*",
             ylims=(0,1.05), ylabel="Suppression factor", legend=false,
             colour=:steelblue, rotation=45)

    pv = bar(labels_v, vals_v, title="Behavioural transition control v_ij*",
             ylabel="Added transition rate (days⁻¹)", legend=false,
             colour=:darkorange, rotation=45)

    return plot(pu, pv, layout=(1,2), size=(900, 400))
end

# 13.  MAIN SCRIPT

function main()
    println("="^60)
    println("  Activity-Structured SEIR Optimal Control")
    println("="^60)

    p     = default_params()
    x0    = initial_conditions(p; ι0_frac=1e-4)
    tspan = (0.0, 150.0)

    R0_val = compute_R0(x0, p)
    @printf("\nBasic reproduction number R0 = %.3f\n\n", R0_val)

    # --- Uncontrolled baseline ---
    println("Solving uncontrolled baseline...")
    θ_baseline = vcat(ones(p.n*(p.n+1)÷2), zeros(p.n*(p.n-1)))
    sol_base   = solve_forward(x0, θ_baseline, p, tspan; dt=1.0)
    attack_base = 1.0 - sum(sol_base.u[end][1:p.n])
    @printf("  Uncontrolled attack rate: %.1f%%\n\n", attack_base * 100)

    # --- Optimal control ---
    println("Optimising controls (LBFGS)...\n")
    θ_opt, loss_hist, sol_opt = optimise_controls(x0, p, tspan;
                                                   dt=1.0,
                                                   n_iter=100,
                                                   verbose=true)

    attack_opt = 1.0 - sum(sol_opt.u[end][1:p.n])
    @printf("\nOptimal attack rate: %.1f%%\n", attack_opt * 100)
    @printf("Attack rate reduction: %.1f percentage points\n\n",
            (attack_base - attack_opt) * 100)

    # --- Extract and display optimal controls ---
    u_opt, v_opt = unpack_controls(θ_opt, p.n)
    u_opt = clamp.(u_opt, p.u_min, 1.0)
    v_opt = max.(v_opt, -p.q0)

    println("Optimal contact suppression matrix u*:")
    display(round.(u_opt, digits=3))
    println("\nOptimal behavioural transition perturbation v*:")
    display(round.(v_opt, digits=4))

    # --- Plots ---
    fig1 = plot_epidemic(sol_base, p; title="Uncontrolled SEIR")
    fig2 = plot_epidemic(sol_opt,  p; title="Optimally Controlled SEIR")
    fig3 = plot_controls(θ_opt, p)

    savefig(fig1, "seir_uncontrolled.png")
    savefig(fig2, "seir_controlled.png")
    savefig(fig3, "optimal_controls.png")
    println("\nFigures saved: seir_uncontrolled.png, seir_controlled.png, optimal_controls.png")

    return θ_opt, sol_base, sol_opt, loss_hist
end

# Uncomment to run:
θ_opt, sol_base, sol_opt, loss_hist = main()