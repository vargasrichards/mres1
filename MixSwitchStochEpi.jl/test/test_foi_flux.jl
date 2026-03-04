# Stochastic FOI flux test
# Run many short simulations with no progression/recovery and no switching
# and verify mean number of S->E events ≈ rate_inf_at_t0 * tmax (small-time approx)

using Test
using Random
using Statistics
using MixSwitchStochEpi

@testset "FOI small-time flux (stochastic)" begin
    # small system to keep tests quick
    n = 3
    class_sizes = [10, 10, 10]
    pop = sum(class_sizes)
    act_levels = [1.0, 1.0, 1.0]

    # uniform mixing (rows equal) -> M[i,j] = 1/3
    M = fill(1/3, n, n)
    W = zeros(Float64, n, n)

    # pick parameters so there's a moderate infection rate at t=0
    β = 5.0
    σ = 0.0   # no E->I progression
    γ = 0.0   # no recovery
    ω = 0.0

    params = SEIRSStochParams(β, σ, γ, ω, M, W, n, class_sizes, act_levels, pop)

    # initial state: place several infectious individuals in class 1
    I = [6, 0, 0]
    Ivec = I
    S = [class_sizes[i] - Ivec[i] for i in 1:n]
    E = zeros(Int, n)
    R = zeros(Int, n)
    initial_state = vcat(Int.(S)', Int.(E)', Int.(Ivec)', Int.(R)')

    # compute initial total infection rate using same formula as simulator
    function initial_rate(params, state)
        n = params.n_activity
        N = [Float64(sum(state[:, j])) for j in 1:n]
        λ = zeros(Float64, n)
        for i in 1:n
            denom = 0.0
            numer = 0.0
            for j in 1:n
                w = params.M[i, j] * N[j]
                denom += w
                if N[j] > 0
                    numer += w * (state[3, j] / N[j])
                end
            end
            λ[i] = denom > 0.0 ? (params.β * params.act_levels[i] * numer / denom) : 0.0
        end
        rate_inf = sum(λ[i] * state[1, i] for i in 1:n)
        return rate_inf
    end

    rate0 = initial_rate(params, initial_state)

    # small time horizon so that probability of multiple infections per trajectory is small
    tmax = 0.02

    # expected mean number of S->E events (first-order small-time approximation)
    expected_mean = rate0 * tmax

    # run many short simulations and record the number of new E at tmax
    rng = MersenneTwister(1234)
    n_sims = 800
    counts = zeros(Int, n_sims)
    for sim in 1:n_sims
        res = simulate_system!(params, initial_state; tmax = tmax, rng = rng, store_states = true)
        final_state = res.states[end]
        # number of S->E events = final E total - initial E total
        counts[sim] = sum(final_state[2, :]) - sum(initial_state[2, :])
    end

    observed_mean = mean(counts)

    # Allow some sampling variability but be strict enough to catch regressions
    @test isapprox(observed_mean, expected_mean; rtol = 0.2, atol = 0.05)
end
