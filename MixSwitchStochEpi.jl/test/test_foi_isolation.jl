# FOI isolation tests — ensure production FOI in `simulate_system!` is used
# Check that the class receiving the first infection follows the analytic probabilities
# computed from the same FOI formula used in the simulator.

using Test
using Random
using Statistics
using LinearAlgebra
using MixSwitchStochEpi

function analytic_first_infection_probs(params, state)
    n = params.n_activity
    N = [Float64(sum(state[:, j])) for j in 1:n]
    Ivec = [state[3, j] for j in 1:n]
    Svec = [state[1, j] for j in 1:n]

    # compute λ as in simulate_system!
    λ = zeros(Float64, n)
    for i in 1:n
        denom = 0.0
        numer = 0.0
        for j in 1:n
            w = params.M[i, j] * N[j]
            denom += w
            if N[j] > 0.0
                numer += w * (Ivec[j] / N[j])
            end
        end
        λ[i] = denom > 0.0 ? (params.β * params.act_levels[i] * numer / denom) : 0.0
    end

    weights = [λ[i] * Svec[i] for i in 1:n]
    tot = sum(weights)
    if tot == 0.0
        return fill(1.0 / n, n)
    else
        return weights ./ tot
    end
end

@testset "FOI isolation: first infection class distribution" begin
    rng = MersenneTwister(20260218)
    n = 4
    class_sizes = [10, 20, 30, 40]
    act_levels = [1.0, 1.5, 2.0, 0.5]
    pop = sum(class_sizes)

    # Infectious seed in multiple classes to exercise mixing
    Iseed = [2, 1, 3, 0]
    Sstart = [class_sizes[i] - Iseed[i] for i in 1:n]
    Estart = zeros(Int, n)
    Rstart = zeros(Int, n)
    initial_state = vcat(Int.(Sstart)', Int.(Estart)', Int.(Iseed)', Int.(Rstart)')

    β = 1.5
    σ = 0.0; γ = 0.0; ω = 0.0

    # Test across several mixing patterns
    mix_set = Dict(
        :identity => Matrix{Float64}(I, n, n),
        :proportionate => make_garnett_contact(0.0, n, act_levels, class_sizes),
        :mixed => make_garnett_contact(0.4, n, act_levels, class_sizes),
    )

    n_sims = 2000
    i = 0
    for (name, M) in mix_set
        i += 1
        @testset "mix: $(name)" begin
            # reset RNG for each mixing case for independent sampling
            case_rng = MersenneTwister(20260218 + i)
            W = zeros(Float64, n, n) # no switching
            params = SEIRSStochParams(β, σ, γ, ω, M, W, n, class_sizes, act_levels, pop)

            probs = analytic_first_infection_probs(params, initial_state)

            counts = zeros(Int, n)
            for sim in 1:n_sims
                res = simulate_system!(params, initial_state; tmax = Inf, rng = case_rng, store_states = true)
                # find first post-initial state (should be states[2])
                @assert length(res.states) >= 2
                first_after = res.states[2]
                deltaE = first_after[2, :] .- initial_state[2, :]
                # locate the class where E increased by 1
                idx = findfirst(==(1), deltaE)
                @assert idx !== nothing
                counts[idx] += 1
            end

            observed = counts ./ sum(counts)

            maxabs = maximum(abs.(observed .- probs))
            @test maxabs <= 0.05
        end
    end
end
