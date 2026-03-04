# FOI zero-activity tests
# Ensure classes with zero activity score experience zero force-of-infection

using Test
using Random
using Statistics
using LinearAlgebra
using MixSwitchStochEpi

function compute_lambda(params, state)
    n = params.n_activity
    N = [Float64(sum(state[:, j])) for j in 1:n]
    λ = zeros(Float64, n)
    for i in 1:n
        denom = 0.0
        numer = 0.0
        for j in 1:n
            w = params.M[i, j] * N[j]
            denom += w
            if N[j] > 0.0
                numer += w * (state[3, j] / N[j])
            end
        end
        λ[i] = denom > 0.0 ? (params.β * params.act_levels[i] * numer / denom) : 0.0
    end
    return λ
end

@testset "FOI zero activity deterministic and stochastic" begin
    rng = MersenneTwister(20260219)
    n = 3
    class_sizes = [10, 10, 10]
    # middle class has zero activity
    act_levels = [1.0, 0.0, 2.0]
    pop = sum(class_sizes)

    # initial infections placed in classes 1 and 3 only
    Iseed = [1, 0, 2]
    Sstart = [class_sizes[i] - Iseed[i] for i in 1:n]
    Estart = zeros(Int, n)
    Rstart = zeros(Int, n)
    initial_state = vcat(Int.(Sstart)', Int.(Estart)', Int.(Iseed)', Int.(Rstart)')

    β = 1.0; σ = 0.0; γ = 0.0; ω = 0.0

    # test a few mixing matrices
    mix_list = [
        ("identity", Matrix{Float64}(I, n, n)),
        ("proportionate", make_garnett_contact(0.0, n, act_levels, class_sizes)),
        ("mixed", make_garnett_contact(0.4, n, act_levels, class_sizes)),
    ]

    for (name, M) in mix_list
        @testset "mix=$(name)" begin
            W = zeros(Float64, n, n)
            params = SEIRSStochParams(β, σ, γ, ω, M, W, n, class_sizes, act_levels, pop)

            λ = compute_lambda(params, initial_state)
            @test isapprox(λ[2], 0.0; atol = 1e-14)

            # stochastic check: run many short sims and ensure no S->E occurs in class 2
            n_sims = 1000
            # record new E in class 2 for each sim
            deltaE2 = zeros(Int, n_sims)
            for sim in 1:n_sims
                res = simulate_system!(params, initial_state; tmax = 1.0, rng = rng, store_states = true)
                deltaE = res.states[end][2, :] .- initial_state[2, :]
                deltaE2[sim] = deltaE[2]
            end

            @test all(deltaE2 .== 0)
        end
    end
end
