# Deterministic FOI snapshot test
# Alexis Vargas Richards, Feb 2026

using Test
using LinearAlgebra
using MixSwitchStochEpi

@testset "FOI snapshot equals analytic formula" begin
    # small 3-class system where M is identity so FOI reduces to β*c_i*I_i/N_i
    n = 3
    class_sizes = [10, 20, 30]
    act_levels = [1.0, 2.0, 3.0]
    pop = sum(class_sizes)

    # contact matrix: identity (purely assortative)
    M = Matrix{Float64}(I, n, n)
    # no switching for this deterministic check
    W = zeros(Float64, n, n)

    β = 0.25
    σ = 1/3
    γ = 1/4
    ω = 0.0

    params = SEIRSStochParams(β, σ, γ, ω, M, W, n, class_sizes, act_levels, pop)

    # states: S, E, I, R (4 x n)
    S = [9, 18, 27]
    E = zeros(Int, n)
    Ivec = [1, 2, 3]
    R = zeros(Int, n)
    state = vcat(Int.(S)', Int.(E)', Int.(Ivec)', Int.(R)')

    # compute FOI following the same algorithm as the simulator
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
                if N[j] > 0
                    numer += w * (state[3, j] / N[j])
                end
            end
            λ[i] = denom > 0.0 ? (params.β * params.act_levels[i] * numer / denom) : 0.0
        end
        return λ
    end

    λ = compute_lambda(params, state)

    # analytic expected values for identity M: λ_i = β * c_i * I_i / N_i
    expected = [β * act_levels[i] * (Ivec[i] / class_sizes[i]) for i in 1:n]

    @test all(isapprox.(λ, expected; atol=1e-12))
end
