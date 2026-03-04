
using QuadGK

"""
    gamma_bin_expectations(; mean::Union{Nothing,Float64}=nothing,
                             variance::Union{Nothing,Float64}=nothing,
                             shape::Union{Nothing,Float64}=nothing,
                             rate::Union{Nothing,Float64}=nothing,
                             n_activity::Int)

Find the activity scores for each activity class using a discretised Gamma distribution.

You can provide either (mean, variance) or (shape, rate) to parameterise the Gamma:

- If `mean` and `variance` are given (backwards-compatible), the Gamma shape-rate
  parameters are computed as
    (alpha = mean^2 / variance, rate = mean / variance)

- Alternatively, you can pass `shape` and `rate` directly. The function uses the
  Gamma(α, θ) distribution where θ is the scale parameter (θ = 1 / rate).

The function returns a vector of length `n_activity` with the conditional
expected activity score in each discretised bin.

# Examples
```jldoctest
julia> using MixSwitchStochEpi

julia> gamma_bin_expectations(mean = 5., variance = 5., n_activity = 5)
5-element Vector{Float64}:
 2.330329878881237
 3.6287950366507467
 4.677938202312152
 5.919141132242913
 8.443795749913242

# or using shape & rate directly (equivalent to mean=5,variance=5 here)
julia> gamma_bin_expectations(shape = 25/5, rate = 5/5, n_activity = 5)
5-element Vector{Float64}:
 2.330329878881237
 3.6287950366507467
 4.677938202312152
 5.919141132242913
 8.443795749913242
```
"""
function gamma_bin_expectations(; mean::Union{Nothing,Float64}=nothing,
                                variance::Union{Nothing,Float64}=nothing,
                                shape::Union{Nothing,Float64}=nothing,
                                rate::Union{Nothing,Float64}=nothing,
                                n_activity::Int)
    # Determine shape (α) and rate (β) from provided keywords.
    if shape !== nothing && rate !== nothing
        α = shape
        β = rate
    elseif mean !== nothing && variance !== nothing
        α = mean^2 / variance
        β = mean / variance
    else
        error("gamma_bin_expectations: provide either (mean and variance) or (shape and rate)")
    end

    # Construct Gamma with scale θ = 1 / rate
    γ = Gamma(α, 1 / β)
    edges = quantile.(Ref(γ), range(0, stop = 1, length = n_activity + 1))
    expectations = Float64[]
    for i = 1:n_activity
        a, b = edges[i], edges[i + 1]
        score, _ = quadgk(x -> x * pdf(γ, x), a, b)
        conditional_score = n_activity * score
        push!(expectations, conditional_score)
    end
    return expectations
end

"""
    uniform_scores(; mean::Float64, n_activity::Int)

Generate uniform activity scores across activity classes. 

This is then a homogeneous population.
Useful when testing the model.

# Examples
```jldoctest
julia> using MixSwitchStochEpi

julia> uniform_scores(mean = 5.,n_activity = 10)
10-element Vector{Float64}:
 5.0
 5.0
 5.0
 5.0
 5.0
 5.0
 5.0
 5.0
 5.0
 5.0
```
"""
function uniform_scores(; mean::Float64, n_activity::Int)
    return (fill(mean, n_activity))
end
