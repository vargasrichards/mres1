# stochastic SEIR(S) model with class switching
#
#
# A. Vargas Richards Nov. 2025

using Random
using Distributions
using LinearAlgebra
using Statistics
using DataFrames

"""
    SEIRSStochParams()

Construct the stochastic parameter struct for simulation
"""
struct SEIRSStochParams
    β::Float64          # transmission 
    σ::Float64          # progression to infectivity rate E->I
    γ::Float64          # recovery rate I->R
    ω::Float64          # immune waning rate R->S (usually set = 0 in our proj.)

    M::Matrix{Float64}  # mixing/contact matrix used for FOI
    W::Matrix{Float64}  # switching/movement rate matrix between activity classes (CTMC generator)
    n_activity::Int # number of activity classes (usually we set this to around 5-10 or so.)
    class_sizes::Vector{Int64}   # initial activity class sizes
    act_levels::Vector{Float64}  # relative activity level per class
    pop_size::Int 
end

"""
    simulate_system!(params::SEIRSStochParams, initial_state::Matrix{Int}, tmat = Inf, rng = Random.GLOBAL_RNG)

Simulate a SEIR(S) pathogen on a time-varying activity-structured population.

Uses Gillespie's Direct simulation algorithm (exact)

"""
function simulate_system!(
    params::SEIRSStochParams,
    initial_state::Matrix{Int};
    tmax::Float64 = Inf,
    rng::AbstractRNG = Random.TaskLocalRNG(),
    store_states::Bool = true,
)
    n = params.n_activity
    @assert size(initial_state) == (4, n)

    # integer compartment vectors (mutable)
    S = collect(initial_state[1, :])
    E = collect(initial_state[2, :])
    I = collect(initial_state[3, :])
    R = collect(initial_state[4, :])

    comps = (S, E, I, R)

    times = Float64[0.0]
    states = store_states ? [copy(initial_state)] : Matrix{Int}[]
    t = 0.0

    # preallocate reusable buffers to avoid repeated allocations
    N = Vector{Float64}(undef, n)
    λ = zeros(Float64, n)

    # local aliases for faster access
    local_M = params.M
    local_W = params.W
    local_act = params.act_levels
    local_β = params.β
    local_σ = params.σ
    local_γ = params.γ
    local_ω = params.ω

    # main loop
    while true
        # compute class totals once
        ssum = 0
        esum = 0
        isum = 0
        rsum = 0
        @inbounds for i in 1:n
            ssum += S[i]
            esum += E[i]
            isum += I[i]
            rsum += R[i]
            N[i] = float(S[i] + E[i] + I[i] + R[i])
        end
        # termination check
        (esum + isum) == 0 && break

        # FOI with dynamic class sizes
        @inbounds for i in 1:n
            denom = 0.0
            numer = 0.0
            for j in 1:n
                # skip self if desired; matrix may have zeros on diagonal
                w = local_M[i, j] * N[j]
                denom += w
                if N[j] > 0.0
                    numer += w * (I[j] / N[j])
                end
            end
            λ[i] = denom > 0.0 ? (local_β * local_act[i] * numer / denom) : 0.0
        end

        # transition rates
        rate_inf = 0.0
        @inbounds for i in 1:n
            rate_inf += λ[i] * S[i] #computing FOI * susceptibles for each class and summing to get total infection rate across all classes
        end
        rate_EI = local_σ * esum
        rate_IR = local_γ * isum
        rate_RS = local_ω * rsum

        rate_move = 0.0
        @inbounds for k in 1:4
            comp = comps[k]
            for i in 1:n, j in 1:n
                if i != j
                    rate_move += local_W[i, j] * comp[i]
                end
            end
        end

        total_rate = rate_inf + rate_EI + rate_IR + rate_RS + rate_move
        total_rate <= 0 && break

        # sample next event time (faster exponential sampling)
        t_new = t - log(rand(rng)) / total_rate
        if t_new > tmax
            t = tmax
            if store_states
                push!(times, t)
                push!(states, vcat(S', E', I', R'))
            end
            break
        end
        t = t_new

        r = rand(rng) * total_rate

        # event selection
        if r < rate_inf
            thresh = rand(rng) * rate_inf
            acc = 0.0
            @inbounds for i in 1:n
                acc += λ[i] * S[i]
                if acc >= thresh
                    if S[i] > 0
                        S[i] -= 1
                        E[i] += 1
                    end
                    break
                end
            end

        elseif r < rate_inf + rate_EI
            thresh = rand(rng) * rate_EI
            acc = 0.0
            @inbounds for i in 1:n
                acc += local_σ * E[i]
                if acc >= thresh
                    E[i] -= 1
                    I[i] += 1
                    break
                end
            end

        elseif r < rate_inf + rate_EI + rate_IR
            thresh = rand(rng) * rate_IR
            acc = 0.0
            @inbounds for i in 1:n
                acc += local_γ * I[i]
                if acc >= thresh
                    I[i] -= 1
                    R[i] += 1
                    break
                end
            end

        elseif r < rate_inf + rate_EI + rate_IR + rate_RS
            thresh = rand(rng) * rate_RS
            acc = 0.0
            @inbounds for i in 1:n
                acc += local_ω * R[i]
                if acc >= thresh
                    R[i] -= 1
                    S[i] += 1
                    break
                end
            end

        else
            thresh = rand(rng) * rate_move
            acc = 0.0
            for k in 1:4
                comp = comps[k]
                @inbounds for i in 1:n
                    for j in 1:n
                        if i != j
                            acc += local_W[i, j] * comp[i]
                            if acc >= thresh
                                comp[i] -= 1
                                comp[j] += 1
                                break
                            end
                        end
                    end
                    if acc >= thresh
                        break
                    end
                end
                if acc >= thresh
                    break
                end
            end
        end

        if store_states
            push!(times, t)
            push!(states, vcat(S', E', I', R'))
        end
    end

    if !store_states
        push!(times, t)
        push!(states, vcat(S', E', I', R'))
    end

    return (times = times, states = states)
end


"""
    many_simulations(params::SEIRSStochParams, initial_state::Matrix{Int}, n_sims::Int; tmax::Float64 = Inf)

Run many simulations of the SEIR(S) model and return the final sizes (proportion of population infected) for each simulation.
Also returns the last time of the simulation. Our termination condition is that there are no infecteds left, so 
Corrects for the nunber of initial cases to return the secondary attack rate (final size of epidemic excluding initial cases).

# Arguments
- `params::SEIRSStochParams`
- `initial_state::Matrix{Int}`: the initial state of the population (4 x
n_activity matrix)
- `num_initial::Int`: the number of initial cases (used to compute secondary attack rate
- `n_sims::Int`: the number of simulations to run
- `tmax::Float64`: the maximum time to run each simulation for (default is Inf, meaning run until extinction)

"""
function many_simulations(
    params::SEIRSStochParams,
    initial_state::Matrix{Int64}, num_initial::Int,
    n_sims::Int;
    tmax::Float64 = Inf,
    rng::AbstractRNG = Random.TaskLocalRNG(),
)
    final_sizes = Float64[]
    end_times = Float64[] #will store the 
    for sim = 1:n_sims
        result = simulate_system!(params, initial_state, tmax = tmax, rng = rng, store_states = false)
        final_time = maximum(result.times)
        final_R = sum(result.states[end][4, :])
        # secondary R is :
        final_R -= num_initial #if we have a couple of infecteds this quantity shouldnt matter much
        push!(final_sizes, final_R / params.pop_size)
        push!(end_times, final_time)
    end
    return (secondary_sizes = final_sizes, durations = end_times)
end



"""
    pr_fadeout(params::SEIRSStochParams, 
                initial_state::Matrix{Int64},
                num_initial::Int,
                n_sims::Int,
                t_cutoff::Float64)

Compute probability of stochastic fadeout

# Arguments




# """
# function pr_fadeout(
#     params::SEIRSStochParams,
#     initial_state::Matrix{Int64};
#     num_initial::Int,
#     n_sims::Int,
#     t_cutoff::Float64 = 30.)

#     results = many_simulations(params, initial_state, num_initial, n_sims) #let tmax be INF here.
#     # has fields: results.secondary_sizes, results.durations
#     #make a histogram of final sizes to check for bimodality etc if desired
#     #println("final sizes are $(results.secondary_sizes) with mean $(mean(results.secondary_sizes))\n durations are $(results.durations) with mean duration $(results.durations)")

#     pr = sum(x -> x < t_cutoff, results.durations) / n_sims
#     # if hist == true
#     #     histogram(results, bins=30,
#     #     xlabel="Final size (secondary cases/population)", 
#     #     ylabel="Freq",  
#     #     title="f_s for $n_sims runs, R0 = $(compute_R0(params)), threshold = $threshold_value");
#     #     savefig("final_size_thresh$threshold_value,pr=$pr.svg")
#     # end
#     return pr
# end

function pr_fadeout(
    params::SEIRSStochParams,
    initial_state::Matrix{Int64};
    num_initial::Int = 1,
    n_sims::Int = 1000,
    t_cutoff::Float64,
    tmax::Union{Nothing, Float64} = nothing,
    threshold_fraction::Float64 = 0.05,
    threshold_value::Union{Nothing, Float64} = nothing,
    rng::AbstractRNG = Random.TaskLocalRNG(),
)
    # Allow callers to pass either `tmax` or `t_cutoff`, and either `threshold_value`
    # (legacy name) or `threshold_fraction` (new name).
    t_use = tmax === nothing ? t_cutoff : tmax
    thr = threshold_value === nothing ? threshold_fraction : threshold_value

    n_extinct = 0
    for sim = 1:n_sims
        result = simulate_system!(params, initial_state, tmax = t_use, rng = rng, store_states = false)

        final_state = result.states[end]
        final_E = sum(final_state[2, :])
        final_I = sum(final_state[3, :])
        final_R = sum(final_state[4, :])

        # True fadeout criteria:
        # 1) chain terminated by t_use (E+I == 0), and
        # 2) cumulative secondary cases (final_R - num_initial) is small relative to pop
        secondary_cases = final_R - num_initial
        if (final_E + final_I == 0) && (secondary_cases < thr * params.pop_size)
            n_extinct += 1
        end
    end
    return n_extinct / n_sims
end

function scan_assortativity(
    params::SEIRSStochParams,
    initial_state::Matrix{Int64},
    n_sims::Int,
    n_values::Int,
    threshold_value::Float64,
    tmax::Float64 = Inf,
)
    mean_final_sizes = Float64[]
    median_final_sizes = Float64[]
    pr_fadeouts = Float64[]
    ε_tested = Float64[]
    for ε in range(0, 1, step = 1 / n_values)
        println("Testing assortativity ε = $ε with $n_sims replicates")
        #remake the contact matrix here.
        curr_params = SEIRSStochParams(
            params.β,
            params.σ,
            params.γ,
            params.ω,
            make_contact_matrix(
                ε,
                params.n_activity,
                params.act_levels,
                params.class_sizes,
            ),
            params.W,
            params.n_activity,
            params.class_sizes,
            params.act_levels,
            params.pop_size,
        )
        final_sizes = many_simulations(curr_params, initial_state, n_sims, tmax = tmax)
        prob_fadeout = pr_fadeout(final_sizes, n_sims, threshold_value) #the probability of an early epidemic fadeout 
        push!(mean_final_sizes, mean(final_sizes))
        push!(median_final_sizes, median(final_sizes))
        #push additional summary statistics here??
        push!(pr_fadeouts, prob_fadeout)
        push!(ε_tested, ε)
    end

    return (DataFrame(
        assortativity = ε_tested,
        mean_final_size = mean_final_sizes,
        median_final_size = median_final_sizes,
        probab_fadeout = pr_fadeouts,
        n_replicates = n_sims,
    ))
end

function class_size_flux(results)
    times = results.times
    n_times = length(times)::Int

    n_classes = size(results.states[1], 2)

    class_sizes = zeros(n_times, n_classes)

    for (i, k) in enumerate(times)
        sysstate = results.states[i]
        class_sizes[i, :] = vec(sum(sysstate, dims = 1))
    end
    return times, class_sizes
end



function plot_class_sizes(result)
    times, class_sizes = class_size_flux(result)
    display(
        plot(times, class_sizes, xlabel = "Time", ylabel = "Class sizes", ylim = (0, Inf)),
    )
    savefig("class_size_singlerun.svg")
end

function scan_switching(
    params::SEIRSStochParams,
    initial_state::Matrix{Int64},
    n_sims::Int,
    n_values::Int,
    threshold_value::Float64,
    tmax::Float64 = Inf,
)
    mean_final_sizes = Float64[]
    median_final_sizes = Float64[]
    pr_fadeouts = Float64[]
    ξ_tested = Float64[]
    for ξ in range(0, 1, step = 1 / n_values)
        println("Testing switching ξ = $ξ with $n_sims replicates")
        W = simple_switching(ξ, params.n_activity)
        curr_params = SEIRSStochParams(
            params.β,
            params.σ,
            params.γ,
            params.ω,
            params.M,
            W,
            params.n_activity,
            params.class_sizes,
            params.act_levels,
            params.pop_size,
        )

        final_sizes = many_simulations(curr_params, initial_state, n_sims, tmax = Inf)
        prob_fadeout = pr_fadeout(final_sizes, n_sims, threshold_value) #the probability of an early epidemic fadeout 
        push!(mean_final_sizes, mean(final_sizes))
        push!(median_final_sizes, median(final_sizes))
        #push additional summary statistics here??
        push!(pr_fadeouts, prob_fadeout)
        push!(ξ_tested, ξ)
    end

    return (DataFrame(
        switching = ξ_tested,
        mean_final_size = mean_final_sizes,
        median_final_size = median_final_sizes,
        probab_fadeout = pr_fadeouts,
    ))
end

function ass_switch(
    params::SEIRSStochParams,
    initial_state::Matrix{Int64},
    n_sims::Int,
    n_values::Int,
    threshold_value::Float64,
    tmax::Float64 = Inf,
)
    mean_final_sizes = Float64[]
    median_final_sizes = Float64[]
    pr_fadeouts = Float64[]
    ξ_tested = Float64[]
    ε_tested = Float64[]
    num_points = length(range(0, 1, step = 1 / n_values))^2
    point_cnt = 0
    for ε in range(0, 1, step = 1 / n_values)
        M = make_contact_matrix(ε, params.n_activity, params.act_levels, params.class_sizes)
        for ξ in range(0, 1, step = 1 / n_values)
            point_cnt += 1
            println(
                "Testing switching ξ = $ξ, assortativity ε = $ε with $n_sims replicates",
            )
            println("at $point_cnt / $num_points")
            W = simple_switching(ξ, params.n_activity)
            curr_params = SEIRSStochParams(
                params.β,
                params.σ,
                params.γ,
                params.ω,
                M,
                W,
                params.n_activity,
                params.class_sizes,
                params.act_levels,
                params.pop_size,
            )

            final_sizes = many_simulations(curr_params, initial_state, n_sims, tmax = Inf)
            prob_fadeout = pr_fadeout(final_sizes, n_sims, threshold_value) #the probability of an early epidemic fadeout 
            push!(mean_final_sizes, mean(final_sizes))
            push!(median_final_sizes, median(final_sizes))
            #push additional summary statistics here??
            push!(pr_fadeouts, prob_fadeout)
            push!(ξ_tested, ξ)
            push!(ε_tested, ε)
        end
    end
    return (DataFrame(
        switching = ξ_tested,
        assortativity = ε_tested,
        mean_final_size = mean_final_sizes,
        median_final_size = median_final_sizes,
        probab_fadeout = pr_fadeouts,
    ))
end
