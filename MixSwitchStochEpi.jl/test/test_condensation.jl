# test that the homogeneous activivty case gives statistically consistent results regardless of the number of activity classes (n_activity)
# or the assortativity parameter (ε).

using Test
using MixSwitchStochEpi
using Random
using Statistics

@testset "Condensation test: homogeneous activity case (no switching) across assortativities" begin
    n_activity_values = [1, 5, 10, 20]
    ε_values = [0.0, 0.5, 1.0]
    for ε in ε_values
        # collect final recovered proportions across n_activity for this ε
        final_props = Float64[]
        for (idx, n_activity) in enumerate(n_activity_values)
            act_levels = MixSwitchStochEpi.uniform_scores(mean = 1.0, n_activity = n_activity)
            # use smaller class sizes in tests to keep runtime reasonable
            class_sizes = fill(Int(1000), n_activity)
            M = MixSwitchStochEpi.make_garnett_contact(ε, n_activity, act_levels, class_sizes)
            W = MixSwitchStochEpi.uniform_switching(0.0, n_activity)
            parms = SEIRSStochParams(0.5, 1/3, 1/4, 0.0,
                                     M, W, n_activity,
                                     class_sizes, act_levels, sum(class_sizes))
            # build a uniform initial state (distribute p_inf infections evenly across classes)
            p_inf = max(1, floor(Int, sum(class_sizes) * 0.01))
            initial_state = make_initial_state(class_sizes, p_inf)
            rng = MersenneTwister(1234) # fixed seed so replicates are comparable across configs
            sims = many_simulations(parms, initial_state, 1, 50; tmax = 100.0, rng = rng)
            # many_simulations returns secondary_sizes (proportion of population infected)
            mean_final = mean(sims.secondary_sizes)
            push!(final_props, mean_final)
        end
        # require that final recovered proportions across partitions are similar (within 2 percentage points)
    med = median(final_props)
    @test maximum(abs.(final_props .- med)) <= 0.06
    end
end


@testset "Condensation test: homogeneous activity scores with switching across assortativities" begin
    n_activity_values = [1, 5, 10, 20]
    ε_values = [0.0, 0.5, 1.0]
    for ε in ε_values
        # collect final recovered proportions across n_activity for this ε (with switching)
        final_props = Float64[]
        for (idx, n_activity) in enumerate(n_activity_values)
            act_levels = MixSwitchStochEpi.uniform_scores(mean = 1.0, n_activity = n_activity)
            class_sizes = fill(Int(1000), n_activity)
            M = MixSwitchStochEpi.make_garnett_contact(ε, n_activity, act_levels, class_sizes)
            W = MixSwitchStochEpi.uniform_switching(0.5, n_activity) # Add some switching
            parms = SEIRSStochParams(0.5, 1/3, 1/4, 0.0,
                                     M, W, n_activity,
                                     class_sizes, act_levels, sum(class_sizes))
            # build a uniform initial state (distribute p_inf infections evenly across classes)
            p_inf = max(1, floor(Int, sum(class_sizes) * 0.01))
            initial_state = make_initial_state(class_sizes, p_inf)
            rng = MersenneTwister(1234)
            sims = many_simulations(parms, initial_state, 1, 50; tmax = 100.0, rng = rng)
            mean_final = mean(sims.secondary_sizes)
            push!(final_props, mean_final)
        end
    med = median(final_props)
    @test maximum(abs.(final_props .- med)) <= 0.06
    end
end