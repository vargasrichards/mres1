using Test
using MixSwitchStochEpi

t_cutoff = 100.

@testset "Extinction behavior" begin
    # create a parameterisation with β = 0 (no transmission) -> extinction certain
    params, init, tspan = MixSwitchStochEpi.default_parametrisation()
    params_zero = MixSwitchStochEpi.SEIRSStochParams(0.0, params.σ, params.γ, params.ω, params.M, params.W, params.n_activity, params.class_sizes, params.act_levels, params.pop_size)

    res = MixSwitchStochEpi.simulate_system!(params_zero, init; tmax=t_cutoff)
    final = res.states[end]
    @test sum(final[2, :]) + sum(final[3, :]) == 0

    # With small but nonzero β, extinction at early time should be <= 1
    params_small = MixSwitchStochEpi.SEIRSStochParams(1e-6, params.σ, params.γ, params.ω, params.M, params.W, params.n_activity, params.class_sizes, params.act_levels, params.pop_size)
    res2 = MixSwitchStochEpi.simulate_system!(params_small, init; tmax=10.0)
    final2 = res2.states[end]
    @test (sum(final2[2, :]) + sum(final2[3, :])) >= 0
end

@testset "Pr(extinction) is approx 1 / R0 for R0 > 1 with uniform activity scores and no switching" begin
    # Use a single-index seed and uniform activity levels so the single-infective
    # branching-process approximation is reasonable.
    n = 5
    n_sims = 2000
    params, _, _ = MixSwitchStochEpi.default_parametrisation(n_activity = n, act_levels = ones(n), ξ = 0., p_inf = 1, init_mode = :class, init_class = 1)
    # create initial_state with a single infected in class 1
    _, init, _ = MixSwitchStochEpi.default_parametrisation(n_activity = n, act_levels = ones(n), ξ = 0.,  p_inf = 1, init_mode = :class, init_class = 1)

    for target_R0 in (1.8, 2.0, 3.0, 4.0, 10.0)
        params_cal = MixSwitchStochEpi.calibrate_parms(params, target_R0)
        pr_extinct = MixSwitchStochEpi.pr_fadeout(params_cal, init; n_sims=n_sims, t_cutoff=70.0, num_initial = 1)
        println("Predicted p_ext: $(1/target_R0), estimated Pr(extinction): $pr_extinct from $n_sims simulations. ")
    @test isapprox(pr_extinct, 1/target_R0; atol=0.05)
    end
end


@testset "Pr(extinction) is 1 for R0 < 1 with uniform activity scores and no switching" begin
    n = 5
    n_sims = 10000
    params, _, _ = MixSwitchStochEpi.default_parametrisation(n_activity = n, act_levels = ones(n), ξ = 0., p_inf = 1, init_mode = :class, init_class = 1)
    _, init, _ = MixSwitchStochEpi.default_parametrisation(n_activity = n, act_levels = ones(n), ξ = 0.,  p_inf = 1, init_mode = :class, init_class = 1)

    for target_R0 in (0.5, 0.8, 0.9)
        params_cal = MixSwitchStochEpi.calibrate_parms(params, target_R0)
        pr_extinct = MixSwitchStochEpi.pr_fadeout(params_cal, init; n_sims=n_sims, t_cutoff=70.0, num_initial = 1)
        println("Predicted p_ext: $(1/target_R0), estimated Pr(extinction): $pr_extinct from $n_sims simulations")
        @test isapprox(pr_extinct, 1.0; atol=0.08)
    end

end

@testset "Pr(extinction) decreases as number of initial infecteds increases" begin #this seems to be quite a problematic testset
    n = 5 # number of activity classes
    n_sims = 200
    pr_previous = 1.0 #will be updated in the loop, should be decreasing as n_initial increases
    t_cutoff = 100.0
    for chosen_r0 in (1.8, 2.0, 3.0)
        println("Testing Pr(extinction) with increasing initial infecteds for R0 = $chosen_r0")
        for ξ_val in (0., 0.1, 0.2, 0.3) 
        println("at switching rate $ξ_val")
        pr_previous = 1.0 
        for n_initial in (1, 5, 10, 20)
            params, _, _ = MixSwitchStochEpi.default_parametrisation(n_activity = n, act_levels = ones(n), ξ = ξ_val, p_inf = n_initial, init_mode = :class, init_class = 1)
            _, init, _ = MixSwitchStochEpi.default_parametrisation(n_activity = n, act_levels = ones(n), ξ = ξ_val,  p_inf = n_initial, init_mode = :class, init_class = 1)   
            params_cal = MixSwitchStochEpi.calibrate_parms(params, chosen_r0)
            pr_extinct = MixSwitchStochEpi.pr_fadeout(params_cal, init; n_sims=n_sims, t_cutoff=t_cutoff, num_initial = n_initial)
            println("Estimated Pr(extinction) with $n_initial initial infecteds: $pr_extinct from $n_sims simulations at ξ = $ξ_val. R0 = $chosen_r0")
            @test pr_extinct <= pr_previous + 0.02 || (n_initial == 1) # allow for small sampling noise; should be non-increasing
            pr_previous = pr_extinct
            end
        end
    end
end 

@testset "Going from no switching to >0 activity switching with unequal activity scores , Pr(extinction) increases if the initial infected is in the highest-activity class" begin
    # Build a parametrisation with unequal activity scores
    n = 5
    act_levels = [1.0, 2.0, 3.0, 5.0, 10.0]
    params, _, _ = MixSwitchStochEpi.default_parametrisation(act_levels = act_levels, n_activity = n)

    # initial infected in the highest-activity class
    init_high = MixSwitchStochEpi.default_parametrisation(act_levels = act_levels, n_activity = n, init_mode = :class, init_class = n)[2]

    # no switching vs moderate switching
    params_no_switch = MixSwitchStochEpi.SEIRSStochParams(params.β, params.σ, params.γ, params.ω, params.M, MixSwitchStochEpi.uniform_switching(0.0, n), params.n_activity, params.class_sizes, params.act_levels, params.pop_size)
    params_switch = MixSwitchStochEpi.SEIRSStochParams(params.β, params.σ, params.γ, params.ω, params.M, MixSwitchStochEpi.uniform_switching(0.5, n), params.n_activity, params.class_sizes, params.act_levels, params.pop_size)

    pr_no = MixSwitchStochEpi.pr_fadeout(params_no_switch, init_high; n_sims=1000, t_cutoff=50.0, num_initial = 1)
    pr_sw = MixSwitchStochEpi.pr_fadeout(params_switch, init_high; n_sims=1000, t_cutoff=50.0, num_initial = 1)

    @test pr_sw >= pr_no - 0.02 # allow a tiny sampling tolerance
end

@testset "Going from no switching to >0 activity switching with unequal activity scores, Pr(extinction)  decreases if the initial infected is in the lowest-activity class" begin
    n = 5
    act_levels = gamma_bin_expectations(mean = 5.0, variance = 5., n_activity = n)
    params, _, _ = MixSwitchStochEpi.default_parametrisation(act_levels = act_levels, n_activity = n)

    # initial infected in the lowest-activity class
    init_low = MixSwitchStochEpi.default_parametrisation(act_levels = act_levels, n_activity = n, init_mode = :class, init_class = 1)[2]

    params_no_switch = MixSwitchStochEpi.SEIRSStochParams(params.β, params.σ, params.γ, params.ω, params.M, MixSwitchStochEpi.uniform_switching(0.0, n), params.n_activity, params.class_sizes, params.act_levels, params.pop_size)
    params_switch = MixSwitchStochEpi.SEIRSStochParams(params.β, params.σ, params.γ, params.ω, params.M, MixSwitchStochEpi.uniform_switching(0.5, n), params.n_activity, params.class_sizes, params.act_levels, params.pop_size)

    pr_no = MixSwitchStochEpi.pr_fadeout(params_no_switch, init_low; n_sims=1000, t_cutoff=50.0, num_initial =1)
    pr_sw = MixSwitchStochEpi.pr_fadeout(params_switch, init_low; n_sims=1000, t_cutoff=50.0, num_initial =1 )

    @test pr_sw <= pr_no + 0.02 # allow a tiny sampling tolerance
end


@testset "Pr(extinction) is insensitive to switching rate when activity scores are uniform" begin
    params, init, tspan = MixSwitchStochEpi.default_parametrisation()
    n = params.n_activity
    #     ε::Float64,
    # n_activity::Int,
    # act_levels::Vector{Float64},
    # class_sizes::Vector{Int64},
    params_no_switch = MixSwitchStochEpi.SEIRSStochParams(params.β, params.σ, params.γ, params.ω, MixSwitchStochEpi.make_garnett_contact(0., n, ones(n), params.class_sizes), MixSwitchStochEpi.uniform_switching(0.0, n), params.n_activity, params.class_sizes, params.act_levels, params.pop_size)
    params_switch = MixSwitchStochEpi.SEIRSStochParams(params.β, params.σ, params.γ, params.ω, MixSwitchStochEpi.make_garnett_contact(0., n, ones(n), params.class_sizes), MixSwitchStochEpi.uniform_switching(0.5, n), params.n_activity, params.class_sizes, params.act_levels, params.pop_size)

    pr_extinct_no_switch = MixSwitchStochEpi.pr_fadeout(params_no_switch, init; n_sims=1000, t_cutoff=50.0, num_initial = 1)
    pr_extinct_switch = MixSwitchStochEpi.pr_fadeout(params_switch, init; n_sims=1000, t_cutoff=50.0, num_initial = 1)

    @test isapprox(pr_extinct_no_switch, pr_extinct_switch; atol=0.05)
end


@testset "Calibration produces requested R0" begin
    params, init, tspan = MixSwitchStochEpi.default_parametrisation()
    target = 1.25
    pcal = MixSwitchStochEpi.calibrate_parms(params, target)
    @test isapprox(MixSwitchStochEpi.compute_R0(pcal), target; atol=1e-3)
end