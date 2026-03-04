# test T0 functions and analysis

using Test
using MixSwitchStochEpi

@testset "extinction threshold calculation runs" begin
    params, init, tspan = MixSwitchStochEpi.default_parametrisation()
    # calibrate to R0=1.8 by default
    pcal = MixSwitchStochEpi.calibrate_parms(params, 1.8)
    # find threshold value for extinction probability of 0.05
    threshold = MixSwitchStochEpi.find_threshold(pcal; n_sims=1000, tmax=60.0, target_pext=0.05)
    @test threshold > 0
end

@testset "extinction threshold is higher (more cases required for P_ext < threshold) for lower-activity classes than higher-activity classes if a single class is infected. " begin
    


end


@testset "gap between extinction thresholds for lowest-activity and highest-activity classes is smaller when switching is elevated from 0" begin
    # Build a parametrisation with unequal activity scores
    n = 5
    act_levels = gamma_bin_expectations(mean = 5.0, variance = 5., n_activity = n)
    params, _, _ = MixSwitchStochEpi.default_parametrisation(act_levels = act_levels, n_activity = n)
    # calibrate to R0=1.8 by default
    pcal = MixSwitchStochEpi.calibrate_parms(params, 1.8)

    # find thresholds for no switching vs moderate switching
    threshold_low_activity = MixSwitchStochEpi.find_threshold(pcal; n_sims=1000, tmax=60.0, target_pext=0.05, ξ=0.0)




    threshold_moderate_switch = MixSwitchStochEpi.find_threshold(pcal; n_sims=1000, tmax=60.0, target_pext=0.05, ξ=0.5)

    println("Threshold with no switching: $threshold_no_switch")
    println("Threshold with moderate switching: $threshold_moderate_switch")

    @test threshold_moderate_switch < threshold_no_switch
end


