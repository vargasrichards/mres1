using Test
using MixSwitchStochEpi
using LinearAlgebra

@testset "NGM and R0 calculations" begin
    params, init, tspan = MixSwitchStochEpi.default_parametrisation()
    K = MixSwitchStochEpi.Kmat(params)
    @test size(K) == (params.n_activity, params.n_activity)

    # compute_R0 should equal spectral_radius(K)
    r0 = MixSwitchStochEpi.compute_R0(params)
    @test isapprox(r0, maximum(real(eigvals(K))))

    @testset "calibrate to a target R0 and check" begin
    target = 1.5
    params_cal = MixSwitchStochEpi.calibrate_parms(params, target)
    @test isapprox(MixSwitchStochEpi.compute_R0(params_cal), target; atol=1e-6)
    end
end

@testset "R0 determines whether an epidemic occurs" begin
    params, init, tspan = MixSwitchStochEpi.default_parametrisation()
    # For stochastic sims: when R0 < 1 extinction probability should be ~1
    params_low = MixSwitchStochEpi.calibrate_parms(params, 0.8) # subcritical (R0=0.8)
    p_ext_low = MixSwitchStochEpi.pr_fadeout(params_low, init; num_initial = 1, n_sims = 2000, t_cutoff = 50.0, threshold_fraction = 0.05)
    @test p_ext_low >= 0.95

    # When R0 > 1 extinction probability should be substantially less than 1
    params_high = MixSwitchStochEpi.calibrate_parms(params, 2.0) # supercritical (R0=2)
    p_ext_high = MixSwitchStochEpi.pr_fadeout(params_high, init; num_initial = 1, n_sims = 2000, t_cutoff = 50.0, threshold_fraction = 0.05)
    @test p_ext_high < 0.9
end

@testset "" begin
    
end