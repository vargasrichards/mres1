using Test
using MixSwitchStochEpi

@testset "initial conditions" begin
    params, init, tspan = MixSwitchStochEpi.default_parametrisation()
    @test size(init) == (4, params.n_activity)

    # Test single-class seeding
    params2, init2, _ = MixSwitchStochEpi.default_parametrisation(init_mode=:class, init_class=2)
    @test sum(init2[2, :]) >= 1
    @test init2[2, 1] == 0 && init2[2, 2] >= 1
end