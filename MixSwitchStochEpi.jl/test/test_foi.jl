# force of infection tests
# Alexis Vargas Richards, Jan 2026

using Test



@testset "FOI total can be decomposed into component FOIs from each class" begin
    # we can calculate the FOI from each class separately and check that they sum to the total FOI.
    params, init, tspan = MixSwitchStochEpi.default_parametrisation()
        C = params.M
        act_levels = params.act_levels
        class_sizes = params.class_sizes
        n = params.n_activity
        # calculate FOI from each class separately
        foi_from_class = zeros(n)
        for j in 1:n
            foi_from_class[j] = sum(C[:, j] .* act_levels[j] .* init[3, j]) / sum(class_sizes .* act_levels)
        end
        total_foi = sum(foi_from_class)
        # calculate total FOI directly
        total_foi_direct = sum(C .* (act_levels .* init[3, :])') / sum(class_sizes .* act_levels)
        @test isapprox(total_foi, total_foi_direct; atol=1e-6)
end

# @testset "FOI "