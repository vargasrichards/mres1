# tests for the DGM and activity score generation

using Test
using Statistics 


@testset "mean, variance of DGM bin scores -> mean, variance of Γ dist as N_act -> Inf" begin
    exact_mean = 10.
    exact_variance = 10.
   gamma_approx = gamma_bin_expectations(mean = exact_mean, variance = exact_variance, n_activity = 1000)
   @test isapprox(mean(gamma_approx), exact_mean, rtol = 1e-3)
   @test isapprox(var(gamma_approx), exact_variance, rtol = 1e-3)

    gamma_approx = gamma_bin_expectations(mean = exact_mean, variance = exact_variance, n_activity = 100000)
   @test isapprox(mean(gamma_approx), exact_mean, rtol = 1e-4)
   @test isapprox(var(gamma_approx), exact_variance, rtol = 1e-4)

end


@testset "appropriate validation of the mean and variance" begin
    
end