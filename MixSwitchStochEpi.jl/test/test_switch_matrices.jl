# tests implemented for the switching matrices
using Test
using LinearAlgebra

@testset "Ensure switching matrix generates equal-sized stationary dist " begin
    @testset "stationary dist tests" begin
        for n in 1:20
            W = uniform_switching(0.1, n)
            stat_Wunif = find_steadystate(W)
            @test allequal(isapprox.(stat_Wunif, first(stat_Wunif), atol=1e-8)) == true 

            Wadj = constrained_switching(0.1, n)

            stat_Wadj = find_steadystate(Wadj)
            @test allequal(isapprox.(stat_Wadj, first(stat_Wadj), atol=1e-8)) == true 
        end
    end
end


@testset "Ensure switching matrix sums to 0" begin 
    @testset "Adj only switching matrix sums to 0 " begin 
        Wunif = constrained_switching(0.1, 10)
        @test validate_generator(Wunif) == true
        
    end 

    @testset "uniformly switching switching matrix sums to 0 " begin
        S = uniform_switching(0.1, 10)
        @test validate_generator(S) == true
    end
end


@testset "W transpose = W (switching matrix)" begin # note that this isn't necessarily a rigorous test
    for n in 1:20
        W = uniform_switching(0.1,n)
        @test transpose(W) == W
        Wadj = constrained_switching(0.1, n)
        @test transpose(Wadj) == Wadj
    end
end