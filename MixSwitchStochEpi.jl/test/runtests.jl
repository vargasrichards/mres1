using Test
using Compat
using TestSetExtensions
using MixSwitchStochEpi
using LinearAlgebra
using TidierData

@testset ExtendedTestSet "MixSwitchStochEpi.jl" begin
    include("test_contact_matrices.jl")
    include("test_switch_matrices.jl")
    include("test_gamma_activity.jl")   
    include("test_initial_conditions.jl")
    include("test_ngm.jl")
    include("test_infection_flow.jl") 
    include("test_extinction.jl")
    include("test_foi.jl")
    # include("test_gpu.jl")
    include("test_condensation.jl")
    include("test_aqua.jl")
end