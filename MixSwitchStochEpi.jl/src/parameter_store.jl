# first define the stochastic parameter struct


"""
The ScalarParms struct defines the parameters needed for the homogeneous reference
    model. 
"""
struct ScalarParms
    β::Float64
    σ::Float64
    γ::Float64
    ω::Float64
end

"""
SEIR parameters from (Britton et al., 2020) Science
approximately for SARS-CoV-2.

Note that β needs to be respecified; the value
of 0.5 is just an arbitrary value.
"""
function get_scalar_parms()
    scalar_parms = ScalarParms(
        0.5,      # β
        1/3,      # σ
        1/4,      # γ
        0.0       # ω 
    )
    return scalar_parms
end


