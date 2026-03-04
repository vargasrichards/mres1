# functions to help in calculating the next-generation matrix
using LinearAlgebra

"""
    ρMat(parms::SEIRSStochParams)

Get the matrix of ρ(i,j) from supplied `parms`
"""
function ρMat(parms::SEIRSStochParams)
    return parms.M
end


"""
    Tmat(parms::SEIRSStochParams)

Get the transmission matrix from supplied `parms` 
"""
function Tmat(parms::SEIRSStochParams)
    n = length(parms.act_levels)
    return [Tmat_entry(parms, i, j) for i = 1:n, j = 1:n]
end


""""
    Tmat_entry(parms::SEIRSStochParams, i::Int, j::Int)

Get entry at position `i`,`j` of transmission matrix for `parms`
"""
function Tmat_entry(parms::SEIRSStochParams, i::Int, j::Int)
    β = parms.β
    M = parms.M
    act_levels = parms.act_levels
    entry = β * act_levels[j] * M[j,i]
    return entry
end

"""
    Sigma_mat(parms::SEIRSStochParams)

Make the transition matrix from `parms`

This is used in the next-generation matrix calculation.
"""
function Sigma_mat(parms::SEIRSStochParams)
    n_activity = parms.n_activity
    Σ = parms.W - I(n_activity) * parms.γ
    return Σ
end


"""
    Kmat(parms::SEIRSStochParams)

Make the Next-generation matrix K of `parms`
"""
function Kmat(parms::SEIRSStochParams)
    Σmat = Sigma_mat(parms)
    Kmat = -Tmat(parms) * inv(Σmat)
    return Kmat
end

"""
    compute_R0(parms::SEIRSStochParams)

Compute R0 as `spectral_radius` of next-generation matrix K of `parms`
"""
function compute_R0(parms::SEIRSStochParams)
    K = Kmat(parms)
    R0 = spectral_radius(K)
    return R0
end

"""
    calibrate_parms(draft_parms::SEIRSStochParams, target_r0::Float64)

Calibrate the `draft_parms` to produce the desired `target_r0`
"""
function calibrate_parms(draft_parms::SEIRSStochParams, target_r0::Float64)
    draft_r0 = compute_R0(draft_parms)
    β_new = draft_parms.β * target_r0 / draft_r0
    parms = SEIRSStochParams(
        β_new,
        draft_parms.σ,
        draft_parms.γ,
        draft_parms.ω,
        draft_parms.M,
        draft_parms.W,
        draft_parms.n_activity,
        draft_parms.class_sizes,
        draft_parms.act_levels,
        draft_parms.pop_size
    )
    @assert compute_R0(parms) ≈ target_r0 "Calibration failed: new parameters do not give target_r0"
    return parms
end