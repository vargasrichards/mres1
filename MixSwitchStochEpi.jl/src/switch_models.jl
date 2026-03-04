# this file contains all the switching models currently implemented.

using LinearAlgebra

"""
    constrained_switching(ξ::Float64, n_activity::Int)

Generate an adjacent-constrained switching matrix.  

The expected residence time in the most- or least-active classes
is half that of the intermediate classes. 


# Examples
```jldoctest
julia> using MixSwitchStochEpi

julia> constrained_switching(0.1, 10)
10×10 Matrix{Float64}:
 -0.05   0.05   0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0
  0.05  -0.1    0.05   0.0    0.0    0.0    0.0    0.0    0.0    0.0
  0.0    0.05  -0.1    0.05   0.0    0.0    0.0    0.0    0.0    0.0
  0.0    0.0    0.05  -0.1    0.05   0.0    0.0    0.0    0.0    0.0
  0.0    0.0    0.0    0.05  -0.1    0.05   0.0    0.0    0.0    0.0
  0.0    0.0    0.0    0.0    0.05  -0.1    0.05   0.0    0.0    0.0
  0.0    0.0    0.0    0.0    0.0    0.05  -0.1    0.05   0.0    0.0
  0.0    0.0    0.0    0.0    0.0    0.0    0.05  -0.1    0.05   0.0
  0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.05  -0.1    0.05
  0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.05  -0.05
```
"""
function constrained_switching(ξ::Float64, n_activity::Int)
    main_diag = fill(-ξ, n_activity)
    main_diag[1] = -ξ/2
    main_diag[end] = -ξ/2
    off_diag = fill(ξ / 2, n_activity - 1)
    return convert(Matrix,Tridiagonal(off_diag, main_diag, off_diag))
end



"""
    uniform_switching(ξ::Float64, n_activity::Int)

Generate an uniform switching matrix.  

The expected residence time in a given activity class
is 1/ξ where xi is the switching rate

# Examples
```jldoctest
julia> using MixSwitchStochEpi

julia> uniform_switching(0.1,10)
10×10 Matrix{Float64}:
 -0.1         0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111
  0.0111111  -0.1         0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111
  0.0111111   0.0111111  -0.1         0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111
  0.0111111   0.0111111   0.0111111  -0.1         0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111
  0.0111111   0.0111111   0.0111111   0.0111111  -0.1         0.0111111   0.0111111   0.0111111   0.0111111   0.0111111
  0.0111111   0.0111111   0.0111111   0.0111111   0.0111111  -0.1         0.0111111   0.0111111   0.0111111   0.0111111
  0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111  -0.1         0.0111111   0.0111111   0.0111111
  0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111  -0.1         0.0111111   0.0111111
  0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111  -0.1         0.0111111
  0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111   0.0111111  -0.1
```
"""
function uniform_switching(ξ::Float64, n_activity::Int)
    W = fill(ξ / (n_activity-1), n_activity, n_activity)
    W[diagind(W)] .= -ξ
    W
end