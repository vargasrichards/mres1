

"""
    make_initial_state(class_sizes::Vector{Int}, p_inf::Float64; init_mode::Symbol = :uniform,
    init_class::Int = 1)

Make the initial state of population to be simulated.

# Examples
```{jldoctest}
julia> using MixSwitchStochEpi
```

"""
function make_initial_state(
    class_sizes::Vector{Int},
    p_inf::Int;
    init_mode::Symbol = :uniform,
    init_class::Int = 1)

    n = length(class_sizes)
    E = zeros(Int, n)

    if p_inf < 0
        error("p_inf must be >= 0 (number of initially infected)")
    end

    if init_mode == :uniform
        # distribute p_inf as evenly as possible across classes into the I compartment
        if p_inf == 0
            # nothing infected
        else
            base = div(p_inf, n)
            rem = p_inf - base * n
            for i = 1:n
                E[i] = base
            end
            for i = 1:rem
                E[i] += 1
            end
        end
    elseif init_mode == :class
        E[init_class] = p_inf
    else
        return error("Unknown init_mode")
    end

    S = class_sizes .- E
    I = zeros(Int, n)
    R = zeros(Int, n)

    mat = zeros(Int, 4, n)
    mat[1, :] .= S
    mat[2, :] .= E
    mat[3, :] .= I
    mat[4, :] .= R
    return mat
end

