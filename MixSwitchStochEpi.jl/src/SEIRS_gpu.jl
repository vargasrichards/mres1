# GPU-accelerated Gillespie SSA for the class-structured SEIRS model.
#
# Uses KernelAbstractions.jl for platform-agnostic GPU compute:
#   - Metal (Apple Silicon) via Metal.jl
#   - CUDA (NVIDIA)         via CUDA.jl
#   - ROCm (AMD)            via AMDGPU.jl
#   - CPU fallback           (always available)
#
# One thread = one complete Gillespie trajectory.  Only final-state
# scalars are returned (no dynamic allocation on device).
#
# All kernel arithmetic is Float32 (required by Metal).
# Maximum n_activity = 16.
#
# A. Vargas Richards -- Feb 2026

using KernelAbstractions

# GPU-compatible xoshiro128+ PRNG  (state = 4 x UInt32)

struct GPURng
    s0::UInt32
    s1::UInt32
    s2::UInt32
    s3::UInt32
end

@inline function _rotl(x::UInt32, k::Int)
    (x << k) | (x >> (32 - k))
end

@inline function gpu_next(rng::GPURng)
    result = rng.s0 + rng.s3
    t = rng.s1 << 9
    s2 = rng.s2 ⊻ rng.s0
    s3 = rng.s3 ⊻ rng.s1
    s1 = rng.s1 ⊻ s2
    s0 = rng.s0 ⊻ s3
    s2 = s2 ⊻ t
    s3 = _rotl(s3, 11)
    return GPURng(s0, s1, s2, s3), result
end

@inline function gpu_rand(rng::GPURng)
    rng2, v = gpu_next(rng)
    return rng2, Float32(v) * Float32(2.3283064f-10)
end

@inline function gpu_randexp(rng::GPURng)
    rng2, u = gpu_rand(rng)
    u = max(u, Float32(1.0f-30))
    return rng2, -log(u)
end

@inline function seed_rng(seed::UInt64)
    z = seed
    z = (z ⊻ (z >> 30)) * UInt64(0xbf58476d1ce4e5b9)
    z = (z ⊻ (z >> 27)) * UInt64(0x94d049bb133111eb)
    z = z ⊻ (z >> 31)
    s0 = UInt32(z & 0xffffffff)
    s1 = UInt32((z >> 32) & 0xffffffff)
    z = (seed + UInt64(0x9E3779B97F4A7C15))
    z = (z ⊻ (z >> 30)) * UInt64(0xbf58476d1ce4e5b9)
    z = (z ⊻ (z >> 27)) * UInt64(0x94d049bb133111eb)
    z = z ⊻ (z >> 31)
    s2 = UInt32(z & 0xffffffff)
    s3 = UInt32((z >> 32) & 0xffffffff)
    return GPURng(s0, s1, s2, s3)
end

const MAX_N_ACT = 16

# KernelAbstractions kernel

@kernel function _ssa_fadeout_kernel!(
    out_extinct,
    out_finalR,
    out_duration,
    scratch_S, scratch_E, scratch_I, scratch_R,
    M_flat, W_flat, act_flat,
    beta::Float32, sigma::Float32, gamma_r::Float32, omega::Float32,
    pop_size::Int32, n_i::Int32,
    tmax::Float32, threshold_frac::Float32, num_initial::Int32,
    seed_base::UInt64, max_events::Int32,
)
    sim = @index(Global)
    rng = seed_rng(seed_base + UInt64(sim))
    t = Float32(0)

    @inbounds begin

    evt = Int32(0)
    while evt < max_events
        evt += Int32(1)

        # class totals
        esum = Int32(0); isum = Int32(0)
        ssum = Int32(0); rsum = Int32(0)
        for c in Int32(1):n_i
            ssum += scratch_S[sim, c]
            esum += scratch_E[sim, c]
            isum += scratch_I[sim, c]
            rsum += scratch_R[sim, c]
        end
        (esum + isum) == Int32(0) && break

        # FOI per class
        rate_inf = Float32(0)
        for c in Int32(1):n_i
            denom_f = Float32(0); numer_f = Float32(0)
            for j in Int32(1):n_i
                Nj = Float32(scratch_S[sim, j] + scratch_E[sim, j] +
                             scratch_I[sim, j] + scratch_R[sim, j])
                w = M_flat[(j - Int32(1)) * n_i + c] * Nj
                denom_f += w
                if Nj > Float32(0)
                    numer_f += w * (Float32(scratch_I[sim, j]) / Nj)
                end
            end
            lam_c = denom_f > Float32(0) ? (beta * act_flat[c] * numer_f / denom_f) : Float32(0)
            rate_inf += lam_c * Float32(scratch_S[sim, c])
        end

        rate_EI = sigma * Float32(esum)
        rate_IR = gamma_r * Float32(isum)
        rate_RS = omega * Float32(rsum)

        # movement rates
        rate_move = Float32(0)
        for c in Int32(1):n_i, j in Int32(1):n_i
            if c != j
                wij = W_flat[(j - Int32(1)) * n_i + c]
                rate_move += wij * Float32(scratch_S[sim, c])
                rate_move += wij * Float32(scratch_E[sim, c])
                rate_move += wij * Float32(scratch_I[sim, c])
                rate_move += wij * Float32(scratch_R[sim, c])
            end
        end

        total_rate = rate_inf + rate_EI + rate_IR + rate_RS + rate_move
        total_rate <= Float32(0) && break

        # sample dt
        rng, exp_val = gpu_randexp(rng)
        dt = exp_val / total_rate
        if t + dt > tmax
            t = tmax; break
        end
        t += dt

        # choose event
        rng, u = gpu_rand(rng)
        r = u * total_rate

        if r < rate_inf
            # infection
            rng, u2 = gpu_rand(rng)
            thresh_inf = u2 * rate_inf
            acc = Float32(0)
            for c in Int32(1):n_i
                denom2 = Float32(0); numer2 = Float32(0)
                for j in Int32(1):n_i
                    Nj2 = Float32(scratch_S[sim, j] + scratch_E[sim, j] +
                                  scratch_I[sim, j] + scratch_R[sim, j])
                    w2 = M_flat[(j - Int32(1)) * n_i + c] * Nj2
                    denom2 += w2
                    if Nj2 > Float32(0)
                        numer2 += w2 * (Float32(scratch_I[sim, j]) / Nj2)
                    end
                end
                lam_c2 = denom2 > Float32(0) ? (beta * act_flat[c] * numer2 / denom2) : Float32(0)
                acc += lam_c2 * Float32(scratch_S[sim, c])
                if acc >= thresh_inf
                    if scratch_S[sim, c] > Int32(0)
                        scratch_S[sim, c] -= Int32(1)
                        scratch_E[sim, c] += Int32(1)
                    end
                    break
                end
            end

        elseif r < rate_inf + rate_EI
            # E -> I
            rng, u2 = gpu_rand(rng)
            thresh_ei = u2 * rate_EI
            acc = Float32(0)
            for c in Int32(1):n_i
                acc += sigma * Float32(scratch_E[sim, c])
                if acc >= thresh_ei
                    scratch_E[sim, c] -= Int32(1)
                    scratch_I[sim, c] += Int32(1)
                    break
                end
            end

        elseif r < rate_inf + rate_EI + rate_IR
            # I -> R
            rng, u2 = gpu_rand(rng)
            thresh_ir = u2 * rate_IR
            acc = Float32(0)
            for c in Int32(1):n_i
                acc += gamma_r * Float32(scratch_I[sim, c])
                if acc >= thresh_ir
                    scratch_I[sim, c] -= Int32(1)
                    scratch_R[sim, c] += Int32(1)
                    break
                end
            end

        elseif r < rate_inf + rate_EI + rate_IR + rate_RS
            # R -> S
            rng, u2 = gpu_rand(rng)
            thresh_rs = u2 * rate_RS
            acc = Float32(0)
            for c in Int32(1):n_i
                acc += omega * Float32(scratch_R[sim, c])
                if acc >= thresh_rs
                    scratch_R[sim, c] -= Int32(1)
                    scratch_S[sim, c] += Int32(1)
                    break
                end
            end

        else
            # movement
            rng, u2 = gpu_rand(rng)
            thresh_mv = u2 * rate_move
            acc = Float32(0)
            done = false
            for comp_id in Int32(1):Int32(4)
                done && break
                for c in Int32(1):n_i
                    done && break
                    for j in Int32(1):n_i
                        if c != j
                            wij = W_flat[(j - Int32(1)) * n_i + c]
                            count_cj = if comp_id == Int32(1)
                                scratch_S[sim, c]
                            elseif comp_id == Int32(2)
                                scratch_E[sim, c]
                            elseif comp_id == Int32(3)
                                scratch_I[sim, c]
                            else
                                scratch_R[sim, c]
                            end
                            acc += wij * Float32(count_cj)
                            if acc >= thresh_mv
                                if comp_id == Int32(1)
                                    scratch_S[sim, c] -= Int32(1)
                                    scratch_S[sim, j] += Int32(1)
                                elseif comp_id == Int32(2)
                                    scratch_E[sim, c] -= Int32(1)
                                    scratch_E[sim, j] += Int32(1)
                                elseif comp_id == Int32(3)
                                    scratch_I[sim, c] -= Int32(1)
                                    scratch_I[sim, j] += Int32(1)
                                else
                                    scratch_R[sim, c] -= Int32(1)
                                    scratch_R[sim, j] += Int32(1)
                                end
                                done = true
                                break
                            end
                        end
                    end
                end
            end
        end

    end  # while events

    # write output
    finalR = Int32(0); finalE = Int32(0); finalI = Int32(0)
    for c in Int32(1):n_i
        finalR += scratch_R[sim, c]
        finalE += scratch_E[sim, c]
        finalI += scratch_I[sim, c]
    end

    secondary = finalR - num_initial
    extinct = ((finalE + finalI) == Int32(0)) &&
              (Float32(secondary) < threshold_frac * Float32(pop_size))

    out_extinct[sim]  = extinct ? UInt8(1) : UInt8(0)
    out_finalR[sim]   = finalR
    out_duration[sim] = t

    end # @inbounds
end


# Host-side helpers

function gpu_backend()
    # Metal (Apple Silicon)
    try
        Metal = Base.require(Base.PkgId(Base.UUID("dde4c033-4e86-420c-a63e-0dd931031962"), "Metal"))
        backend = Core.eval(Metal, :(MetalBackend()))
        @info "GPU backend: Metal (Apple Silicon)"
        return backend
    catch e
        @debug "Metal backend not available" exception = (e, catch_backtrace())
    end
    # CUDA
    try
        CUDA = Base.require(Base.PkgId(Base.UUID("052768ef-5323-5732-b1bb-66c8b64840ba"), "CUDA"))
        if Core.eval(CUDA, :(functional()))
            backend = Core.eval(CUDA, :(CUDABackend()))
            @info "GPU backend: CUDA"
            return backend
        end
    catch e
        @debug "CUDA backend not available" exception = (e, catch_backtrace())
    end
    # ROCm / AMDGPU
    try
        AMDGPU = Base.require(Base.PkgId(Base.UUID("21141c5a-9bdb-4563-92ae-f87d6854732e"), "AMDGPU"))
        if Core.eval(AMDGPU, :(functional()))
            backend = Core.eval(AMDGPU, :(ROCBackend()))
            @info "GPU backend: ROCm (AMD)"
            return backend
        end
    catch e
        @debug "AMDGPU backend not available" exception = (e, catch_backtrace())
    end
    @info "GPU backend: CPU (no GPU detected; using KernelAbstractions CPU backend)"
    return KernelAbstractions.CPU()
end

function _get_array_constructor(backend)
    t = string(typeof(backend))
    if occursin("Metal", t)
        Metal = Base.require(Base.PkgId(Base.UUID("dde4c033-4e86-420c-a63e-0dd931031962"), "Metal"))
        return x -> Base.invokelatest(Metal.MtlArray, x)
    elseif occursin("CUDA", t)
        CUDA = Base.require(Base.PkgId(Base.UUID("052768ef-5323-5732-b1bb-66c8b64840ba"), "CUDA"))
        return x -> Base.invokelatest(CUDA.CuArray, x)
    elseif occursin("ROC", t) || occursin("AMDGPU", t)
        AMDGPU = Base.require(Base.PkgId(Base.UUID("21141c5a-9bdb-4563-92ae-f87d6854732e"), "AMDGPU"))
        return x -> Base.invokelatest(AMDGPU.ROCArray, x)
    else
        return x -> Array(x)
    end
end


# Public API

"""
    pr_fadeout_gpu(params::SEIRSStochParams, initial_state::Matrix{Int}; n_sims, tmax, ...)

GPU-accelerated calculation of the probability of stochastic fadeout.

# Examples
```julia
using MixSwitchStochEpi, KernelAbstractions
params, init, _ = MixSwitchStochEpi.default_parametrisation()
# Run 100 simulations using the CPU backend for a quick test
res = pr_fadeout_gpu(params, init; n_sims=100, tmax=50.0, backend=KernelAbstractions.CPU())
res.p_extinct >= 0.0
```
"""
function pr_fadeout_gpu(
    params::SEIRSStochParams,
    initial_state::Matrix{Int};
    n_sims::Int = 1000,
    tmax::Float64 = 50.0,
    threshold_fraction::Float64 = 0.05,
    num_initial::Int = 1,
    seed_base::UInt64 = UInt64(12345),
    max_events::Int = 2_000_000,
    backend = nothing,
)
    n = params.n_activity
    @assert size(initial_state) == (4, n)
    @assert n <= MAX_N_ACT "n_activity must be <= $MAX_N_ACT for GPU kernel"

    be = backend === nothing ? gpu_backend() : backend
    Arr = _get_array_constructor(be)

    M_flat   = Arr(Float32.(vec(params.M)))
    W_flat   = Arr(Float32.(vec(params.W)))
    act_flat = Arr(Float32.(params.act_levels))

    scratch_S = Arr(repeat(Int32.(initial_state[1, :])', n_sims, 1))
    scratch_E = Arr(repeat(Int32.(initial_state[2, :])', n_sims, 1))
    scratch_I = Arr(repeat(Int32.(initial_state[3, :])', n_sims, 1))
    scratch_R = Arr(repeat(Int32.(initial_state[4, :])', n_sims, 1))

    out_extinct  = Arr(zeros(UInt8,   n_sims))
    out_finalR   = Arr(zeros(Int32,   n_sims))
    out_duration = Arr(zeros(Float32, n_sims))

    kernel! = _ssa_fadeout_kernel!(be)
    Base.invokelatest(kernel!,
        out_extinct, out_finalR, out_duration,
        scratch_S, scratch_E, scratch_I, scratch_R,
        M_flat, W_flat, act_flat,
        Float32(params.β), Float32(params.σ), Float32(params.γ), Float32(params.ω),
        Int32(params.pop_size), Int32(n),
        Float32(tmax), Float32(threshold_fraction), Int32(num_initial),
        seed_base, Int32(max_events);
        ndrange = n_sims,
    )
    Base.invokelatest(KernelAbstractions.synchronize, be)

    h_extinct  = Base.invokelatest(Array, out_extinct)
    h_finalR   = Base.invokelatest(Array, out_finalR)
    h_duration = Base.invokelatest(Array, out_duration)

    n_ext = Int(sum(h_extinct))
    p_ext = n_ext / n_sims

    return (
        p_extinct  = p_ext,
        n_extinct  = n_ext,
        n_sims     = n_sims,
        final_R    = h_finalR,
        durations  = h_duration,
    )
end

"""
    simulate_system_gpu(params::SEIRSStochParams, initial_state::Matrix{Int}; n_sims, tmax, ...)

GPU-accelerated SSA simulations returning final states and durations.

# Examples
```julia
using MixSwitchStochEpi, KernelAbstractions
params, init, _ = MixSwitchStochEpi.default_parametrisation()
res = simulate_system_gpu(params, init; n_sims=10, tmax=30.0, backend=KernelAbstractions.CPU())
size(res.final_states) == (10, 4, params.n_activity)
```
"""
function simulate_system_gpu(
    params::SEIRSStochParams,
    initial_state::Matrix{Int};
    n_sims::Int = 1,
    tmax::Float64 = Inf,
    seed_base::UInt64 = UInt64(42),
    max_events::Int = 5_000_000,
    backend = nothing,
)
    n = params.n_activity
    @assert size(initial_state) == (4, n)
    @assert n <= MAX_N_ACT

    be = backend === nothing ? gpu_backend() : backend
    Arr = _get_array_constructor(be)

    M_flat   = Arr(Float32.(vec(params.M)))
    W_flat   = Arr(Float32.(vec(params.W)))
    act_flat = Arr(Float32.(params.act_levels))

    scratch_S = Arr(repeat(Int32.(initial_state[1, :])', n_sims, 1))
    scratch_E = Arr(repeat(Int32.(initial_state[2, :])', n_sims, 1))
    scratch_I = Arr(repeat(Int32.(initial_state[3, :])', n_sims, 1))
    scratch_R = Arr(repeat(Int32.(initial_state[4, :])', n_sims, 1))

    out_extinct  = Arr(zeros(UInt8,   n_sims))
    out_finalR   = Arr(zeros(Int32,   n_sims))
    out_duration = Arr(zeros(Float32, n_sims))

    tmax32 = Float32(min(tmax, Float64(typemax(Float32)) / 2))

    kernel! = _ssa_fadeout_kernel!(be)
    Base.invokelatest(kernel!,
        out_extinct, out_finalR, out_duration,
        scratch_S, scratch_E, scratch_I, scratch_R,
        M_flat, W_flat, act_flat,
        Float32(params.β), Float32(params.σ), Float32(params.γ), Float32(params.ω),
        Int32(params.pop_size), Int32(n),
        tmax32, Float32(0.05), Int32(1),
        seed_base, Int32(max_events);
        ndrange = n_sims,
    )
    Base.invokelatest(KernelAbstractions.synchronize, be)

    hS = Base.invokelatest(Array, scratch_S)
    hE = Base.invokelatest(Array, scratch_E)
    hI = Base.invokelatest(Array, scratch_I)
    hR = Base.invokelatest(Array, scratch_R)
    h_dur    = Base.invokelatest(Array, out_duration)
    h_finalR = Base.invokelatest(Array, out_finalR)

    final_states = Array{Int32}(undef, n_sims, 4, n)
    final_states[:, 1, :] = hS
    final_states[:, 2, :] = hE
    final_states[:, 3, :] = hI
    final_states[:, 4, :] = hR

    return (final_states = final_states, durations = h_dur, final_R = h_finalR)
end

function many_simulations_gpu(
    params::SEIRSStochParams,
    initial_state::Matrix{Int},
    num_initial::Int,
    n_sims::Int;
    tmax::Float64 = Inf,
    seed_base::UInt64 = UInt64(42),
    max_events::Int = 5_000_000,
    backend = nothing,
)
    result = simulate_system_gpu(
        params, initial_state;
        n_sims = n_sims, tmax = tmax,
        seed_base = seed_base, max_events = max_events,
        backend = backend,
    )

    secondary_sizes = Float64.((result.final_R .- Int32(num_initial))) ./ params.pop_size
    durations = Float64.(result.durations)

    return (secondary_sizes = secondary_sizes, durations = durations)
end

"""
    find_threshold_gpu(parms::SEIRSStochParams; n_sims=1000, tmax=50.0, threshold_value=0.05, backend=nothing)

GPU-accelerated search for the minimum initial infected count required for a fadeout probability below `threshold_value`.

# Examples
```julia
using MixSwitchStochEpi, KernelAbstractions
params, _, _ = MixSwitchStochEpi.default_parametrisation(Ntot=100)
pcal = MixSwitchStochEpi.calibrate_parms(params, 2.0)
t = find_threshold_gpu(pcal; n_sims=100, tmax=40.0, threshold_value=0.1, backend=KernelAbstractions.CPU())
t > 0
```
"""
function find_threshold_gpu(parms::SEIRSStochParams; n_sims::Int = 1000, tmax::Float64 = 50.0, threshold_value::Float64 = 0.05, backend = nothing)
    pop_size = parms.pop_size
    # We can potentially parallelise this loop better in the future, 
    # but for now we reuse the parallelized pr_fadeout_gpu.
    for init_infected in 1:pop_size
        init = make_initial_state(parms.class_sizes, init_infected; init_mode = :class, init_class = 1)
    res = pr_fadeout_gpu(parms, init; n_sims=n_sims, tmax=tmax, num_initial=init_infected,
                threshold_fraction=threshold_value, backend=backend)
        if res.p_extinct < threshold_value
            return init_infected
        end
    end
    return pop_size 
end
