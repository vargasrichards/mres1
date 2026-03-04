"""
Simple benchmarking utilities to collect timing samples during scans.

Provides:
- `PerfSamples()` - container for timing vectors
- `record!(perf, key, seconds)` - append a timing
- `timeit(perf, key, fn)` - run fn(), record elapsed seconds and return (result, seconds)
- `summary_df(perf)` - return a DataFrame with mean/median/min/max/count per key
- `save_perf(perf, path)` - write summary DataFrame to CSV

This is intentionally lightweight and depends only on base + DataFrames/CSV/Statistics.
"""

using DataFrames
using CSV
using Statistics

struct PerfSamples
    samples::Dict{String, Vector{Float64}}
    lock::ReentrantLock
end

PerfSamples() = PerfSamples(Dict{String, Vector{Float64}}(), ReentrantLock())

function record!(p::PerfSamples, key::AbstractString, dt::Float64)
    lock(p.lock) do
        v = get!(p.samples, String(key), Float64[])
        push!(v, dt)
    end
    return nothing
end

function timeit(p::PerfSamples, key::AbstractString, f::Function)
    t0 = time()
    val = f()
    dt = time() - t0
    record!(p, key, dt)
    return (val, dt)
end

function summary_df(p::PerfSamples)
    keys = String[]
    means = Float64[]
    medians = Float64[]
    mins = Float64[]
    maxs = Float64[]
    counts = Int[]
    for (k, v) in p.samples
        push!(keys, k)
        push!(means, mean(v))
        push!(medians, median(v))
        push!(mins, minimum(v))
        push!(maxs, maximum(v))
        push!(counts, length(v))
    end
    return DataFrame(metric = keys, mean = means, median = medians, min = mins, max = maxs, count = counts)
end

function save_perf(p::PerfSamples, path::AbstractString)
    df = summary_df(p)
    CSV.write(path, df)
    return nothing
end


"""
    bench!(p::PerfSamples, key::AbstractString, f::Function; samples::Int=5)

Run a small micro-benchmark of `f()` using BenchmarkTools if it is available.
If BenchmarkTools is not installed, fall back to running `f()` `samples` times and
recording elapsed times. The measured median time is recorded under `key` in `p`.
Returns the timing in seconds (median) and the raw BenchmarkTools.Trial when available.
"""
function bench!(p::PerfSamples, key::AbstractString, f::Function; samples::Int = 5)
    # Lightweight fallback benchmarking: run f() `samples` times and record timings.
    times = Float64[]
    for i in 1:samples
        t0 = time(); f(); push!(times, time() - t0)
    end
    med = median(times)
    record!(p, key, med)
    return (med, nothing)
end


"""
    bench_full!(p::PerfSamples, key::AbstractString, f::Function; samples::Int=50)

Run a fuller microbenchmark using BenchmarkTools if it is available. If BenchmarkTools
is present this attempts to run `@benchmark` and `@belapsed` via runtime `eval` so the
macro expansion happens after `using BenchmarkTools`. If BenchmarkTools is not
available the function falls back to the simple timing loop.

The function records the median time (seconds) into `p` under `key` and also returns
the median and any BenchmarkTools trial object (or `nothing` when not available).
"""
function bench_full!(p::PerfSamples, key::AbstractString, f::Function; samples::Int = 50)
    # Try to use BenchmarkTools when present
    try
        # import BenchmarkTools at runtime
        eval(:(using BenchmarkTools))

        # create a global binding for the function so macros evaluated in module scope can call it
        fname = gensym(:_benchfun)
        @eval const $(fname) = $f

        # run @belapsed to get a fast median estimate (returns seconds)
        belapsed_expr = Meta.parse("@belapsed $(fname)() samples=$samples")
        t = eval(belapsed_expr)

        # run @benchmark to get a Trial object (may be large); capture if desired
        try
            bench_expr = Meta.parse("@benchmark $(fname)() samples=$samples")
            trial = eval(bench_expr)
        catch
            trial = nothing
        end

        record!(p, key, float(t))
        return (float(t), trial)
    catch err
        @warn "BenchmarkTools not available or failed; falling back to simple timing: $err"
        times = Float64[]
        for i in 1:samples
            t0 = time(); f(); push!(times, time() - t0)
        end
        med = median(times)
        record!(p, key, med)
        return (med, nothing)
    end
end
