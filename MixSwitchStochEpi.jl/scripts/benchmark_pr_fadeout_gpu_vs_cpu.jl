#!/usr/bin/env julia

using Random
using Statistics
using Dates
using MixSwitchStochEpi
using CSV
using DataFrames

function parseargs()
    n_sims = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 5000
    reps = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 3
    return n_sims, reps
end

function make_example_params(; n_activity=6, pop_per_class=Int(1_000), ε=0.2, ξ=0.1)
    # basic scalar parms (get_scalar_parms is not exported; call via module)
    scalar = MixSwitchStochEpi.get_scalar_parms()
    β, σ, γ, ω = scalar.β, scalar.σ, scalar.γ, scalar.ω

    act_levels = collect(range(1.0, stop=2.0, length=n_activity))
    class_sizes = fill(pop_per_class, n_activity)
    pop_size = sum(class_sizes)

    M = make_garnett_contact(ε, n_activity, act_levels, class_sizes)
    W = uniform_switching(ξ, n_activity)

    params = SEIRSStochParams(β, σ, γ, ω, M, W, n_activity, class_sizes, act_levels, pop_size)
    return params
end

function make_benchmark_initial_state(params::SEIRSStochParams; num_initial=1)
    n = params.n_activity
    init = zeros(Int, 4, n)
    # place susceptibles across classes
    for i in 1:n
        init[1, i] = params.class_sizes[i]
    end
    # seed a single infected in the middle class
    mid = cld(n, 2)
    init[1, mid] -= num_initial
    init[3, mid] += num_initial
    return init
end

function time_call(f)
    GC.gc()
    t0 = time_ns()
    f()
    t1 = time_ns()
    return (t1 - t0) / 1e9
end

function run_benchmark(n_sims::Int, reps::Int)
    params = make_example_params()
    init = make_benchmark_initial_state(params; num_initial=1)

    results = DataFrame(backend=String[], rep=Int[], elapsed_s=Float64[], n_sims=Int[], n_activity=Int[], timestamp=String[])

    @info "Starting benchmark: n_sims=$n_sims, reps=$reps, n_activity=$(params.n_activity)"

    # CPU runs
    for rep in 1:reps
        rng = MersenneTwister(1234 + rep)
        elapsed = time_call(()-> pr_fadeout(params, init; num_initial=1, n_sims=n_sims, t_cutoff=30.0, rng=rng))
        push!(results, ("cpu", rep, elapsed, n_sims, params.n_activity, string(Dates.now())))
        @info "CPU rep $rep => $(round(elapsed, digits=4)) s"
    end

    # GPU runs: choose backend autodetect via pr_fadeout_gpu default
    for rep in 1:reps
        seed_base = UInt64(4321 + rep)
        elapsed = time_call(()-> pr_fadeout_gpu(params, init; n_sims=n_sims, tmax=30.0, num_initial=1, seed_base=seed_base))
        push!(results, ("gpu", rep, elapsed, n_sims, params.n_activity, string(Dates.now())))
        @info "GPU rep $rep => $(round(elapsed, digits=4)) s"
    end

    # summary
    for backend in unique(results.backend)
        times = results.elapsed_s[results.backend .== backend]
        @info "$backend summary: mean=$(round(mean(times), digits=4)) s; std=$(round(std(times), digits=4)) s; min=$(minimum(times)) s"
    end

    # ensure output dir exists
    outdir = joinpath(@__DIR__, "../output")
    mkpath(outdir)
    outfile = joinpath(outdir, "benchmark_pr_fadeout.csv")
    CSV.write(outfile, results)
    @info "Wrote results to $outfile"
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    n_sims, reps = parseargs()
    run_benchmark(n_sims, reps)
end
