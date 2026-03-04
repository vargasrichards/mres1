


"""

"""
function run_sample(params::SEIRSStochParams, initial_state, n_reps::Int)
    sample_n = 10
    t_early = 50.0
    
    out_dir = "output"
    isdir(out_dir) || mkdir(out_dir)

    p = plot(xlabel="Time / days", ylabel="E + I", title="Sample epidemic paths")

    for i in 1:sample_n
        rng = MersenneTwister(rand(UInt))
        res = simulate_system!(params, initial_state; rng = rng)
        times = res.times
        EIcounts = [sum(s[2, :]) + sum(s[3, :]) for s in res.states]
        plot!(p, times, EIcounts, label = "sim $i")
    end
    svg_path = joinpath(out_dir, "samplepaths.svg")
    savefig(p, svg_path)
    println("Saved sample stochastic paths to $svg_path")

    # Estimate early extinction probability at t_early
    n_ext = 0
    for rep in 1:n_reps
        rng = MersenneTwister(rand(UInt))
        res = simulate_system!(params, initial_state; rng = rng, tmax = t_early)
        final_state = res.states[end]
        total_inf = sum(final_state[2, :]) + sum(final_state[3, :])
        if total_inf == 0
            n_ext += 1
        end
    end
    p_ext = n_ext / n_reps
    println("Estimated early extinction probability at t=$t_early days: $p_ext (n_reps=$n_reps)")

    df = DataFrame()
    df[!, :t_early] = [t_early]
    df[!, :n_reps] = [n_reps]
    df[!, :p_ext] = [p_ext]
    CSV.write(joinpath(out_dir, "sample_extinction_summary.csv"), df)
    println("Saved summary to output/sample_extinction_summary.csv")
    return p_ext
end

# # run when included

