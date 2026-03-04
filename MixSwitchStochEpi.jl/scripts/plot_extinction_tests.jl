# Plot sample epidemic paths for the extinction-related test scenarios.
# Saves plots to ./output/

using MixSwitchStochEpi
using Plots
using Random
using Statistics
using DataFrames

function plot_paths(params, initial_state; sample_n=10, tmax=Inf, outpath::String)
    p = plot(xlabel="Time / days", ylabel="E + I", title=basename(outpath))
    for i in 1:sample_n
        rng = MersenneTwister(rand(UInt))
        res = MixSwitchStochEpi.simulate_system!(params, initial_state; rng=rng, tmax=tmax)
        times = res.times
        EIcounts = [sum(s[2, :]) + sum(s[3, :]) for s in res.states]
        plot!(p, times, EIcounts, label = "sim $i")
    end
    savefig(p, outpath)
    @info "Saved plot to: $outpath"
end

function run_all()
    outdir = joinpath(@__DIR__, "..", "output")
    isdir(outdir) || mkpath(outdir)

    # 1) Beta = 0 (extinction certain)
    params, init, tspan = MixSwitchStochEpi.default_parametrisation()
    params_zero = MixSwitchStochEpi.SEIRSStochParams(0.0, params.σ, params.γ, params.ω, params.M, params.W, params.n_activity, params.class_sizes, params.act_levels, params.pop_size)
    plot_paths(params_zero, init; sample_n=10, tmax=100.0, outpath=joinpath(outdir, "extinction_beta0.svg"))

    # 2) Branching-process scenario: single infected in class 1, uniform activity
    n = 5
    params_bp, init_bp, _ = MixSwitchStochEpi.default_parametrisation(n_activity=n, act_levels = ones(n), p_inf = 1, init_mode = :class, init_class = 1)
    # calibrate to some R0 to get interesting dynamics
    params_bp_cal = MixSwitchStochEpi.calibrate_parms(params_bp, 1.8)
    plot_paths(params_bp_cal, init_bp; sample_n=12, tmax=100.0, outpath=joinpath(outdir, "extinction_branching.svg"))

    # 3) Initial infected in highest-activity class (unequal activity)
    act_levels = [1.0, 2.0, 3.0, 5.0, 10.0]
    params_uneq, _, _ = MixSwitchStochEpi.default_parametrisation(act_levels = act_levels, n_activity = length(act_levels))
    init_high = MixSwitchStochEpi.default_parametrisation(act_levels = act_levels, n_activity = length(act_levels), init_mode = :class, init_class = length(act_levels))[2]
    params_uneq_cal = MixSwitchStochEpi.calibrate_parms(params_uneq, 1.8)
    plot_paths(params_uneq_cal, init_high; sample_n=12, tmax=300.0, outpath=joinpath(outdir, "extinction_high_activity.svg"))

    # 4) Initial infected in lowest-activity class
    init_low = MixSwitchStochEpi.default_parametrisation(act_levels = act_levels, n_activity = length(act_levels), init_mode = :class, init_class = 1)[2]
    plot_paths(params_uneq_cal, init_low; sample_n=12, tmax=300.0, outpath=joinpath(outdir, "extinction_low_activity.svg"))

    # 5) No-switch vs moderate-switch comparison (plot both sets in a single figure)
    params_uniform, init_u, _ = MixSwitchStochEpi.default_parametrisation(n_activity=5, act_levels = ones(5), p_inf = 5, init_mode=:uniform)
    params_no_switch = MixSwitchStochEpi.SEIRSStochParams(params_uniform.β, params_uniform.σ, params_uniform.γ, params_uniform.ω, params_uniform.M, MixSwitchStochEpi.uniform_switching(0.0, params_uniform.n_activity), params_uniform.n_activity, params_uniform.class_sizes, params_uniform.act_levels, params_uniform.pop_size)
    params_switch = MixSwitchStochEpi.SEIRSStochParams(params_uniform.β, params_uniform.σ, params_uniform.γ, params_uniform.ω, params_uniform.M, MixSwitchStochEpi.uniform_switching(0.5, params_uniform.n_activity), params_uniform.n_activity, params_uniform.class_sizes, params_uniform.act_levels, params_uniform.pop_size)

    p = plot(layout=(2,1), size=(800,800))
    for i in 1:8
        rng = MersenneTwister(rand(UInt))
        res1 = MixSwitchStochEpi.simulate_system!(params_no_switch, init_u; rng=rng, tmax=300.0)
        times1 = res1.times
        EI1 = [sum(s[2,:]) + sum(s[3,:]) for s in res1.states]
        plot!(p[1], times1, EI1, label = i==1 ? "no switch" : "")

        rng2 = MersenneTwister(rand(UInt))
        res2 = MixSwitchStochEpi.simulate_system!(params_switch, init_u; rng=rng2, tmax=300.0)
        times2 = res2.times
        EI2 = [sum(s[2,:]) + sum(s[3,:]) for s in res2.states]
        plot!(p[2], times2, EI2, label = i==1 ? "switch" : "")
    end
    plot!(p[1], xlabel="Time / days", ylabel="E+I", title="No switching sample paths")
    plot!(p[2], xlabel="Time / days", ylabel="E+I", title="Moderate switching sample paths")
    savefig(p, joinpath(outdir, "extinction_switch_comparison.svg"))
    @info "Saved switch comparison to output/extinction_switch_comparison.svg"

    @info "All plots saved to $outdir"
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_all()
end
