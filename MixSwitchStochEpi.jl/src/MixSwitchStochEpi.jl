module MixSwitchStochEpi


include("activity_scoring.jl")
include("contact_models.jl")
include("switch_models.jl")
include("parameter_store.jl")
include("SEIRS_stochastic.jl")
include("initial_conditions.jl")
include("plot_model.jl")
include("helper_functions.jl")
include("ngm.jl")
include("benchmarking.jl")
include("extinction_scan_grid.jl")
include("extinction_threshold.jl")
include("SEIRS_gpu.jl")
include("haslemere_analysis.jl")

export make_garnett_contact,
       constrained_switching,
       find_steadystate,
       validate_generator,
       uniform_switching,
       gamma_bin_expectations,
       calibrate_parms,
       ρMat,
       Tmat,
       Tmat_entry,
       Sigma_mat,
       Kmat,
       compute_R0,
       SEIRSStochParams,
    simulate_system!,
    many_simulations,
    pr_fadeout,
    make_initial_state,
    PerfSamples,
    record!,
    timeit,
    summary_df,
    save_perf,
    bench!,
    bench_full!,
    default_parametrisation,
    make_parset,
    scan_switch_assortativity,
    plot_pext_initial,
    find_threshold,
    run_haslemere_analysis,
    gpu_backend,
    pr_fadeout_gpu,
    simulate_system_gpu,
    many_simulations_gpu


end # module MixSwitchStochEpi
