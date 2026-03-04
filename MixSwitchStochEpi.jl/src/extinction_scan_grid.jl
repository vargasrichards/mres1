# Scan over switch rate (ξ) and assortativity (ε) and compute R0 and early-time extinction probability.
# X-axis: assortativity ε (0..1), Y-axis: residence_time = 1/ξ (use large value for ξ=0).
 
using MixSwitchStochEpi
using Random
using Base.Threads
using DataFrames
using CSV
using Plots
using Statistics

# Guard so a lightweight JIT warm-up runs at most once per Julia process.
const _warmup_done_lock = ReentrantLock()
const _warmup_done = Ref(false)



"""
    examine_extremal(n_initial::Int, class_initial, num_trials)

Scan Pr(extinction) for both most- and least-active classes.

This then naturally allows computing the difference between the two scan_switch_assortativity
-> a reflection of the sensitivity to the initially infected individual(s)


# Arguments

# Returns
A dataframe with the results. Also contains the diffs of the extremal cases. That is,
the Pr(extinction | Initial infected is lowest activity class) - Pr(extinction | Initial infected is highest activity class) 
for each point in parameter space examined.


"""
function examine_extremal(num_initial::Int, 
                            target_R0::Float64, 
                            num_classes::Int, 
                            n_sims::Int)
    lowest_activity = Int(1)
    highest_activity = num_classes 

    low_run = scan_switch_assortativity(
    n_activity = num_classes,
    Ntot = 50000,
    target_R0 = target_R0,
    ε_ref = 0.,
    n_eps = 25,
    n_xi = 25,
    n_sims = n_sims,
    t_cutoff = 100.0,
    mean_activity = 5.0,
    variance_activity = 5.0,
    outdir = "output_low",
    init_class = lowest_activity,
    num_initial  = num_initial
)


end





"""
    scan_switch_assortativity()

# Description 
Scan over switch rate (ξ) and assortativity (ε) and compute R0 and early-time extinction probability.

Note that we display as Expected residence time in days 



Default number of classes is 5
"""
function scan_switch_assortativity(; 
    n_activity::Int = 5,
    Ntot::Int = 50000, #check whether this is appropriate
    target_R0::Float64, 
    ε_ref::Float64 = 0.,
    n_eps::Int = 10,
    n_xi::Int = 10,
    n_sims::Int = 20,
    t_cutoff::Float64 = 100.0,
    mean_activity = 5.0,
    variance_activity = 5.0,
    act_levels::Union{Nothing, Vector{Float64}} = nothing,
    outdir::String = "output",
    init_class::Int = 1,
    init_mode::Symbol = :class,
    num_initial::Int = 1, #initial dose of exposed individuals
    seed_base::UInt32 = UInt32(1234567)
)
    # ensure output directory exists
    isdir(outdir) || mkpath(outdir)
    if act_levels === nothing
        act_levels = gamma_bin_expectations(
            mean = mean_activity,
            variance = variance_activity,
            n_activity = n_activity
        )
        @info "activity levels [scores] being computed as $act_levels"
    else
        @info "using supplied activity levels (length=$(length(act_levels)))"
    end

    params_base, init_base, tspan =
        default_parametrisation(
            n_activity = n_activity,
            Ntot = Ntot,
            ξ = 0., ε = ε_ref,
            act_levels = act_levels,
            init_mode = init_mode,
            init_class = init_class)
    # @info "base parameter set is $params_base"

    act_levels = params_base.act_levels
    class_sizes = params_base.class_sizes
    ξ_ref = 0.0



    M_ref = make_garnett_contact(ε_ref, n_activity, act_levels, class_sizes)
    W_ref = uniform_switching(ξ_ref, n_activity)
    @info "constructing reference contact and switching matrices at calibration point ASSORT $ε_ref SWITCH $ξ_ref."
    params_ref = SEIRSStochParams(
        params_base.β,
        params_base.σ,
        params_base.γ,
        params_base.ω,
        M_ref,
        W_ref,
        n_activity,
        class_sizes,
        act_levels,
        params_base.pop_size
    )
    β_cal = calibrate_parms(params_ref, target_R0).β
    @info "Calibrated β = $β_cal ΑΤ ε=$ε_ref, ξ=$ξ_ref to target R0=$target_R0"

    epsilons = collect(range(0.0, 1.0, length = n_eps)) #n_eps determines the 
    residences = vcat(range(30.0, 1.0, length = n_xi - 1), Inf) #n_xi is num of switch rates to look at
    xis = map(rt -> isinf(rt) ? 0.0 : 1.0 / rt, residences)
    @info "Switch rates to be scanned are $xis"
    @info "Estimating pr(extinction) from $n_sims replicates per parm point"
    R0_mat = Array{Float64}(undef, n_xi, n_eps)
    pext_mat = Array{Float64}(undef, n_xi, n_eps)

    # deterministic per-point seed generator (independent of thread count)
    # Use a small SplitMix64-style mixer to derive a reproducible 32-bit seed from integers.
    splitmix64(x::UInt64) = begin
        z = x + 0x9e3779b97f4a7c15
        z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
        z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
        return z ⊻ (z >> 31)
    end

    deterministic_seed(seed_base::UInt32, iξ::Int, iε::Int, init_class::Int) = begin
        input = (UInt64(seed_base) << 32) ⊻ (UInt64(iξ) << 16) ⊻ (UInt64(iε) << 8) ⊻ UInt64(init_class)
        z = splitmix64(input)
        return UInt32(z & 0xffffffff)
    end

    # --- Warm up JIT on all threads (run at most once per process) ---
    lock(_warmup_done_lock)
    already = _warmup_done[]
    if !already
        # mark as done while we perform the warmup so concurrent calls skip it
        _warmup_done[] = true
        unlock(_warmup_done_lock)

        @info "Warming up JIT (lightweight) on up to 4 threads..."
        # Limit the number of concurrent warm-up tasks to avoid doing heavy work
        # at startup (especially for large Ntot). Use a no-infections initial state
        # so simulate_system! / pr_fadeout return immediately while still
        # triggering JIT compilation of the hot paths.
        warmup_tasks = min(Threads.nthreads(), 4)
        init_zero = zeros(Int, 4, n_activity)
        # put everyone in S initially
        init_zero[1, :] .= class_sizes

        @sync for _ in 1:warmup_tasks
            @spawn begin
                try
                    tmp_W = uniform_switching(0.0, n_activity)
                    tmp_M = make_garnett_contact(ε_ref, n_activity, act_levels, class_sizes)
                    tmp_params = SEIRSStochParams(
                        β_cal,
                        params_base.σ,
                        params_base.γ,
                        params_base.ω,
                        tmp_M,
                        tmp_W,
                        n_activity,
                        class_sizes,
                        act_levels,
                        params_base.pop_size,
                    )
                    rng_w = MersenneTwister(UInt32(1))
                    # small calls that should return immediately because init_zero has no E/I
                    compute_R0(tmp_params)
                    simulate_system!(tmp_params, init_zero; tmax = 0.0, rng = rng_w)
                catch e
                    @debug "Warmup task failed: $e"
                end
            end
        end
    else
        unlock(_warmup_done_lock)
        @info "Warm-up already performed in this process; skipping."
    end


        # Prepare output CSV for incremental appends so long runs can be monitored/resumed
        out_csv = joinpath(outdir, "extinction_scan_grid.csv")
        # open writer task that appends CSV lines (avoid reopening file per-row)
        header_line = "ξ,ε,residence_time,R0,p_ext,iξ,iε,seed,runtime_seconds\n"

        # channel for lines to write
        write_ch = Channel{String}(1024)

        # start writer task (single writer to avoid locking overhead)
        writer = @spawn begin
            open(out_csv, "w") do io
                write(io, header_line)
                flush(io)
                for line in write_ch
                    write(io, line)
                    flush(io)
                end
            end
        end

        # Channel of work items: (iξ, iε)
        pairs = [(iξ, iε) for iξ in 1:length(xis) for iε in 1:length(epsilons)]
        ch = Channel{Tuple{Int,Int}}(length(pairs)) do c
            for p in pairs
                put!(c, p)
            end
        end

        # benchmarking collector
        perf = PerfSamples()

        @sync begin
            for worker_id in 1:Threads.nthreads()
                @spawn begin
                    for pair in ch
                        iξ, iε = pair
                        ξ = xis[iξ]
                        ε = epsilons[iε]
                        # compute contact & switching matrices
                        t_m0 = time()
                        W = uniform_switching(ξ, n_activity)
                        M = make_garnett_contact(ε, n_activity, act_levels, class_sizes)
                        dt_m = time() - t_m0
                        record!(perf, "make_garnett_contact", dt_m)

                        # seed per (iξ,iε) so results are independent of how many threads are used
                        seed = deterministic_seed(seed_base, iξ, iε, init_class)
                        rng = MersenneTwister(seed)

                        params_point = SEIRSStochParams(
                            β_cal,
                            params_base.σ,
                            params_base.γ,
                            params_base.ω,
                            M,
                            W,
                            n_activity,
                            class_sizes,
                            act_levels,
                            params_base.pop_size
                        )

                        # measure compute_R0 and pr_fadeout separately
                        t0 = time()
                        R0_val = compute_R0(params_point)
                        dt_R0 = time() - t0
                        record!(perf, "compute_R0", dt_R0)

                        t_start = time()
                        pext = pr_fadeout(params_point, init_base; num_initial = num_initial, n_sims = n_sims, t_cutoff = t_cutoff, rng = rng)
                        t_elapsed = time() - t_start
                        record!(perf, "pr_fadeout", t_elapsed)

                        # total time for the point (including matrix builds)
                        record!(perf, "point_total", t_elapsed + dt_R0 + dt_m)

                        # prepare CSV line (simple, no quoting since values are numeric)
                        line = string(ξ, ",", ε, ",", residences[iξ], ",", R0_val, ",", pext, ",", iξ, ",", iε, ",", seed, ",", t_elapsed, "\n")
                        put!(write_ch, line)

                        # store in matrices for final plotting/return
                        R0_mat[iξ, iε] = R0_val
                        pext_mat[iξ, iε] = pext
                    end
                end
            end
        end

        # close writer channel and wait for writer to finish
        close(write_ch)
        wait(writer)

    rows = NamedTuple[]
    for (iξ, ξ) in enumerate(xis), (iε, ε) in enumerate(epsilons)
        push!(rows, (
            ξ = ξ,
            ε = ε,
            residence_time = residences[iξ],
            R0 = R0_mat[iξ, iε],
            p_ext = pext_mat[iξ, iε]
        ))
    end

    df = DataFrame(rows)
    CSV.write(joinpath(outdir, "extinction_scan_grid.csv"), df)
    res_plot = copy(residences)
    finite_res = res_plot[isfinite.(res_plot)]
    res_inf_proxy = maximum(finite_res) * 1.2
    res_plot[.!isfinite.(res_plot)] .= res_inf_proxy

    # sort y-values and permute matrices
    perm = sortperm(res_plot)
    res_plot = res_plot[perm]
    R0_plot = R0_mat[perm, :]
    pext_plot = pext_mat[perm, :]

yticks_vals = vcat(sort(finite_res), res_inf_proxy)

yticks_labels = vcat(
    string.(round.(sort(finite_res); digits = 1)),
    "∞"
)
    p1 = heatmap(
        epsilons,
        res_plot,
        R0_plot;
        xlabel = "Assortativity ε",
        ylabel = "Expected (Residence time) /days",
        yticks = (yticks_vals, yticks_labels),
        title = "R₀ (calibrated at ε=0, ξ=0)",
        colorbar_title = "R₀",
        c = :viridis
    )

    savefig(p1, joinpath(outdir, "heatmap_R0.png"))
    savefig(p1, joinpath(outdir, "heatmap_R0.pdf"))

    p2 = heatmap(
        epsilons,
        res_plot,
        pext_plot;
        xlabel = "Assortativity ε",
        ylabel = "Residence time (1/ξ)",
        yticks = (yticks_vals, yticks_labels),
        title = "Pr(extinction at t=$t_cutoff days)",
        colorbar_title = "P(ext)",
        c = :viridis
    )
    fname = "heatmap_pr_ext"
    savefig(p2, joinpath(outdir, "$fname.png"))
    savefig(p2, joinpath(outdir, "$fname.pdf"))
    println("saved heatmap to $fname")

    # save benchmarking summary
    try
        save_perf(perf, joinpath(outdir, "perf_summary.csv"))
        @info "Wrote performance summary to $(joinpath(outdir, "perf_summary.csv"))"
    catch e
        @warn "Could not write perf summary: $e"
    end

    return (
        R0 = R0_mat,
        p_ext = pext_mat,
        epsilons = epsilons,
        xis = xis,
        residence_times = residences
    )
end

