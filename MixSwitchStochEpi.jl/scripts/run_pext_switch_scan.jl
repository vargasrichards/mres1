# run_pext_switch_scan.jl
#
# For each target R0 (calibrated at ξ=0), and for each switching rate
# ξ ∈ {0, 0.1, 0.2}  (uniform switching), sweep over increasing numbers
# of initially-exposed individuals and record Pr(extinction).
#
# Two init_class values are used:
#   init_class = 1  (lowest-activity class)
#   init_class = n_activity  (highest-activity class)
#
# Uses Haslemere contact patterns and activity scores.
# Results are written to output/pext_switch_scan.csv
#
# Multi-threaded for speedup across (R0, ξ, init_class) combinations.
#
# Usage:
#   julia -t 10 --project=. scripts/run_pext_switch_scan.jl

using MixSwitchStochEpi
using Random
using DataFrames
using CSV
using Base.Threads
using Statistics

# simulation parameters
const n_activity   = 5
const NTOT         = 50_000
const N_SIMS       = 1000     # replicates per (E0, init_class, ξ, R0) point
const T_CUTOFF     = 500.0      # time horizon for pr_fadeout (must be long enough
                                # for subcritical chains seeded with many E0 to die)
const THRESHOLD_FR = 0.05       # secondary-case fraction for fadeout definition

const target_r0s   = [1.5, 2.0, 2.5]
const switch_rates = [0.0, 0.1, 0.2]
const INIT_CLASSES = [1, n_activity]   # lowest & highest activity classes

# Haslemere activity scores and contact matrix
const HASLEMERE_ACT_LEVELS = [0.3146078, 1.1479220, 2.3361036, 4.2211382, 9.3466695]
const HASLEMERE_CONTACT_MATRIX = [
    0.09678707 0.1164464 0.1585146 0.1789701 0.4492819;
    0.05356591 0.2486653 0.2015198 0.2065853 0.2896636;
    0.04424087 0.1228340 0.2407142 0.2317050 0.3605060;
    0.05088035 0.1183536 0.2307819 0.2154986 0.3844855;
    0.05177194 0.0595667 0.1716733 0.2396574 0.4773307
]
const MAX_E0       = 50     
const USE_CAFFEINATE = true   

function hom_actual_r0(activity_levels, parms::SEIRSStochParams)
    return mean(activity_levels) * parms.β / parms.γ
end

function main()
    @info "Running with $(Threads.nthreads()) threads"

    caffeinate_io = nothing
    if USE_CAFFEINATE && Sys.isapple()
        try
            caffeinate_io = open(`caffeinate -dimsu`, "r")
            @info "Started caffeinate to prevent sleep during the run."
        catch e
            @warn "Could not start caffeinate: $e"
        end
    end

    atexit() do
        if caffeinate_io !== nothing
            try
                close(caffeinate_io)
                @info "Stopped caffeinate."
            catch e
                @warn "Error closing caffeinate: $e"
            end
        end
    end

    act_levels_inhom = HASLEMERE_ACT_LEVELS
    M_inhom          = HASLEMERE_CONTACT_MATRIX
    class_sizes      = fill(Int(NTOT / n_activity), n_activity)

    # homogeneous parametrisation: mean activity, proportionate mixing
    mean_act          = mean(HASLEMERE_ACT_LEVELS)
    act_levels_hom    = fill(mean_act, n_activity)
    proportionate_mix = fill(1 / n_activity, n_activity, n_activity)

    W_zero = MixSwitchStochEpi.uniform_switching(0.0, n_activity)
    σ = 1/3
    γ = 1/4

    # calibrate β to hit each target R0 using the inhomogeneous model,
    # then copy that β to the homogeneous model
    calibrated_betas_inhom = Dict{Float64, Float64}()
    calibrated_betas_hom   = Dict{Float64, Float64}()
    for r0 in target_r0s
        draft = SEIRSStochParams(0.5, σ, γ, 0.0,
                                 M_inhom, W_zero, n_activity,
                                 class_sizes, act_levels_inhom, NTOT)
        β_cal = calibrate_parms(draft, r0).β
        calibrated_betas_inhom[r0] = β_cal
        calibrated_betas_hom[r0]   = β_cal
        @info "calibrated β=$β_cal for R0=$r0"
    end

    tasks = Tuple[]
    for target_r0 in target_r0s, ξ in switch_rates, ic in INIT_CLASSES
        push!(tasks, (r0=target_r0, xi=ξ, init_class=ic, model_type=:inhom))
        push!(tasks, (r0=target_r0, xi=ξ, init_class=ic, model_type=:hom))
    end
    @info "Processing $(length(tasks)) parameter combinations in parallel"

    all_rows    = Vector{Vector{NamedTuple}}(undef, length(tasks))
    task_counter = Threads.Atomic{Int}(0)

    Threads.@threads for task_idx in 1:length(tasks)
        task       = tasks[task_idx]
        target_r0  = task.r0
        ξ          = task.xi
        ic         = task.init_class
        model_type = task.model_type

        W = MixSwitchStochEpi.uniform_switching(ξ, n_activity)
        if model_type == :inhom
            β_cal   = calibrated_betas_inhom[target_r0]
            M_use   = M_inhom
            act_use = act_levels_inhom
        else
            β_cal   = calibrated_betas_hom[target_r0]
            M_use   = proportionate_mix
            act_use = act_levels_hom
        end

        parms     = SEIRSStochParams(β_cal, 1/3, 1/4, 0.0,
                                     M_use, W, n_activity,
                                     class_sizes, act_use, NTOT)
        actual_r0 = compute_R0(parms)

        tid       = Threads.threadid()
        completed = Threads.atomic_add!(task_counter, 1)
        @info "[Thread $tid] task $(completed+1)/$(length(tasks)): R0=$target_r0 ξ=$ξ ic=$ic model=$model_type (actual_R0=$(round(actual_r0; digits=3)))"

        local_rows = NamedTuple{(:E0, :p_extinction, :init_class, :switch_rate, :R0, :model_type),
                                Tuple{Int, Float64, Int, Float64, Float64, Symbol}}[]

        for e0 in 1:min(MAX_E0, NTOT)
            init = make_initial_state(class_sizes, e0; init_mode=:class, init_class=ic)
            pext = pr_fadeout(parms, init;
                              num_initial        = e0,
                              n_sims             = N_SIMS,
                              t_cutoff           = T_CUTOFF,
                              threshold_fraction = THRESHOLD_FR,
                              rng                = Random.default_rng())
            push!(local_rows, (E0=e0, p_extinction=pext, init_class=ic,
                               switch_rate=ξ, R0=target_r0, model_type=model_type))
            if pext < THRESHOLD_FR
                @info "[Thread $tid]   stopped at E0=$e0 (p_ext=$(round(pext; digits=4)))"
                break
            end
        end

        all_rows[task_idx] = local_rows
        @info "[Thread $tid] done: R0=$target_r0 ξ=$ξ ic=$ic model=$model_type"
    end

    rows = reduce(vcat, all_rows)
    @info "Collected $(length(rows)) total data points"

    outdir  = joinpath(@__DIR__, "..", "output")
    isdir(outdir) || mkpath(outdir)
    outfile = joinpath(outdir, "pext_switch_scan_quicker.csv")
    CSV.write(outfile, DataFrame(rows))
    @info "Results written to $outfile"
end

main()

# Rscript scripts/plotters/plot_pext_switch_scan.R