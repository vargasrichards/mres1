using MixSwitchStochEpi
using Random
using DataFrames
using CSV
using Base.Threads

const N_ACTIVITY   = 5
const NTOT         = 50_000
const N_SIMS       = 1000
const T_CUTOFF     = 500.0
const THRESHOLD_FR = 0.05

const TARGET_R0S = [1.5, 2.0, 2.5]
const MAX_E0     = 50

const HASLEMERE_ACT_LEVELS = [0.3146078, 1.1479220, 2.3361036, 4.2211382, 9.3466695]

# Calibrated homogeneous model: equal contact matrix, activity levels set to the
# mean of the Haslemere distribution so β is comparable with the inhomogeneous scan.
function construct_homogeneous_case(r0, n_activity)
    mean_act = mean(HASLEMERE_ACT_LEVELS)
    hom_parms = SEIRSStochParams(r0, 1/3, 1/4, 0.0,
                                 ones(n_activity, n_activity) ./ n_activity,
                                 MixSwitchStochEpi.uniform_switching(0.0, n_activity),
                                 n_activity,
                                 fill(Int(NTOT ÷ n_activity), n_activity),
                                 fill(mean_act, n_activity),
                                 NTOT)
    return MixSwitchStochEpi.calibrate_parms(hom_parms, r0)
end

# Sweep E0 = 1..max_e0 and record Pr(extinction) for each value.
function scan_pr_ext_over_E0(parms::SEIRSStochParams, target_r0;
                              n_sims::Int    = N_SIMS,
                              t_cutoff::Float64 = T_CUTOFF,
                              max_e0::Int    = MAX_E0)
    rows = NamedTuple{(:E0, :p_extinction, :init_class, :switch_rate, :R0),
                      Tuple{Int, Float64, Int, Float64, Float64}}[]
    for e0 in 1:min(max_e0, parms.pop_size)
        init = MixSwitchStochEpi.make_initial_state(parms.class_sizes, e0;
                                                    init_mode = :class, init_class = 1)
        pext = MixSwitchStochEpi.pr_fadeout(parms, init;
                                            num_initial        = e0,
                                            n_sims             = n_sims,
                                            t_cutoff           = t_cutoff,
                                            threshold_fraction = THRESHOLD_FR,
                                            rng                = Random.default_rng())
        push!(rows, (E0 = e0, p_extinction = pext,
                     init_class = 1, switch_rate = 0.0, R0 = target_r0))
    end
    return rows
end

function main(r0_values = TARGET_R0S)
    @info "Running homogeneous case scan with $(Threads.nthreads()) threads"

    all_rows = NamedTuple{(:E0, :p_extinction, :init_class, :switch_rate, :R0),
                          Tuple{Int, Float64, Int, Float64, Float64}}[]

    for r0 in r0_values
        @info "Processing R0=$r0..."
        parms = construct_homogeneous_case(r0, N_ACTIVITY)
        rows  = scan_pr_ext_over_E0(parms, r0)
        append!(all_rows, rows)
    end

    outdir = joinpath(@__DIR__, "..", "output")
    isdir(outdir) || mkpath(outdir)

    outfile = joinpath(outdir, "pext_homogeneous_scan.csv")
    CSV.write(outfile, DataFrame(all_rows))
    @info "Results written to $outfile"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
