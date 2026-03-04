"""
    run_haslemere_analysis(; kwargs...)

Run the Haslemere extinction scans from inside the package. This wraps the
existing `scripts/Haslemere_analysis.jl` behaviour in a callable function so it
can be precompiled as part of the package or invoked programmatically.

Keyword arguments mirror the script defaults and are documented in the
wrapper's signature below.
"""
function run_haslemere_analysis(
    ;
    n_activity::Int = 5,
    Ntot::Int = 500000,
    n_eps::Int = 10,
    n_xi::Int = 10,
    n_sims::Int = 1000,
    t_cutoff::Float64 = 50.0,
    num_initial::Int = 1,
    r0_list::Vector{Float64} = [1.1, 1.5, 2.5, 3.5],
    act_levels::Union{Nothing, Vector{Float64}} = nothing,
    output_base::Union{Nothing, String} = nothing,
    use_caffeinate::Bool = Sys.isapple(),
)

    # default activity scores if not supplied
    if act_levels === nothing
        act_levels = gamma_bin_expectations(shape = 0.9065905, rate = 0.2610179, n_activity = n_activity)
        @info "Computed default Haslemere activity levels (length=$(length(act_levels)))"
    else
        @info "Using supplied activity levels (length=$(length(act_levels)))"
    end

    # output directory
    out_base = output_base === nothing ? joinpath(@__DIR__, "..", "output", "haslemere_stochastic_extinction") : output_base
    mkpath(out_base)

    # start caffeinate early if requested (macOS only)
    caffeinate_io = nothing
    if use_caffeinate && Sys.isapple()
        try
            caffeinate_io = open(`caffeinate -dimsu`, "r")
            @info "Started caffeinate to prevent sleep for the duration of run_haslemere_analysis."
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

    @assert Threads.nthreads() > 1 "This function benefits from multiple threads. Start Julia with multiple threads."

    # run scans for each r0 with low/high initial seeds
    for r0 in r0_list
        out_low = joinpath(out_base, "r0_$(r0)_low")
        out_high = joinpath(out_base, "r0_$(r0)_high")
        mkpath(out_low); mkpath(out_high)

        @info "Running r0=$r0; low-activity seed -> $out_low"
        scan_switch_assortativity(
            n_activity = n_activity,
            Ntot = Ntot,
            target_R0 = r0,
            ε_ref = 0.,
            n_eps = n_eps,
            n_xi = n_xi,
            n_sims = n_sims,
            t_cutoff = t_cutoff,
            act_levels = act_levels,
            outdir = out_low,
            init_class = 1,
            num_initial = num_initial,
        )

        @info "Running r0=$r0; high-activity seed -> $out_high"
        scan_switch_assortativity(
            n_activity = n_activity,
            Ntot = Ntot,
            target_R0 = r0,
            ε_ref = 0.,
            n_eps = n_eps,
            n_xi = n_xi,
            n_sims = n_sims,
            t_cutoff = t_cutoff,
            act_levels = act_levels,
            outdir = out_high,
            init_class = n_activity,
            num_initial = num_initial,
        )
    end

    # combine outputs like the script does
    all_dfs = DataFrame[]
    for r0 in r0_list
        for seed in ("low", "high")
            dir = joinpath(out_base, "r0_$(r0)_$(seed)")
            csvf = joinpath(dir, "extinction_scan_grid.csv")
            if isfile(csvf)
                df = CSV.read(csvf, DataFrame)
                df[!, :r0] = fill(r0, nrow(df))
                df[!, :seed] = fill(seed, nrow(df))
                push!(all_dfs, df)
            else
                @warn "Missing expected output file: $csvf"
            end
        end
    end

    if !isempty(all_dfs)
        combined = vcat(all_dfs...)
        combined_file = joinpath(out_base, "combined_extinction_scans.csv")
        CSV.write(combined_file, combined)
        @info "Wrote combined results to $combined_file"
    else
        @warn "No per-run CSVs found to combine."
    end

    return true
end
