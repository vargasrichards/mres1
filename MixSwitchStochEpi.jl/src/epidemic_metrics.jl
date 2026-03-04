# Compute the epidemic metrics
# HIT, final sizes, peak size, peak time. 
# A. Vargas Richards Nov. 2025. 


"""
    extract_trajectories(sol, params::SEIRSStochParams)

Extract the epidemic trajectories from the simulation
"""
function extract_trajectories(sol, params::SEIRSStochParams)
    n_times = length(sol.t)
    n_classes = params.n_activity
    
    S_matrix = hcat([sol.u[i][1, :] for i in 1:n_times]...)'
    E_matrix = hcat([sol.u[i][2, :] for i in 1:n_times]...)'
    I_matrix = hcat([sol.u[i][3, :] for i in 1:n_times]...)'
    R_matrix = hcat([sol.u[i][4, :] for i in 1:n_times]...)'
    
    df = DataFrame(time = sol.t)
    
    for j in 1:n_classes
        df[!, Symbol("S_$j")] = S_matrix[:, j]
        df[!, Symbol("E_$j")] = E_matrix[:, j]
        df[!, Symbol("I_$j")] = I_matrix[:, j]
        df[!, Symbol("R_$j")] = R_matrix[:, j]
        df[!, Symbol("EI_$j")] = E_matrix[:, j] + I_matrix[:, j]
        df[!, Symbol("class_$(j)_total")] = S_matrix[:, j] + E_matrix[:, j] + I_matrix[:, j] + R_matrix[:, j]
    end
    
    df.S_total = sum(S_matrix, dims=2)[:]
    df.E_total = sum(E_matrix, dims=2)[:]
    df.I_total = sum(I_matrix, dims=2)[:]
    df.EI_total = df.E_total + df.I_total
    df.R_total = sum(R_matrix, dims=2)[:]
    df.pop_total = sum(df[!, Symbol("class_$(j)_total")] for j in 1:n_classes)
    
    return df
end


"""
Function which computes the metrics of the epidemics 
"""
function compute_epidemic_metrics(sol, params::SEIRSParams)
    df = extract_trajectories(sol, params)

    total_pop = sum(params.class_sizes)

    I_tot = df.I_total
    EI_tot = df.EI_total
    S_tot = df.S_total
    R_tot = df.R_total

    idx_ptI  = argmax(I_tot)
    idx_ptEI = argmax(EI_tot)

    overall_psI  = I_tot[idx_ptI]
    overall_ptI  = df.time[idx_ptI]

    overall_psEI = EI_tot[idx_ptEI]
    overall_ptEI = df.time[idx_ptEI]

    overall_fs   = R_tot[end] / total_pop
    overall_hit  = 1.0 - S_tot[idx_ptEI] / total_pop

    first_I = findfirst(x -> x > 0, I_tot)
    last_I  = findlast(x -> x > 1, I_tot)
    overall_duration = df.time[last_I] - df.time[first_I]

    overall_row = (
        class = -1,
        psI = overall_psI,
        ptI = overall_ptI,
        psEI = overall_psEI,
        ptEI = overall_ptEI,
        fs = overall_fs,
        hit = overall_hit,
        duration = overall_duration
    )

    class_metrics = DataFrame([overall_row])

    for j in 1:params.n_activity
        S = df[!, Symbol("S_$j")]
        E = df[!, Symbol("E_$j")]
        I = df[!, Symbol("I_$j")]
        R = df[!, Symbol("R_$j")]

        EI = E .+ I

        idx_ptI_j  = argmax(I)
        idx_ptEI_j = argmax(EI)

        class_psI  = I[idx_ptI_j]
        class_ptI  = df.time[idx_ptI_j]

        class_psEI = EI[idx_ptEI_j]
        class_ptEI = df.time[idx_ptEI_j]

        final_size = R[end] / params.class_sizes[j]
        class_hit  = 1.0 - S[idx_ptEI_j] / params.class_sizes[j]

        first_I_j = findfirst(x -> x > 0, I)
        last_I_j  = findlast(x -> x > 1, I)
        class_duration = df.time[last_I_j] - df.time[first_I_j]

        push!(class_metrics,
            (j, class_psI, class_ptI,
             class_psEI, class_ptEI,
             final_size, class_hit,
             class_duration)
        )
    end

    return class_metrics
end


function retrieve_popstate(sol, time_point::Float64)
    idx = findfirst(t -> t >= time_point, sol.t)
    if idx === nothing
        error("Time point $time_point exceeds the solution time span.")
    end

    u_at_time = sol.u[idx]
    S_vector = u_at_time[1, :]
    E_vector = u_at_time[2, :]
    I_vector = u_at_time[3, :]
    R_vector = u_at_time[4, :]

    return S_vector, E_vector, I_vector, R_vector
end



