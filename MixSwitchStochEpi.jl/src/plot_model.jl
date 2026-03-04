using ColorSchemes

"""
    plot_compartments(sol, class_idx = nothing)

Plot occupancy of infection states through time
"""
function plot_compartments(sol, class_idx = nothing)
    t = sol.times
    states = sol.states
        S_total = [sum(states[i][1, :]) for i = 1:length(states)]
        E_total = [sum(states[i][2, :]) for i = 1:length(sol)]
        I_total = [sum(states[i][3, :]) for i = 1:length(sol)]
        R_total = [sum(states[i][4, :]) for i = 1:length(sol)]

        plot(
            t,
            S_total,
            label = "Susceptible",
            lw = 2,
            xlabel = "Time",
            ylabel = "Population",
            title = "SEIRS Model - Total Population",
            legend = :right,
        )
        plot!(t, E_total, label = "Exposed", lw = 2)
        plot!(t, I_total, label = "Infectious", lw = 2)
        plot!(t, R_total, label = "Recovered", lw = 2)

end

function plot_infectious_by_class(sol, params::SEIRSStochParams)
    n_activity = params.n_activity
    t = sol.t

    p = plot(
        xlabel = "Time",
        ylabel = "Infectious",
        title = "Infectious by Activity Class",
        legend = :topright,
    )

    for i = 1:n_activity
        I = [sol[j][3, i] for j = 1:length(sol)]
        plot!(t, I, label = "Class $i", lw = 2)
    end

    return p
end

function plot_SEIR(sol, filepath)
    t = sol.times
    U = sol.states

    S = [sum(Ui[1, :]) for Ui in U]
    E = [sum(Ui[2, :]) for Ui in U]
    I = [sum(Ui[3, :]) for Ui in U]
    R = [sum(Ui[4, :]) for Ui in U]

    plot(t, S, label = "S", lw = 2)
    plot!(t, E, label = "E", lw = 2)
    plot!(t, I, label = "I", lw = 2)
    plot!(t, R, label = "R", lw = 2)

    xlabel!("Time")
    ylabel!("Population")
    title!("SEIR dynamics")
    savefig(filepath)
end
