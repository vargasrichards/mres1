# extension work on observation model using Turing.jl plot_paths
# https://turing.ml/dev/tutorials/04-hidden-markov-model/ was used to begin this script - thx
# Load libraries.
using Turing, StatsPlots, Random

# Turn off progress monitor.
Turing.setprogress!(false);

# Set a random seed and use the forward_diff AD mode.
Random.seed!(12345678);

using Turing
using Distributions
using LinearAlgebra
using Plots
using StatsPlots

"""
Bayesian model for inferring interpersonal proximity from binary detection observations.

The model assumes:
- True proximity d(t) varies smoothly over time
- Detection probability P(detect | d) = p_detect if d < radius, 0 otherwise
- We observe binary detection events y(t) ∈ {0, 1}
- We want to infer the latent true proximity trajectory

Parameters:
- observations: Vector of binary observations (1 = detected, 0 = not detected)
- detection_radius: Maximum distance for detection
- p_detect: Probability of detection when within radius
- dt: Time step between observations
"""

@model function proximity_inference(observations, detection_radius, p_detect; dt=1.0)
    N = length(observations)
    
    # Priors for initial conditions
    d0 ~ truncated(Normal(detection_radius/2, detection_radius/4), 0, Inf)  # Initial distance
    v0 ~ Normal(0, 0.5)  # Initial velocity (m/s or similar units)
    
    # Process noise (how much the acceleration can vary)
    σ_accel ~ truncated(Normal(0.2, 0.1), 0.01, 1.0)
    
    # Arrays to store latent states
    distances = Vector{Float64}(undef, N)
    velocities = Vector{Float64}(undef, N)
    
    # Initial state
    distances[1] = d0
    velocities[1] = v0
    
    # Dynamic model: random walk with velocity
    for t in 2:N
        # Acceleration drawn from zero-mean Gaussian
        accel ~ Normal(0, σ_accel)
        
        # Update velocity and distance
        velocities[t] = velocities[t-1] + accel * dt
        distances[t] = max(0.0, distances[t-1] + velocities[t] * dt)
    end
    
    # Observation model
    for t in 1:N
        # Detection probability based on distance
        p_obs = distances[t] < detection_radius ? p_detect : 0.0
        
        # Observe binary detection
        observations[t] ~ Bernoulli(p_obs)
    end
    
    return distances
end


"""
Alternative model with Gaussian process prior on distances for smoother inference.
"""
@model function proximity_inference_gp(observations, detection_radius, p_detect; 
                                       dt=1.0, length_scale=5.0)
    N = length(observations)
    
    # Mean distance prior
    μ_dist ~ truncated(Normal(detection_radius/2, detection_radius/3), 0, Inf)
    
    # GP hyperparameters
    σ_gp ~ truncated(Normal(detection_radius/4, detection_radius/8), 0.1, detection_radius)
    ℓ ~ truncated(Normal(length_scale, 2.0), 1.0, 20.0)  # Length scale
    
    # Construct GP covariance matrix (squared exponential kernel)
    times = collect(1:N) .* dt
    K = [σ_gp^2 * exp(-0.5 * ((times[i] - times[j])/ℓ)^2) for i in 1:N, j in 1:N]
    K += 1e-6 * I  # Add jitter for numerical stability
    
    # Sample distances from GP
    raw_distances ~ MvNormal(fill(μ_dist, N), K)
    distances = max.(0.0, raw_distances)  # Enforce non-negativity
    
    # Observation model
    for t in 1:N
        p_obs = distances[t] < detection_radius ? p_detect : 0.0
        observations[t] ~ Bernoulli(p_obs)
    end
    
    return distances
end


"""
Simplified model with discrete state (within/outside radius) and smooth transitions.
"""
@model function proximity_inference_discrete(observations, detection_radius, p_detect)
    N = length(observations)
    
    # Prior on initial state: probability of being within radius
    p_initial ~ Beta(2, 2)
    
    # Transition probabilities
    p_stay_in ~ Beta(8, 2)   # Probability of staying inside if already inside
    p_stay_out ~ Beta(8, 2)  # Probability of staying outside if already outside
    
    # Initial state
    states = Vector{Int}(undef, N)
    states[1] ~ Bernoulli(p_initial)
    
    # State transitions (Hidden Markov Model)
    for t in 2:N
        if states[t-1] == 1  # Was inside
            states[t] ~ Bernoulli(p_stay_in)
        else  # Was outside
            states[t] ~ Bernoulli(1 - p_stay_out)
        end
    end
    
    # Observations
    for t in 1:N
        if states[t] == 1  # Inside radius
            observations[t] ~ Bernoulli(p_detect)
        else  # Outside radius
            observations[t] ~ Bernoulli(0.0)
        end
    end
    
    return states
end


"""
Generate synthetic data for testing.
"""
function generate_synthetic_data(N, detection_radius, p_detect; dt=1.0)
    # True distance trajectory (sine wave with drift)
    t = collect(0:N-1) .* dt
    true_distances = detection_radius .* (0.7 .+ 0.5 .* sin.(2π .* t ./ 20) .+ 0.1 .* randn(N))
    true_distances = max.(0.0, true_distances)
    
    # Generate observations
    observations = zeros(Int, N)
    for i in 1:N
        if true_distances[i] < detection_radius
            observations[i] = rand() < p_detect ? 1 : 0
        else
            observations[i] = 0
        end
    end
    
    return observations, true_distances
end


"""
Run inference and visualize results.
"""
function run_proximity_inference_example()
    # Parameters
    N = 100
    detection_radius = 2.0  # meters
    p_detect = 0.8  # 80% detection probability when within radius
    dt = 1.0  # 1 second between observations
    
    # Generate synthetic data
    observations, true_distances = generate_synthetic_data(N, detection_radius, p_detect; dt=dt)
    
    println("Running inference with Random Walk model...")
    model_rw = proximity_inference(observations, detection_radius, p_detect; dt=dt)
    chain_rw = sample(model_rw, NUTS(), 1000)
    
    println("\nRunning inference with Discrete State model...")
    model_discrete = proximity_inference_discrete(observations, detection_radius, p_detect)
    chain_discrete = sample(model_discrete, NUTS(), 1000)
    
    # Extract posterior samples
    distances_samples = Array(group(chain_rw, :distances))
    
    # Compute posterior statistics
    distances_mean = mean(distances_samples, dims=1)[:]
    distances_lower = [quantile(distances_samples[:, i], 0.025) for i in 1:N]
    distances_upper = [quantile(distances_samples[:, i], 0.975) for i in 1:N]
    
    # Plotting
    p1 = plot(1:N, true_distances, label="True Distance", lw=2, 
              xlabel="Time", ylabel="Distance (m)", 
              title="Proximity Inference from Binary Observations")
    plot!(p1, 1:N, distances_mean, ribbon=(distances_mean .- distances_lower, 
                                           distances_upper .- distances_mean),
          label="Inferred (95% CI)", lw=2, alpha=0.3)
    hline!(p1, [detection_radius], label="Detection Radius", ls=:dash, lw=2, color=:red)
    
    # Plot observations as scatter
    detection_times = findall(observations .== 1)
    scatter!(p1, detection_times, zeros(length(detection_times)), 
             label="Detection", marker=:circle, markersize=4, color=:green)
    
    p2 = histogram(chain_rw[:σ_accel], xlabel="σ_accel", ylabel="Frequency",
                   title="Posterior: Process Noise", legend=false, bins=30)
    
    # Discrete state analysis
    states_samples = Array(group(chain_discrete, :states))
    states_prob = mean(states_samples, dims=1)[:]
    
    p3 = plot(1:N, states_prob, label="P(Within Radius)", lw=2,
              xlabel="Time", ylabel="Probability",
              title="Discrete State Model: Within Radius Probability")
    scatter!(p3, detection_times, ones(length(detection_times)), 
             label="Detection", marker=:circle, markersize=4, color=:green, alpha=0.5)
    
    plot(p1, p2, p3, layout=(3, 1), size=(800, 900))
end


# Example usage:
run_proximity_inference_example()