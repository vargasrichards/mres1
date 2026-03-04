library(testthat)
library(MixSwitchEpi)

test_that("in the limit of maximum switching the population behaves homogeneously", {
  # 1. Define Heterogeneous model with FAST switching
  n_act <- 3
  N <- 1e6
  class_sizes <- c(2e5, 3e5, 5e5)
  act_scores <- c(0.5, 1, 2) # heterogeneous

  beta_base <- 0.2
  sigma <- 0.1
  gamma <- 0.05

  # Uniform mixing (Garnett with epsilon=0)
  # This makes mixing proportional to activity
  rho <- garnett_mixing(n_act, 0, act_scores, class_sizes)

  # Fast switching (res time = 0.005 days)
  # 200 is fast enough compared to gamma=0.05 to show homogenization
  swch <- simple_switching(200, n_act)

  pars_hetero <- list(
    n_activity = n_act,
    class_sizes = class_sizes,
    activity_scores = act_scores,
    m = rho,
    beta = beta_base,
    sigma = sigma,
    gamma = gamma,
    omega = 0,
    swch = swch,
    S0 = rep(0.999, n_act),
    E0 = rep(0.001, n_act),
    N = N,
    t_end = 500
  )

  # 2. Run Heterogeneous simulation
  t_vals <- seq(0, 500, by = 1)
  # Create system with increased max_steps to avoid "too many steps"
  sys_hetero <- dust2::dust_system_create(seirs_act(), pars_hetero,
    ode_control = dust2::dust_ode_control(max_steps = 100000)
  )
  dust2::dust_system_set_state_initial(sys_hetero)
  out_hetero <- dust2::dust_system_simulate(sys_hetero, t_vals)
  res_hetero <- neaten_state(sys_hetero, out_hetero, t_vals)

  # 3. Calculate metrics for the heterogeneous model
  # Final size = 1 - S_final/N
  final_size_hetero <- 1 - (res_hetero[.N, s_tot] / res_hetero[.N, pop_tot])
  peak_I_hetero <- res_hetero[, max(i_tot)] / N

  # 4. Run Homogeneous reference
  # Check R0 computed by NGM at this switching rate.
  r0_ngm <- compute_r0(pars_hetero)

  # Simulation for homogeneous model
  times_ode <- seq(0, 500, by = 0.1)
  parameters_hom <- c(
    gamma = gamma,
    sigma = sigma,
    beta  = r0_ngm * gamma
  )
  state_hom <- c(
    S = 0.999,
    E = 0.001,
    I = 0,
    R = 0
  )
  out_hom <- deSolve::ode(
    y = state_hom, times = times_ode,
    func = seir_model, parms = parameters_hom
  )
  out_hom_df <- as.data.frame(out_hom)

  final_size_hom <- 1 - out_hom_df[nrow(out_hom_df), "S"]
  peak_I_hom <- max(out_hom_df$I)

  # 5. Compare
  # Convergence in SEIR models can be slow due to the latent period spreading,
  # but they should be reasonably close as switching rate increases.
  expect_equal(final_size_hetero, final_size_hom, tolerance = 0.05)
  expect_equal(peak_I_hetero, peak_I_hom, tolerance = 0.15)
})
