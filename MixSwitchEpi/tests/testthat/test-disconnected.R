library(testthat)
library(MixSwitchEpi)
library(finalsize)

test_that("Totally assortative independent epidemics match single-population expectations", {
  skip_if_not_installed("dust2")
  skip_if_not_installed("odin2")

  n_act <- 3
  N <- 1000000
  class_sizes <- rep(N / n_act, n_act)

  # R0 = beta * Act / gamma
  # beta=0.05, gamma=0.1. Act: 1 (R0=0.5), 4 (R0=2.0), 10 (R0=5.0)
  act_scores <- c(1, 4, 10)

  # Assortative Mixing (epsilon = 1)
  epsilon_val <- 1.0
  m <- garnett_mixing(
    n_activity = n_act,
    epsilon_assort = epsilon_val,
    act_vector = act_scores,
    class_size_vector = class_sizes
  )

  expect_equal(m, diag(n_act))

  scalars <- basic_scalar_parms()
  beta <- 0.05
  gamma <- 0.1

  # Seed large enough to be resolved by ODE solver easily
  e0_val <- 1e-4
  s0_val <- 1 - e0_val

  pars <- list(
    n_activity = n_act,
    class_sizes = class_sizes,
    activity_scores = act_scores,
    m = m,
    swch = matrix(0, n_act, n_act),
    beta = beta,
    gamma = gamma,
    sigma = scalars$sigma,
    omega = 0,
    t_end = 500,
    S0 = rep(s0_val, n_act),
    E0 = rep(e0_val, n_act),
    N = N
  )

  # Simulation
  sys <- dust2::dust_system_create(seirs_act(), pars, ode_control = dust2::dust_ode_control(max_steps = 1e6))
  dust2::dust_system_set_state_initial(sys)
  t <- seq(0, 500, by = 1)
  out <- dust2::dust_system_simulate(sys, t)

  res <- neaten_state(sys, out, t)

  # Debug: show last row to see epidemic outcomes
  # print(tail(res, 1))

  last_row <- res[nrow(res), ]

  for (i in 1:n_act) {
    r0_i <- beta * act_scores[i] / gamma

    # Extract final size from the last row (Removed fraction)
    # Neaten state columns for arrays are named 'R1', 'R2', etc.
    col_name <- paste0("R", i)
    if (!col_name %in% names(last_row)) {
      stop(paste("Column", col_name, "missing from output. Available:", paste(names(last_row), collapse = ", ")))
    }
    simulated_fs <- last_row[[col_name]] / class_sizes[i]

    # Analytical prediction
    # Model 99.99% susceptible initially.
    expected_fs_obj <- finalsize::final_size(
      r0 = r0_i,
      contact_matrix = matrix(1),
      demography_vector = 1,
      susceptibility = matrix(c(1, 0), nrow = 1),
      p_susceptibility = matrix(c(s0_val, e0_val), nrow = 1)
    )

    expected_fs <- sum(expected_fs_obj$p_infected)

    expect_equal(simulated_fs, expected_fs,
      tolerance = 0.005,
      label = paste("Class", i, "Final Size (R0 =", r0_i, ")")
    )
  }
})
