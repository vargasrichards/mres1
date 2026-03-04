# Tests for confinement of infection in a population with switching and without switchinf

library(testthat)
library(MixSwitchEpi)

test_that("perfectly assortative contact and no switching -> confinement of infection to the class(es) in which seeded initially ", {
  skip_if_not_installed("odin2")

  n_activity <- 3
  N <- 3000
  class_sizes <- rep(N / n_activity, n_activity)

  # Perfect assortativity (identity mixing)
  m_identity <- diag(n_activity)

  # No switching
  swch_zero <- matrix(0, n_activity, n_activity)

  # Loop through each class to test confinement
  for (focal_class in 1:n_activity) {
    # Infect only focal_class
    ic <- make_initial_specific(n_activity = n_activity, initial_class = focal_class, fraction_exposed = 0.01)

    pars <- list(
      n_activity = n_activity,
      class_sizes = class_sizes,
      activity_scores = rep(1, n_activity), # uniform activity
      S0 = ic$S0,
      E0 = ic$E0,
      N = N,
      beta = 0.8,
      sigma = 1 / 3,
      gamma = 0.2, # R0 = 0.8/0.2 = 4
      omega = 0,
      m = m_identity,
      swch = swch_zero,
      t_end = 200
    )

    # Run simulation
    sys <- dust2::dust_system_create(
      seirs_act(),
      pars = pars,
      ode_control = dust2::dust_ode_control(max_steps = 1e6)
    )
    dust2::dust_system_set_state_initial(sys)
    t <- seq(0, pars$t_end, by = 1)
    out <- dust2::dust_system_simulate(sys, t)
    nt <- neaten_state(system_used = sys, dust_state = out, times_simulated = t)

    # Check Focal Class had an epidemic
    focal_R_col <- paste0("R", focal_class)
    final_R_focal <- nt[[focal_R_col]][nrow(nt)]
    expect_true(final_R_focal > 100, label = paste("Class", focal_class, "should have epidemic"))

    # Check ALL other classes remained uninfected
    for (other_class in setdiff(1:n_activity, focal_class)) {
      other_R_col <- paste0("R", other_class)
      other_I_col <- paste0("I", other_class)
      other_E_col <- paste0("E", other_class)

      expect_equal(nt[[other_R_col]], rep(0, nrow(nt)),
        label = paste("Class", other_class, "should be uninfected when class", focal_class, "is seeded")
      )
      expect_equal(nt[[other_I_col]], rep(0, nrow(nt)))
      expect_equal(nt[[other_E_col]], rep(0, nrow(nt)))
    }
  }
})

test_that("switching breaks confinement even with perfect assortativity", {
  skip_if_not_installed("odin2")

  n_activity <- 2
  N <- 2000
  class_sizes <- rep(N / n_activity, n_activity)

  # Perfect assortativity (identity mixing)
  m_identity <- diag(n_activity)

  # With switching!
  swch <- simple_switching(switch_rate = 0.1, num_classes = n_activity)

  # Infect only class 1
  ic <- make_initial_specific(n_activity = n_activity, initial_class = 1, fraction_exposed = 0.01)
  S0 <- ic$S0
  E0 <- ic$E0

  pars <- list(
    n_activity = n_activity,
    class_sizes = class_sizes,
    activity_scores = rep(1, n_activity),
    S0 = S0,
    E0 = E0,
    N = N,
    beta = 0.8,
    sigma = 1 / 3,
    gamma = 0.2,
    omega = 0,
    m = m_identity,
    swch = swch,
    t_end = 200
  )

  sys <- dust2::dust_system_create(
    seirs_act(),
    pars = pars,
    ode_control = dust2::dust_ode_control(max_steps = 1e6)
  )
  dust2::dust_system_set_state_initial(sys)
  t <- seq(0, pars$t_end, by = 1)
  out <- dust2::dust_system_simulate(sys, t)
  nt <- neaten_state(system_used = sys, dust_state = out, times_simulated = t)

  # Class 2 should get infected due to switching
  final_R2 <- nt$R2[nrow(nt)]
  expect_true(final_R2 > 10) # Should be substantial
})
