library(testthat)
library(MixSwitchEpi)

test_that("R0 > 1 leads to growth and R0 < 1 leads to decay in the odin model", {
  # Case 1: R0 < 1
  pars_low <- list(
    n_activity = 1,
    class_sizes = 1e5,
    activity_scores = 1,
    m = matrix(1, 1, 1),
    beta = 0.05,
    sigma = 0.2,
    gamma = 0.1,
    omega = 0,
    swch = matrix(0, 1, 1),
    S0 = 1,
    E0 = 0
  )
  r0_low <- MixSwitchEpi::compute_r0(pars_low) # 0.05 / 0.1 = 0.5
  expect_lt(r0_low, 1)

  sys_low <- dust2::dust_system_create(seirs_act(), pars_low)
  initial_state_low <- c(1e5, 100, 100, 0) # S, E, I, R
  dust2::dust_system_set_state(sys_low, initial_state_low)

  out_low <- dust2::dust_system_simulate(sys_low, seq(0, 50, by = 10))
  up_low <- dust2::dust_unpack_state(sys_low, out_low)
  total_inf_low <- up_low$E + up_low$I
  expect_lt(as.numeric(total_inf_low[length(total_inf_low)]), as.numeric(total_inf_low[1]))

  # Case 2: R0 > 1
  pars_high <- pars_low
  pars_high$beta <- 0.5
  r0_high <- MixSwitchEpi::compute_r0(pars_high) # 0.5 / 0.1 = 5
  expect_gt(r0_high, 1)

  sys_high <- dust2::dust_system_create(seirs_act(), pars_high)
  dust2::dust_system_set_state(sys_high, initial_state_low)

  out_high <- dust2::dust_system_simulate(sys_high, seq(0, 5, by = 1))
  up_high <- dust2::dust_unpack_state(sys_high, out_high)
  total_inf_high <- up_high$E + up_high$I
  expect_gt(as.numeric(total_inf_high[length(total_inf_high)]), as.numeric(total_inf_high[1]))
})

test_that("Asymmetric mixing correctly spreads infection in odin model", {
  n_act <- 2
  class_sizes <- c(1e5, 1e5)
  # Strictly one-way mixing: Class 2 contacts Class 1, but Class 1 is isolated from 2
  rho <- matrix(c(
    1.0, 0.0, # Row 1: Class 1 only contacts Class 1
    1.0, 0.0 # Row 2: Class 2 only contacts Class 1 (!!! This is a very weird forced mixing)
  ), nrow = 2, byrow = TRUE)

  pars <- list(
    n_activity = n_act,
    class_sizes = class_sizes,
    activity_scores = c(1, 1),
    m = rho,
    beta = 0.5,
    sigma = 0.2,
    gamma = 0.1,
    omega = 0,
    swch = matrix(0, n_act, n_act),
    S0 = c(1, 1),
    E0 = c(0, 0)
  )

  # scenario: Only Class 2 infected.
  # Since Class 2 ONLY contacts Class 1, and Class 1 has NO infections,
  # Class 2 should NOT be able to infect anyone (even themselves).
  # Class 1 is isolated and should also not be infected.
  I_counts <- c(0, 1000)
  sys <- dust2::dust_system_create(seirs_act(), pars)
  initial_state <- c(class_sizes, rep(0, n_act), I_counts, rep(0, n_act))
  dust2::dust_system_set_state(sys, initial_state)

  out <- dust2::dust_system_simulate(sys, seq(0, 10, by = 2))
  up <- dust2::dust_unpack_state(sys, out)

  # Class 1 remained 0
  expect_equal(as.numeric(up$E[1, ]), rep(0, ncol(up$E)))
  expect_equal(as.numeric(up$I[1, ]), rep(0, ncol(up$I)))
  # Class 2 infections should decay
  expect_lt(as.numeric(up$I[2, ncol(up$I)]), as.numeric(up$I[2, 1]))
})
