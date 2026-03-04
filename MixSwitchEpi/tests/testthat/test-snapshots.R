library(testthat)
library(MixSwitchEpi)

test_that("basic scalar parameters are stable", {
  expect_snapshot(basic_scalar_parms())
})

test_that("computed NGM and R0 are consistent", {
  scalars <- basic_scalar_parms()
  pars <- list(
    n_activity = 3,
    class_sizes = c(333333, 333333, 333334),
    activity_scores = c(1, 10, 20),
    beta = scalars$beta,
    gamma = scalars$gamma,
    sigma = scalars$sigma,
    m = diag(3), # Identity mixing for simplicity
    swch = matrix(0, 3, 3),
    S0 = c(0.9, 0.9, 0.9),
    E0 = c(0.1, 0.1, 0.1),
    N = 1000000,
    omega = 0,
    t_end = 100
  )
  pars$m <- garnett_mixing(
    n_activity = pars$n_activity,
    epsilon_assort = 0.5,
    act_vector = pars$activity_scores,
    class_size_vector = pars$class_sizes
  )

  expect_snapshot(compute_r0(pars))
  expect_snapshot(compute_ngm(pars))
})

test_that("neaten_state output structure is consistent", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("dust2")

  scalars <- basic_scalar_parms()
  pars <- list(
    n_activity = 2,
    class_sizes = c(500, 500),
    activity_scores = c(1, 10),
    beta = scalars$beta,
    gamma = scalars$gamma,
    sigma = scalars$sigma,
    m = diag(2),
    swch = matrix(0, 2, 2),
    S0 = c(0.99, 1),
    E0 = c(0.01, 0),
    N = 1000,
    omega = 0,
    t_end = 10
  )

  sys <- dust2::dust_system_create(seirs_act(), pars, ode_control = dust2::dust_ode_control(max_steps = 1000))
  dust2::dust_system_set_state_initial(sys)
  t <- seq(0, 10, by = 1)
  out <- dust2::dust_system_simulate(sys, t)
  nt <- neaten_state(sys, out, t)

  # Snapshot the first few rows and column names
  expect_snapshot(names(nt))
  expect_snapshot(head(nt))
})
