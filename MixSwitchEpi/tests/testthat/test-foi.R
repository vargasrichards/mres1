library(testthat)
library(MixSwitchEpi)

# Helper to measure FOI from a dust system by simulating a tiny step
measure_foi <- function(sys, pars, dt = 1e-7) {
  n <- pars$n_activity
  # Reset time to 0 so we can simulate from 0 to dt
  sys$methods$set_time(sys$ptr, 0)
  t <- c(0, dt)
  out <- dust2::dust_system_simulate(sys, t)

  # State unpack
  s0 <- out[, 1]
  s1 <- out[, 2]

  S_at_t0 <- s0[1:n]
  S_at_t1 <- s1[1:n]

  # lambda = -delta_S / (S * dt)
  measured_lambda <- -(S_at_t1 - S_at_t0) / (pmax(S_at_t0, 1e-10) * dt)
  as.numeric(measured_lambda)
}

test_that("FOI is uniform when activity and mixing are uniform (In-model test)", {
  n_act <- 3
  N <- 1e6
  class_sizes <- rep(N / n_act, n_act)
  act_scores <- rep(1, n_act)
  beta <- 0.1

  # Uniform mixing
  rho <- garnett_mixing(n_act, 0, act_scores, class_sizes)

  pars <- list(
    n_activity = n_act,
    class_sizes = class_sizes,
    activity_scores = act_scores,
    m = rho,
    beta = beta,
    sigma = 0, gamma = 0, omega = 0,
    swch = matrix(0, n_act, n_act),
    S0 = rep(0.9, n_act),
    E0 = rep(0, n_act)
  )

  sys <- dust2::dust_system_create(seirs_act(), pars)
  # Set uniform I
  I_counts <- rep(1000, n_act)
  initial_state <- c(pars$S0 * class_sizes, rep(0, n_act), I_counts, rep(0, n_act))
  dust2::dust_system_set_state(sys, initial_state)

  foi <- measure_foi(sys, pars)

  expect_equal(sd(foi), 0, tolerance = 1e-4)
  expect_true(all(foi > 0))
})

test_that("Class with zero activity receives zero FOI (In-model test)", {
  n_act <- 2
  N <- 1e5
  class_sizes <- c(5e4, 5e4)
  act_scores <- c(0, 1) # Class 1 has NO activity
  beta <- 0.1

  rho <- garnett_mixing(n_act, 0.5, act_scores, class_sizes)

  pars <- list(
    n_activity = n_act,
    class_sizes = class_sizes,
    activity_scores = act_scores,
    m = rho,
    beta = beta,
    sigma = 0, gamma = 0, omega = 0,
    swch = matrix(0, n_act, n_act),
    S0 = c(1, 1),
    E0 = c(0, 0)
  )

  sys <- dust2::dust_system_create(seirs_act(), pars)
  # Give all infection to class 2
  I_counts <- c(0, 1000)
  initial_state <- c(class_sizes, rep(0, n_act), I_counts, rep(0, n_act))
  dust2::dust_system_set_state(sys, initial_state)

  foi <- measure_foi(sys, pars)

  # Class 1 should have 0 FOI because its activity score is 0
  expect_equal(foi[1], 0, tolerance = 1e-10)
  # Class 2 should have FOI > 0
  expect_true(foi[2] > 0)
})

test_that("FOI matches formula: lambda_i = beta * score_i * sum_j(rho_ij * I_j / N_j) (In-model test)", {
  n_act <- 3
  class_sizes <- c(10000, 20000, 30000)
  act_scores <- c(1, 2, 5)
  beta <- 0.1

  rho <- matrix(c(
    0.7, 0.2, 0.1,
    0.1, 0.8, 0.1,
    0.05, 0.05, 0.9
  ), nrow = 3, byrow = TRUE)

  pars <- list(
    n_activity = n_act,
    class_sizes = class_sizes,
    activity_scores = act_scores,
    m = rho,
    beta = beta,
    sigma = 0, gamma = 0, omega = 0,
    swch = matrix(0, n_act, n_act),
    S0 = c(0.9, 0.8, 0.7),
    E0 = c(0, 0, 0)
  )

  I_counts <- c(500, 1000, 2000)

  sys <- dust2::dust_system_create(seirs_act(), pars)
  initial_state <- c(pars$S0 * class_sizes, rep(0, n_act), I_counts, rep(0, n_act))
  dust2::dust_system_set_state(sys, initial_state)

  measured_foi <- measure_foi(sys, pars)

  expected_foi <- numeric(n_act)
  for (i in 1:n_act) {
    expected_foi[i] <- beta * act_scores[i] * sum(rho[i, ] * (I_counts / class_sizes))
  }

  expect_equal(measured_foi, expected_foi, tolerance = 1e-4)
})

test_that("FOI matches formula at an intermediate time point (post-run)", {
  n_act <- 3
  class_sizes <- c(3e5, 3e5, 4e5)
  act_scores <- c(1, 10, 20)
  beta <- 0.05
  sigma <- 0.1
  gamma <- 0.05

  rho <- garnett_mixing(n_act, 0.2, act_scores, class_sizes)
  swch <- matrix(0, n_act, n_act)

  pars <- list(
    n_activity = n_act,
    class_sizes = class_sizes,
    activity_scores = act_scores,
    m = rho,
    beta = beta,
    sigma = sigma,
    gamma = gamma,
    omega = 0,
    swch = swch,
    S0 = rep(0.999, n_act),
    E0 = rep(0.001, n_act),
    N = sum(class_sizes),
    t_end = 100
  )

  sys <- dust2::dust_system_create(seirs_act(), pars)
  dust2::dust_system_set_state_initial(sys)

  # Run the model to an intermediate point (t=20)
  out_20 <- dust2::dust_system_simulate(sys, c(0, 20))
  state_20 <- out_20[, 2]

  unpacked_20 <- dust2::dust_unpack_state(sys, state_20)
  I_20 <- as.numeric(unpacked_20$I)

  # Manually calculate FOIs using the state at t=20
  expected_foi <- numeric(n_act)
  for (i in 1:n_act) {
    expected_foi[i] <- beta * act_scores[i] * sum(rho[i, ] * (I_20 / class_sizes))
  }

  # Measure FOI from model by taking a tiny step from the state at t=20
  dust2::dust_system_set_state(sys, state_20)
  measured_foi <- measure_foi(sys, pars)

  expect_equal(measured_foi, expected_foi, tolerance = 1e-4)
})

test_that("FOI accurately captures cross-class infection (Strict Asymmetry)", {
  n_act <- 2
  class_sizes <- c(1e5, 1e5)
  act_scores <- c(1, 1)
  beta <- 0.1

  # Strictly asymmetric/directional mixing:
  # Class 1 contacts Class 2, but Class 2 contacts nobody (only itself)
  rho <- matrix(c(
    0.1, 0.9, # Class 1 spends 90% time with Class 2
    0,   1.0 # Class 2 spends 100% time with Class 2 (isolated from 1)
  ), nrow = 2, byrow = TRUE)

  pars <- list(
    n_activity = n_act,
    class_sizes = class_sizes,
    activity_scores = act_scores,
    m = rho,
    beta = beta,
    sigma = 0, gamma = 0, omega = 0,
    swch = matrix(0, n_act, n_act),
    S0 = c(1, 1),
    E0 = c(0, 0)
  )

  # Scenario A: Only Class 2 is infectious.
  I_counts_A <- c(0, 1000)
  sys <- dust2::dust_system_create(seirs_act(), pars)
  initial_state_A <- c(class_sizes, rep(0, n_act), I_counts_A, rep(0, n_act))
  dust2::dust_system_set_state(sys, initial_state_A)
  foi_A <- measure_foi(sys, pars)

  # Expected:
  # lambda_1 = beta * a1 * (m11*I1/N1 + m12*I2/N2) = 0.1 * 1 * (0.1*0 + 0.9*1000/1e5) = 0.0009
  expect_equal(foi_A[1], 0.0009, tolerance = 1e-5)
  # lambda_2 = beta * a2 * (m21*I1/N1 + m22*I2/N2) = 0.1 * 1 * (0*0 + 1.0*1000/1e5) = 0.001
  expect_equal(foi_A[2], 0.001, tolerance = 1e-5)

  # Scenario B: Only Class 1 is infectious.
  I_counts_B <- c(1000, 0)
  initial_state_B <- c(class_sizes, rep(0, n_act), I_counts_B, rep(0, n_act))
  dust2::dust_system_set_state(sys, initial_state_B)
  foi_B <- measure_foi(sys, pars)

  # Expected:
  # lambda_1 = beta * a1 * (m1,1 * I1/N1 + m1,2 * I2/N2) = 0.1 * 1 * (0.1 * 0.01 + 0.9 * 0) = 0.0001
  expect_equal(foi_B[1], 0.0001, tolerance = 1e-5)
  # lambda_2 = beta * a2 * (m2,1 * I1/N1 + m2,2 * I2/N2) = 0.1 * 1 * (0 * 0.01 + 1.0 * 0) = 0
  expect_equal(foi_B[2], 0, tolerance = 1e-10)
})

test_that("FOI scales correctly with activity_scores (asymmetric scores)", {
  n_act <- 2
  class_sizes <- rep(1e5, 2)
  beta <- 0.1
  rho <- matrix(0.5, 2, 2)
  act_scores <- c(10, 1) # Class 1 is 10x more active

  pars <- list(
    n_activity = n_act,
    class_sizes = class_sizes,
    activity_scores = act_scores,
    m = rho,
    beta = beta,
    sigma = 0, gamma = 0, omega = 0,
    swch = matrix(0, n_act, n_act),
    S0 = c(1, 1),
    E0 = c(0, 0)
  )

  I_counts <- c(1000, 1000)
  sys <- dust2::dust_system_create(seirs_act(), pars)
  initial_state <- c(class_sizes, rep(0, n_act), I_counts, rep(0, n_act))
  dust2::dust_system_set_state(sys, initial_state)

  foi <- measure_foi(sys, pars)

  expect_equal(as.numeric(foi[1] / foi[2]), 10, tolerance = 1e-3)
})
