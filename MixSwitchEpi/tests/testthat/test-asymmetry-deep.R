library(testthat)
library(MixSwitchEpi)

test_that("Mixing Matrix Reciprocity and Asymmetry", {
  # Test that asymmetric initial conditions behave correctly with symmetric mixing flows
  n_act <- 2
  class_sizes <- c(1000, 10000) # Massive difference in class sizes
  act_scores <- c(1, 1)
  beta <- 0.1

  # Garnett mixing with zero assortativity (proportionate)
  # m[i, j] = a_j N_j / sum(a N)
  m <- garnett_mixing(2, 0, act_scores, class_sizes)

  # m[1, 2] = 10000 / 11000 = 0.90909
  # m[2, 1] = 1000 / 11000 = 0.09090
  expect_equal(m[1, 2], 10 / 11)
  expect_equal(m[2, 1], 1 / 11)

  pars <- list(
    n_activity = n_act,
    class_sizes = class_sizes,
    activity_scores = act_scores,
    m = m,
    beta = beta,
    sigma = 1, gamma = 0.1, omega = 0,
    swch = matrix(0, 2, 2),
    S0 = c(1, 1),
    E0 = c(0, 0)
  )

  sys <- dust2::dust_system_create(seirs_act(), pars)

  # Scenario: 100 people infected in Class 1 (the small class, 10% of it)
  # Scenario: 100 people infected in Class 2 (the large class, 1% of it)

  # A: Class 1 infected
  init_A <- c(900, 10000, 0, 0, 100, 0, 0, 0)
  sys$methods$set_time(sys$ptr, 0)
  dust2::dust_system_set_state(sys, init_A)
  out_A <- dust2::dust_system_simulate(sys, c(0, 1))
  up_A <- dust2::dust_unpack_state(sys, out_A)

  # B: Class 2 infected
  init_B <- c(1000, 9900, 0, 0, 0, 100, 0, 0)
  sys$methods$set_time(sys$ptr, 0)
  dust2::dust_system_set_state(sys, init_B)
  out_B <- dust2::dust_system_simulate(sys, c(0, 1))
  up_B <- dust2::dust_unpack_state(sys, out_B)

  # Force of infection from 100 people in class 1 on class 2:
  # lambda_2 = beta * a2 * m[2, 1] * I1/N1 = 0.1 * 1 * (1/11) * 100/1000 = 0.1 * 1/11 * 1/10 = 0.001 / 1.1 = 0.000909

  # Force of infection from 100 people in class 2 on class 1:
  # lambda_1 = beta * a1 * m[1, 2] * I2/N2 = 0.1 * 1 * (10/11) * 100/10000 = 0.1 * 10/11 * 1/100 = 0.01 / 11 = 0.000909

  # They should be EQUAL! (Reciprocity: flows are equal if classes have same activity)
  # delta_S1_B = -lambda_1 * S1 * dt = -0.000909 * 1000 * 1 = -0.909
  # delta_S2_A = -lambda_2 * S2 * dt = -0.000909 * 10000 * 1 = -9.09

  S1_at_1_B <- up_B$S[1, 2]
  S2_at_1_A <- up_A$S[2, 2]

  # Note: Because of exponents, it's safer to compare lambda directly
  lambda_1_B <- -(S1_at_1_B - 1000) / (1000 * 1)
  lambda_2_A <- -(S2_at_1_A - 10000) / (10000 * 1)

  expect_equal(lambda_1_B, lambda_2_A, tolerance = 1e-4)
})

test_that("R0 Invariance under Symmetric Population and Matrix", {
  # If parameters are symmetric, R0 should not change when classes are swapped
  n_act <- 2
  class_sizes <- c(1000, 2000)
  act_scores <- c(1, 5)

  pars_1 <- list(
    n_activity = n_act,
    class_sizes = class_sizes,
    activity_scores = act_scores,
    beta = 0.1, gamma = 0.1, sigma = 0.2,
    swch = matrix(0, 2, 2),
    epsilon = 0.2
  )
  r0_1 <- compute_r0(pars_1)

  # Swap them
  pars_2 <- pars_1
  pars_2$class_sizes <- rev(class_sizes)
  pars_2$activity_scores <- rev(act_scores)
  r0_2 <- compute_r0(pars_2)

  expect_equal(r0_1, r0_2, tolerance = 1e-10)
})

test_that("compute_rt properly scales with susceptible depletion in asymmetric classes", {
  n_act <- 2
  class_sizes <- c(1e4, 1e4)
  act_scores <- c(10, 1) # Class 1 is 10x more active

  pars <- list(
    n_activity = n_act,
    class_sizes = class_sizes,
    activity_scores = act_scores,
    beta = 0.1, gamma = 0.1, sigma = 0.2,
    swch = matrix(0, 2, 2),
    epsilon = 0 # Proportionate
  )

  R0 <- compute_r0(pars)

  # Case A: Class 1 (Hyper-active) is 50% depleted
  state_A <- list(s_sizes = c(5000, 10000))
  Rt_A <- compute_rt(pars, state_A)

  # Case B: Class 2 (Less active) is 50% depleted
  state_B <- list(s_sizes = c(10000, 5000))
  Rt_B <- compute_rt(pars, state_B)

  # Rt_A should be MUCH lower than R0 because the hyper-active class is gone.
  # Rt_B should be only SLIGHTLY lower than R0 because only the low-activity class is depleted.
  expect_true(Rt_A < Rt_B)
  expect_true(Rt_B < R0)
  expect_true(abs(Rt_B - Rt_A) > abs(R0 - Rt_B))

  # message(paste("Full R0:", R0))
  # message(paste("Rt (Hyper-active depleted):", Rt_A))
  # message(paste("Rt (Less active depleted):", Rt_B))
})
