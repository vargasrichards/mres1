library(testthat)
library(MixSwitchEpi)

# ---------------------------------------------------------------------------
# Pure R0 / NGM tests (no dust2 simulation required)
# ---------------------------------------------------------------------------

test_that("equal activity gives the same R0 across assortativities", {
  n_activity <- 4
  N <- 1e6
  beta <- 0.3
  gamma <- 0.1
  no_swch <- simple_switching(switch_rate = 0, num_classes = n_activity)
  unif_act_vector <- rep(1, n_activity)
  epsilon_values <- c(0, 0.3, 0.6, 0.9)
  i <- 0
  for (epsilon in epsilon_values) {
    m <- MixSwitchEpi::generate_mixingmat(
      n_activity = n_activity,
      epsilon_assort = epsilon,
      act_vector = unif_act_vector,
      class_size_vector = rep(N / n_activity, n_activity)
    )

    pars <- list(
      m = m,
      n_activity = n_activity,
      beta = beta,
      sigma = 0.1,
      gamma = gamma,
      swch = no_swch,
      class_sizes = rep(N / n_activity, n_activity),
      activity_scores = unif_act_vector
    )

    r0 <- MixSwitchEpi::compute_r0(pars)
    if (i == 0) {
      # Store the first R0 as reference
      r0_ref <- r0
    } else {
      # For subsequent epsilon values, check that R0 matches the reference
      expect_equal(r0, r0_ref, tolerance = 1e-6)
    }
    i <- i + 1
  }
})

test_that("beta = 0 gives R0 = 0 and Rt = 0 for any susceptible state", {
  n_activity <- 4
  pars <- list(
    n_activity = n_activity,
    class_sizes = rep(1000, n_activity),
    activity_scores = rep(1, n_activity),
    beta = 0,
    gamma = 1 / 7,
    swch = matrix(0, n_activity, n_activity),
    sigma = 1 / 5
  )

  r0 <- MixSwitchEpi::compute_r0(pars)
  expect_equal(r0, 0)

  popstate_counts <- list(s_sizes = pars$class_sizes)
  popstate_frac <- list(s_sizes = rep(1, n_activity))

  rt_counts <- MixSwitchEpi::compute_rt(pars, popstate_counts)
  rt_frac <- MixSwitchEpi::compute_rt(pars, popstate_frac)

  expect_equal(rt_counts, 0)
  expect_equal(rt_frac, 0)
})

test_that("compute_r0 reflects threshold behaviour (R0 < 1 vs R0 > 1)", {
  n_activity <- 3
  pars_base <- list(
    n_activity = n_activity,
    class_sizes = rep(100, n_activity),
    activity_scores = rep(1, n_activity),
    beta = 0.05,
    gamma = 0.1,
    sigma = 1,
    swch = matrix(0, n_activity, n_activity)
  )

  r0_base <- MixSwitchEpi::compute_r0(pars_base)
  expect_true(r0_base < 1)

  pars_high <- pars_base
  pars_high$beta <- 0.5
  r0_high <- MixSwitchEpi::compute_r0(pars_high)
  expect_true(r0_high > 1)
})


# ---------------------------------------------------------------------------
# dust2 simulation tests
# ---------------------------------------------------------------------------

# Helper: create system with generous ODE control and zero switching by default
create_sys <- function(params) {
  sys <- dust2::dust_system_create(
    seirs_act(),
    pars = params,
    ode_control = dust2::dust_ode_control(max_steps = 1e6)
  )
  dust2::dust_system_set_state_initial(sys)
  sys
}

test_that("SEIRS activity model builds and initialises correctly", {
  skip_if_not_installed("odin2")

  n_activity <- 3

  params <- list(
    n_activity = n_activity,
    S0 = rep(0.9 / n_activity, n_activity),
    E0 = rep(0.1 / n_activity, n_activity),
    N = 1000,
    beta = 0.5,
    sigma = 1 / 5,
    gamma = 1 / 7,
    omega = 0.0,
    m = matrix(1, n_activity, n_activity),
    swch = matrix(0, n_activity, n_activity),
    class_sizes = rep(1000 / n_activity, n_activity),
    activity_scores = rep(1, n_activity)
  )

  sys <- create_sys(params)

  y0 <- dust2::dust_system_state(sys)

  expect_length(y0, 4 * n_activity) # S,E,I,R
  expect_true(all(y0[1:n_activity] > 0)) # S
  expect_true(all(y0[(n_activity + 1):(2 * n_activity)] > 0)) # E
  expect_true(all(y0[(2 * n_activity + 1):(3 * n_activity)] == 0)) # I
  expect_true(all(y0[(3 * n_activity + 1):(4 * n_activity)] == 0)) # R
})

test_that("Population is conserved when omega = 0 and no switching", {
  skip_if_not_installed("odin2")

  n_activity <- 3
  N <- 1200
  params <- list(
    n_activity = n_activity,
    S0 = rep(0.9, n_activity),
    E0 = rep(0.1, n_activity),
    N = N,
    beta = 0.5,
    sigma = 1 / 5,
    gamma = 1 / 7,
    omega = 0,
    m = matrix(1, n_activity, n_activity),
    swch = matrix(0, n_activity, n_activity),
    class_sizes = rep(400, n_activity),
    activity_scores = rep(1, n_activity)
  )

  sys <- create_sys(params)
  t <- seq(0, 50, by = 1)
  out <- dust2::dust_system_simulate(sys, t)
  up <- dust2::dust_unpack_state(sys, out)

  pop <- colSums(up$S) + colSums(up$E) + colSums(up$I) + colSums(up$R)
  expect_equal(pop, rep(N, length(t)), tolerance = 1e-6)
})

test_that("Model runs forward without NaN or Inf", {
  skip_if_not_installed("odin2")

  set.seed(42)
  n_activity <- 4
  # Use a proper rate matrix for switching (rows sum to 0)
  swch <- simple_switching(switch_rate = 0.05, num_classes = n_activity)

  params <- list(
    n_activity = n_activity,
    S0 = rep(0.98 / n_activity, n_activity),
    E0 = rep(0.02 / n_activity, n_activity),
    N = 10000,
    beta = 0.8,
    sigma = 1 / 4,
    gamma = 1 / 6,
    omega = 1 / 200,
    m = matrix(runif(n_activity^2, 0, 1), n_activity, n_activity),
    swch = swch,
    class_sizes = rep(2500, n_activity),
    activity_scores = runif(n_activity, 0.5, 2)
  )

  sys <- create_sys(params)
  t <- seq(0, 200, by = 1)
  out <- dust2::dust_system_simulate(sys, t)

  expect_false(any(is.na(out)))
  expect_false(any(!is.finite(out)))
})

test_that("Mass conservation holds with switching", {
  skip_if_not_installed("odin2")

  set.seed(42)
  n_activity <- 5
  N <- 2000
  swch <- simple_switching(switch_rate = 0.1, num_classes = n_activity)

  params <- list(
    n_activity = n_activity,
    S0 = rep(0.95, n_activity),
    E0 = rep(0.05, n_activity),
    N = N,
    beta = 0.3,
    sigma = 1 / 5,
    gamma = 1 / 7,
    omega = 0,
    m = matrix(runif(n_activity^2, 0, 1), n_activity, n_activity),
    swch = swch,
    class_sizes = rep(N / n_activity, n_activity),
    activity_scores = rep(1, n_activity)
  )

  sys <- create_sys(params)
  t <- seq(0, 100, by = 1)
  out <- dust2::dust_system_simulate(sys, t)
  up <- dust2::dust_unpack_state(sys, out)

  pop <- colSums(up$S) + colSums(up$E) + colSums(up$I) + colSums(up$R)
  expect_equal(pop, rep(N, length(t)), tolerance = 1e-6)
})

test_that("Uniform mixing + symmetric initial conditions produce identical trajectories", {
  skip_if_not_installed("odin2")

  n_activity <- 4
  N <- 1000

  m <- matrix(1, n_activity, n_activity)
  class_sizes <- rep(N / n_activity, n_activity)
  activity_scores <- rep(1, n_activity)

  S0 <- rep(0.9 / n_activity, n_activity)
  E0 <- rep(0.1 / n_activity, n_activity)

  params <- list(
    n_activity = n_activity,
    S0 = S0,
    E0 = E0,
    N = N,
    beta = 0.4,
    sigma = 1 / 4,
    gamma = 1 / 6,
    omega = 0.0,
    m = m,
    swch = matrix(0, n_activity, n_activity),
    class_sizes = class_sizes,
    activity_scores = activity_scores
  )

  sys <- create_sys(params)
  t <- seq(0, 80, by = 2)
  out <- dust2::dust_system_simulate(sys, t)
  up <- dust2::dust_unpack_state(sys, out)

  for (cmp in c("S", "E", "I", "R")) {
    mat <- up[[cmp]] # n_activity x n_times
    ref <- mat[1, ]
    for (k in 2:n_activity) {
      expect_equal(as.numeric(mat[k, ]), as.numeric(ref), tolerance = 1e-8)
    }
  }
})

test_that("Proportionate mixing and no switching gives expected final size", {
  skip_if_not_installed("odin2")

  N <- 1e5
  n_activity <- 3
  grp_size <- N / n_activity
  m_default <- matrix(1 / n_activity,
    nrow = n_activity,
    ncol = n_activity
  )
  default_parms <- list(
    m = m_default, n_activity = n_activity,
    S0 = rep(0.9 / n_activity, n_activity),
    E0 = rep(0.1 / n_activity, n_activity),
    beta = 0.3,
    gamma = 0.1,
    sigma = 0.2,
    N = N,
    omega = 0,
    swch = matrix(0, n_activity, n_activity),
    class_sizes = rep(grp_size, n_activity),
    activity_scores = rep(1, n_activity)
  )

  t <- seq(0, 500, by = 1)
  sys <- create_sys(default_parms)
  out <- dust2::dust_system_simulate(sys, t)
  up <- dust2::dust_unpack_state(sys, out)

  final_removeds <- sum(up$R[, ncol(up$R)])

  R0 <- default_parms$beta / default_parms$gamma
  final_size_analytical <- N * (1 - exp(-R0 * final_removeds / N))
  expect_equal(final_removeds,
    final_size_analytical,
    tolerance = 0.05 * N
  )
})

test_that("switching redistributes hosts", {
  skip_if_not_installed("odin2")

  N <- 1e5
  n_activity <- 5
  class_sizes <- rep(N / n_activity, n_activity)

  cmat <- generate_mixingmat(
    n_activity = n_activity,
    epsilon_assort = 0.1,
    act_vector = 1:n_activity,
    class_size_vector = class_sizes
  )

  some_switching <- simple_switching(
    switch_rate = 0.1,
    num_classes = n_activity
  )

  test_switch_parms <- list(
    m = cmat, n_activity = n_activity,
    S0 = c(0.99, 0, 0, 0, 0),
    E0 = c(0.01, 0, 0, 0, 0),
    beta = 0.3,
    gamma = 0.1,
    sigma = 0.2,
    N = N,
    omega = 0,
    swch = some_switching,
    class_sizes = class_sizes,
    activity_scores = 1:n_activity
  )

  t <- seq(0, 500, by = 1)
  sys <- create_sys(test_switch_parms)
  out <- dust2::dust_system_simulate(sys, t)
  nt <- neaten_state(system_used = sys, dust_state = out, times_simulated = t)

  final_class_sizes <- as.numeric(
    nt[nrow(nt), paste0("class", 1:n_activity), with = FALSE]
  )

  expected_class_sizes <- rep(N / n_activity, n_activity)
  expect_equal(final_class_sizes,
    expected_class_sizes,
    tolerance = 0.05 * N
  )
})

test_that("Asymmetric switching leads to correct stationary distribution", {
  skip_if_not_installed("odin2")

  N <- 1e5
  n_activity <- 2
  # Rate matrix: Rows sum to 0.
  # Row 1: 1 -> 2 at rate 0.4. Row 2: 2 -> 1 at rate 0.1.
  swch <- matrix(c(
    -0.4, 0.4,
    0.1, -0.1
  ), nrow = 2, byrow = TRUE)

  # Expected stationary: N1 * 0.4 = N2 * 0.1 => N1 = 0.25 * N2.
  # N1 + N2 = 1.25 N2 = 100,000 => N2 = 80,000, N1 = 20,000.

  params <- list(
    n_activity = n_activity,
    S0 = c(1, 1), # Start with equal fractions
    E0 = c(0, 0),
    N = N,
    beta = 0, # No epidemic
    sigma = 0.2,
    gamma = 0.1,
    omega = 0,
    m = matrix(0.5, 2, 2),
    swch = swch,
    class_sizes = c(5e4, 5e4), # Start with equal sizes
    activity_scores = c(1, 1)
  )

  sys <- create_sys(params)
  # Run long enough to reach equilibrium
  t <- seq(0, 100, by = 10)
  out <- dust2::dust_system_simulate(sys, t)
  up <- dust2::dust_unpack_state(sys, out)

  final_N <- up$S[, ncol(up$S)] + up$E[, ncol(up$E)] + up$I[, ncol(up$I)] + up$R[, ncol(up$R)]

  expect_equal(as.numeric(final_N), c(20000, 80000), tolerance = 1e-2 * N)
})

test_that("Infection spread is limited by asymmetric mixing", {
  skip_if_not_installed("odin2")

  # Class 1 is 'Protected' (contacts only Class 1)
  # Class 2 is 'General' (contacts Class 1 and Class 2)
  # Class 3 is 'Source' (only contacted by Class 2, but Class 3 contacts Class 2)

  n_activity <- 3
  class_sizes <- rep(1e4, 3)

  # m[i, j] is contacts of i with j.
  m <- matrix(c(
    1.0, 0.0, 0.0, # 1 only sees 1
    0.5, 0.5, 0.0, # 2 sees 1 and 2
    0.0, 1.0, 0.0 # 3 only sees 2
  ), nrow = 3, byrow = TRUE)

  params <- list(
    n_activity = n_activity,
    S0 = c(1, 1, 1),
    E0 = c(0, 0, 0),
    N = 3e4,
    beta = 0.5,
    sigma = 1, # Fast transition to I
    gamma = 0.01, # Slow recovery
    omega = 0,
    m = m,
    swch = matrix(0, 3, 3),
    class_sizes = class_sizes,
    activity_scores = c(1, 1, 1)
  )

  sys <- create_sys(params)
  init_state <- c(
    9900, 10000, 10000, # S
    0, 0, 0, # E
    100, 0, 0, # I
    0, 0, 0
  ) # R
  dust2::dust_system_set_state(sys, init_state)

  t <- seq(0, 20, by = 2)
  out <- dust2::dust_system_simulate(sys, t)
  up <- dust2::dust_unpack_state(sys, out)

  # 1. Class 2 should show growth (it contacts Class 1 which is infected)
  expect_true(up$I[2, 11] > up$I[2, 1])
  # 2. Class 3 should show growth LATER (it contacts Class 2)
  expect_true(up$I[3, 11] > 0)

  # Re-verify Scenario: Only Class 3 infected.
  # Class 1 and 2 should stay 0.
  I_counts_rev <- c(0, 0, 100)
  S_counts_rev <- c(10000, 10000, 9900)
  init_state_rev <- c(S_counts_rev, 0, 0, 0, I_counts_rev, 0, 0, 0)
  sys$methods$set_time(sys$ptr, 0)
  dust2::dust_system_set_state(sys, init_state_rev)
  out_rev <- dust2::dust_system_simulate(sys, t)
  up_rev <- dust2::dust_unpack_state(sys, out_rev)

  expect_equal(as.numeric(up_rev$I[1, ]), rep(0, length(t)))
  expect_equal(as.numeric(up_rev$I[2, ]), rep(0, length(t)))
  expect_true(up_rev$I[3, 1] > 0)
})
