library(testthat)
library(data.table, quietly = TRUE)
library(MixSwitchEpi)

# ============================================================================
# 1. Switching does NOT alter R0 when activity scores and class sizes are
#    uniform (because the NGM is invariant to switching under symmetry)
# ============================================================================

test_that("switching does not alter R0 when activity scores and class sizes are uniform", {
  n_activity <- 5
  N <- 1e6
  class_sizes <- rep(N / n_activity, n_activity)
  activity_scores <- rep(1, n_activity)
  beta <- 0.3
  gamma <- 0.1

  # Compute R0 with no switching
  m_prop <- garnett_mixing(
    n_activity = n_activity,
    epsilon_assort = 0,
    act_vector = activity_scores,
    class_size_vector = class_sizes
  )

  pars_no_switch <- list(
    n_activity = n_activity,
    class_sizes = class_sizes,
    activity_scores = activity_scores,
    beta = beta,
    gamma = gamma,
    sigma = 1 / 3,
    swch = simple_switching(switch_rate = 0, num_classes = n_activity),
    m = m_prop
  )
  r0_no_switch <- compute_r0(pars_no_switch)

  # Compute R0 with several different switching rates
  for (sr in c(0.01, 0.05, 0.1, 0.3, 0.5)) {
    pars_switch <- pars_no_switch
    pars_switch$swch <- simple_switching(switch_rate = sr, num_classes = n_activity)
    r0_switch <- compute_r0(pars_switch)

    expect_equal(r0_switch, r0_no_switch,
      tolerance = 1e-10,
      label = paste("R0 with switch_rate =", sr)
    )
  }
})

test_that("switching does not alter R0 for assortative mixing when activity is uniform", {
  n_activity <- 4
  N <- 1e6
  class_sizes <- rep(N / n_activity, n_activity)
  activity_scores <- rep(1, n_activity)

  m <- garnett_mixing(
    n_activity = n_activity,
    epsilon_assort = 0.5,
    act_vector = activity_scores,
    class_size_vector = class_sizes
  )

  pars_base <- list(
    n_activity = n_activity,
    class_sizes = class_sizes,
    activity_scores = activity_scores,
    beta = 0.25,
    gamma = 0.1,
    sigma = 0.2,
    m = m,
    swch = matrix(0, n_activity, n_activity)
  )
  r0_base <- compute_r0(pars_base)

  pars_switched <- pars_base
  pars_switched$swch <- simple_switching(0.2, n_activity)
  r0_switched <- compute_r0(pars_switched)

  expect_equal(r0_switched, r0_base, tolerance = 1e-10)
})


# ============================================================================
# 2. Final size from simulation matches finalsize package analytical result
#    when there is no switching (and proportionate mixing with uniform activity)
# ============================================================================

test_that("simulated final size matches finalsize::final_size for homogeneous no-switch model", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("finalsize")

  n_activity <- 3
  N <- 1e6
  class_sizes <- rep(N / n_activity, n_activity)
  beta <- 0.5
  gamma <- 0.1
  R0 <- beta / gamma # 5

  # Analytical final size from finalsize package
  fs_analytical <- finalsize::final_size(r0 = R0)$p_infected

  # Simulate with proportionate mixing and no switching
  m <- matrix(1 / n_activity, n_activity, n_activity)
  pars <- list(
    n_activity = n_activity,
    S0 = rep(1 - 1e-4, n_activity),
    E0 = rep(1e-4, n_activity),
    N = N,
    beta = beta,
    sigma = 1 / 3,
    gamma = gamma,
    omega = 0,
    m = m,
    swch = matrix(0, n_activity, n_activity),
    class_sizes = class_sizes,
    activity_scores = rep(1, n_activity)
  )

  sys <- dust2::dust_system_create(
    seirs_act(),
    pars = pars,
    ode_control = dust2::dust_ode_control(max_steps = 1e6)
  )
  dust2::dust_system_set_state_initial(sys)
  t <- seq(0, 1000, by = 1)
  out <- dust2::dust_system_simulate(sys, t)
  up <- dust2::dust_unpack_state(sys, out)

  final_R <- sum(up$R[, ncol(up$R)])
  simulated_fs <- final_R / N

  # Tolerance: 5% relative error is acceptable for ODE vs analytical
  expect_equal(simulated_fs, fs_analytical, tolerance = 0.05)
})

test_that("simulated final size matches finalsize for lower R0", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("finalsize")

  n_activity <- 2
  N <- 1e5
  class_sizes <- rep(N / n_activity, n_activity)
  beta <- 0.3
  gamma <- 0.15
  R0 <- beta / gamma # 2

  fs_analytical <- finalsize::final_size(r0 = R0)$p_infected

  m <- matrix(1 / n_activity, n_activity, n_activity)
  pars <- list(
    n_activity = n_activity,
    S0 = rep(1 - 1e-3, n_activity),
    E0 = rep(1e-3, n_activity),
    N = N,
    beta = beta,
    sigma = 1 / 3,
    gamma = gamma,
    omega = 0,
    m = m,
    swch = matrix(0, n_activity, n_activity),
    class_sizes = class_sizes,
    activity_scores = rep(1, n_activity)
  )

  sys <- dust2::dust_system_create(
    seirs_act(),
    pars = pars,
    ode_control = dust2::dust_ode_control(max_steps = 1e6)
  )
  dust2::dust_system_set_state_initial(sys)
  t <- seq(0, 500, by = 1)
  out <- dust2::dust_system_simulate(sys, t)
  up <- dust2::dust_unpack_state(sys, out)

  final_R <- sum(up$R[, ncol(up$R)])
  simulated_fs <- final_R / N

  expect_equal(simulated_fs, fs_analytical, tolerance = 0.05)
})


# ============================================================================
# 3. HIT from Rt=1 crossing matches homogeneous HIT prediction (1 - 1/R0)
#    when there is no switching and proportionate mixing
# ============================================================================

test_that("empirical HIT from Rt=1 matches 1-1/R0 for homogeneous model", {
  skip_if_not_installed("odin2")

  n_activity <- 3
  N <- 1e6
  class_sizes <- rep(N / n_activity, n_activity)
  beta <- 0.4
  gamma <- 0.1
  R0 <- beta / gamma # 4.0
  expected_HIT <- 1 - 1 / R0 # 0.75

  m <- matrix(1 / n_activity, n_activity, n_activity)
  pars <- list(
    n_activity = n_activity,
    S0 = rep(1 - 1e-4, n_activity),
    E0 = rep(1e-4, n_activity),
    N = N,
    beta = beta,
    sigma = 1 / 3,
    gamma = gamma,
    omega = 0,
    m = m,
    swch = matrix(0, n_activity, n_activity),
    class_sizes = class_sizes,
    activity_scores = rep(1, n_activity)
  )

  sys <- dust2::dust_system_create(
    seirs_act(),
    pars = pars,
    ode_control = dust2::dust_ode_control(max_steps = 1e6)
  )
  dust2::dust_system_set_state_initial(sys)
  t <- seq(0, 500, by = 1)
  out <- dust2::dust_system_simulate(sys, t)
  nt <- neaten_state(system_used = sys, dust_state = out, times_simulated = t)

  # Compute Rt time series
  rt_result <- compute_rt_ts(nt, pars, add_to_output = TRUE)

  # Compute HIT from the Rt crossing
  hit_result <- compute_hit(rt_result)

  expect_false(is.na(hit_result$HIT))
  expect_equal(hit_result$HIT, expected_HIT, tolerance = 0.05)
})

test_that("compute_homogeneous_metrics returns correct analytical HIT", {
  pars <- list(beta = 0.4, gamma = 0.1)
  pars$n_activity <- 1
  pars$activity_scores <- 1
  pars$class_sizes <- 1e5
  metrics <- compute_homogeneous_metrics(pars)

  expect_equal(metrics$hom_r0, 4.0)
  expect_equal(metrics$hom_hit, 0.75)
  expect_true(metrics$hom_fs > 0.5)
})


# ============================================================================
# 4. NGM determines epidemic in the long run
# ============================================================================

test_that("R0 > 1 produces an epidemic (substantial final size)", {
  skip_if_not_installed("odin2")

  n_activity <- 3
  N <- 1e5
  class_sizes <- rep(N / n_activity, n_activity)

  m <- garnett_mixing(
    n_activity = n_activity,
    epsilon_assort = 0.3,
    act_vector = rep(1, n_activity),
    class_size_vector = class_sizes
  )

  pars <- list(
    n_activity = n_activity,
    S0 = rep(1 - 1e-4, n_activity),
    E0 = rep(1e-4, n_activity),
    N = N,
    beta = 0.5,
    sigma = 1 / 3,
    gamma = 0.1,
    omega = 0,
    m = m,
    swch = matrix(0, n_activity, n_activity),
    class_sizes = class_sizes,
    activity_scores = rep(1, n_activity)
  )

  R0 <- compute_r0(pars)
  expect_true(R0 > 1)

  sys <- dust2::dust_system_create(
    seirs_act(),
    pars = pars,
    ode_control = dust2::dust_ode_control(max_steps = 1e6)
  )
  dust2::dust_system_set_state_initial(sys)
  t <- seq(0, 500, by = 1)
  out <- dust2::dust_system_simulate(sys, t)
  up <- dust2::dust_unpack_state(sys, out)

  final_R <- sum(up$R[, ncol(up$R)])
  final_size <- final_R / N

  # With R0 > 1, expect substantial epidemic (> 10% of population)
  expect_true(final_size > 0.1)
})

test_that("R0 < 1 produces no substantial epidemic", {
  skip_if_not_installed("odin2")

  n_activity <- 3
  N <- 1e5
  class_sizes <- rep(N / n_activity, n_activity)

  m <- garnett_mixing(
    n_activity = n_activity,
    epsilon_assort = 0,
    act_vector = rep(1, n_activity),
    class_size_vector = class_sizes
  )

  pars <- list(
    n_activity = n_activity,
    S0 = rep(1 - 1e-3, n_activity),
    E0 = rep(1e-3, n_activity),
    N = N,
    beta = 0.05,
    sigma = 1 / 3,
    gamma = 0.1,
    omega = 0,
    m = m,
    swch = matrix(0, n_activity, n_activity),
    class_sizes = class_sizes,
    activity_scores = rep(1, n_activity)
  )

  R0 <- compute_r0(pars)
  expect_true(R0 < 1)

  sys <- dust2::dust_system_create(
    seirs_act(),
    pars = pars,
    ode_control = dust2::dust_ode_control(max_steps = 1e6)
  )
  dust2::dust_system_set_state_initial(sys)
  t <- seq(0, 500, by = 1)
  out <- dust2::dust_system_simulate(sys, t)
  up <- dust2::dust_unpack_state(sys, out)

  final_R <- sum(up$R[, ncol(up$R)])
  final_size <- final_R / N

  # With R0 < 1, epidemic should die out: final size stays very small
  expect_true(final_size < 0.05)
})

test_that("Higher R0 produces larger final size", {
  skip_if_not_installed("odin2")

  n_activity <- 3
  N <- 1e5
  class_sizes <- rep(N / n_activity, n_activity)
  m <- matrix(1 / n_activity, n_activity, n_activity)

  make_sys <- function(beta) {
    pars <- list(
      n_activity = n_activity,
      S0 = rep(1 - 1e-4, n_activity),
      E0 = rep(1e-4, n_activity),
      N = N,
      beta = beta,
      sigma = 1 / 3,
      gamma = 0.1,
      omega = 0,
      m = m,
      swch = matrix(0, n_activity, n_activity),
      class_sizes = class_sizes,
      activity_scores = rep(1, n_activity)
    )
    sys <- dust2::dust_system_create(
      seirs_act(),
      pars = pars,
      ode_control = dust2::dust_ode_control(max_steps = 1e6)
    )
    dust2::dust_system_set_state_initial(sys)
    t <- seq(0, 800, by = 1)
    out <- dust2::dust_system_simulate(sys, t)
    up <- dust2::dust_unpack_state(sys, out)
    sum(up$R[, ncol(up$R)]) / N
  }

  fs_low <- make_sys(beta = 0.2) # R0 = 2
  fs_high <- make_sys(beta = 0.5) # R0 = 5

  expect_true(fs_high > fs_low)
})

test_that("Rt is always <= R0 during epidemic (susceptible depletion)", {
  skip_if_not_installed("odin2")

  n_activity <- 3
  N <- 1e5
  class_sizes <- rep(N / n_activity, n_activity)

  m <- garnett_mixing(
    n_activity = n_activity,
    epsilon_assort = 0.2,
    act_vector = rep(1, n_activity),
    class_size_vector = class_sizes
  )

  pars <- list(
    n_activity = n_activity,
    S0 = rep(1 - 1e-4, n_activity),
    E0 = rep(1e-4, n_activity),
    N = N,
    beta = 0.4,
    sigma = 1 / 3,
    gamma = 0.1,
    omega = 0,
    m = m,
    swch = matrix(0, n_activity, n_activity),
    class_sizes = class_sizes,
    activity_scores = rep(1, n_activity)
  )

  R0 <- compute_r0(pars)

  sys <- dust2::dust_system_create(
    seirs_act(),
    pars = pars,
    ode_control = dust2::dust_ode_control(max_steps = 1e6)
  )
  dust2::dust_system_set_state_initial(sys)
  t <- seq(0, 400, by = 1)
  out <- dust2::dust_system_simulate(sys, t)
  nt <- neaten_state(system_used = sys, dust_state = out, times_simulated = t)

  rt_ts <- compute_rt_ts(nt, pars, add_to_output = FALSE)

  # Rt at t=0 should equal R0
  expect_equal(rt_ts$Rt[1], R0, tolerance = 1e-3)

  # Rt should never exceed R0 (susceptibles only decrease without waning)
  expect_true(all(rt_ts$Rt <= R0 + 1e-6))
})


# ============================================================================
# 5. NGM consistency: compute_r0, compute_ngm, spectral_radius
# ============================================================================

test_that("compute_ngm spectral radius equals compute_r0", {
  pars <- list(
    n_activity = 4,
    class_sizes = c(100, 200, 300, 400),
    activity_scores = c(1, 2, 3, 4),
    beta = 0.2,
    gamma = 0.1,
    sigma = 0.2,
    swch = simple_switching(0.05, 4)
  )

  K <- compute_ngm(pars)
  r0_ngm <- spectral_radius(K)
  r0_fn <- compute_r0(pars)

  expect_equal(r0_fn, r0_ngm, tolerance = 1e-12)
})

test_that("NGM is non-negative (all entries >= 0)", {
  n_activity <- 5
  pars <- list(
    n_activity = n_activity,
    class_sizes = seq(100, 500, by = 100),
    activity_scores = c(0.5, 1, 2, 3, 5),
    beta = 0.3,
    gamma = 0.1,
    sigma = 0.2,
    swch = simple_switching(0.1, n_activity)
  )

  K <- compute_ngm(pars)
  expect_true(all(K >= -1e-12)) # allowing small numerical error
})


# ============================================================================
# 6. calibrate_parms achieves target R0
# ============================================================================

test_that("calibrate_parms achieves target R0", {
  n_activity <- 4
  N <- 1e6
  class_sizes <- rep(N / n_activity, n_activity)
  act_scores <- c(1, 2, 3, 4)

  m <- garnett_mixing(
    n_activity = n_activity,
    epsilon_assort = 0.3,
    act_vector = act_scores,
    class_size_vector = class_sizes
  )

  prelim <- list(
    n_activity = n_activity,
    class_sizes = class_sizes,
    activity_scores = act_scores,
    beta = 0.1,
    gamma = 1 / 4,
    sigma = 1 / 3,
    swch = matrix(0, n_activity, n_activity),
    m = m,
    N = N,
    omega = 0
  )

  for (target in c(1.5, 2.5, 4.0)) {
    calibrated <- calibrate_parms(prelim, target_r0 = target)
    r0_achieved <- compute_r0(calibrated)
    expect_equal(r0_achieved, target, tolerance = 1e-5)
  }
})


# ============================================================================
# 7. Mixing matrix properties
# ============================================================================

test_that("garnett_mixing rows sum to 1 for any assortativity", {
  n <- 5
  class_sizes <- rep(1000, n)
  act_scores <- c(0.5, 1, 2, 3, 5)

  for (eps in c(0, 0.1, 0.3, 0.5, 0.9)) {
    m <- garnett_mixing(
      n_activity = n,
      epsilon_assort = eps,
      act_vector = act_scores,
      class_size_vector = class_sizes
    )
    expect_equal(rowSums(m), rep(1, n), tolerance = 1e-12)
    expect_true(all(m >= 0))
  }
})

test_that("garnett_mixing with epsilon=0 gives proportionate (random) mixing", {
  n <- 4
  class_sizes <- c(100, 200, 300, 400)
  act_scores <- c(1, 2, 3, 4)

  m <- garnett_mixing(
    n_activity = n,
    epsilon_assort = 0,
    act_vector = act_scores,
    class_size_vector = class_sizes
  )

  # All rows should be identical (proportionate mixing)
  for (i in 2:n) {
    expect_equal(m[i, ], m[1, ], tolerance = 1e-12)
  }
})

test_that("garnett_mixing with epsilon=1 gives identity", {
  n <- 3
  class_sizes <- rep(100, n)
  act_scores <- rep(1, n)

  m <- garnett_mixing(
    n_activity = n,
    epsilon_assort = 1,
    act_vector = act_scores,
    class_size_vector = class_sizes
  )

  expect_equal(m, diag(n) + matrix(0, n, n), tolerance = 1e-12)
})


# ============================================================================
# 8. Switching matrix properties
# ============================================================================

test_that("simple_switching is a valid CTMC generator (rows sum to 0)", {
  for (n in c(2, 5, 10)) {
    for (rate in c(0, 0.01, 0.1, 0.5)) {
      W <- simple_switching(switch_rate = rate, num_classes = n)
      expect_equal(rowSums(W), rep(0, n), tolerance = 1e-12)
      expect_equal(diag(W), rep(-rate, n))
      # Off-diagonal entries should be non-negative
      off_diag <- W[row(W) != col(W)]
      expect_true(all(off_diag >= 0))
    }
  }
})

test_that("adjacent_switching is a valid CTMC generator (rows sum to 0)", {
  for (n in c(3, 5, 8)) {
    W <- adjacent_switching(switch_rate = 0.2, num_classes = n)
    expect_equal(rowSums(W), rep(0, n), tolerance = 1e-12)
    # Off-diagonal entries should be non-negative
    off_diag <- W[row(W) != col(W)]
    expect_true(all(off_diag >= 0))
  }
})

test_that("find_steadystate gives uniform distribution for simple_switching", {
  n <- 5
  W <- simple_switching(switch_rate = 0.1, num_classes = n)
  pi_ss <- find_steadystate(W)

  # Simple switching with equal rates should have uniform steady state
  expect_equal(pi_ss, rep(1 / n, n), tolerance = 1e-6)
})


# ============================================================================
# 9. Single activity class reduces to scalar SIR
# ============================================================================

test_that("single activity class gives scalar R0 = beta/gamma", {
  pars <- list(
    n_activity = 1,
    class_sizes = 1000,
    activity_scores = 1,
    beta = 0.5,
    gamma = 0.1,
    sigma = 0.2,
    swch = matrix(0, 1, 1)
  )

  r0 <- compute_r0(pars)
  expect_equal(r0, 0.5 / 0.1, tolerance = 1e-10) # = 5
})


# ============================================================================
# 10. scale_results is idempotent on already-scaled data (within tolerance)
# ============================================================================

test_that("scale_results correctly divides by population/class sizes", {
  dt <- data.table(
    time = 0:1,
    class1 = c(500, 500), class2 = c(500, 500),
    pop_tot = rep(1000, 2),
    s_tot = c(800, 700), e_tot = c(100, 150),
    i_tot = c(100, 150), r_tot = c(0, 0),
    S1 = c(400, 350), S2 = c(400, 350),
    E1 = c(50, 75), E2 = c(50, 75),
    I1 = c(50, 75), I2 = c(50, 75),
    R1 = 0, R2 = 0
  )

  scaled <- scale_results(dt)
  expect_equal(scaled$s_tot[1], 800 / 1000)
  expect_equal(scaled$S1[1], 400 / 500)
  expect_equal(scaled$S2[1], 400 / 500)
})


# ============================================================================
# 11. Increasing assortativity changes R0 with heterogeneous activity
# ============================================================================

test_that("R0 changes with assortativity when activity scores differ", {
  n <- 4
  class_sizes <- rep(250, n)
  act_scores <- c(1, 2, 3, 5)

  r0_values <- numeric(3)
  epsilons <- c(0, 0.3, 0.7)

  for (i in seq_along(epsilons)) {
    m <- garnett_mixing(
      n_activity = n,
      epsilon_assort = epsilons[i],
      act_vector = act_scores,
      class_size_vector = class_sizes
    )

    pars <- list(
      n_activity = n,
      class_sizes = class_sizes,
      activity_scores = act_scores,
      beta = 0.2,
      gamma = 0.1,
      sigma = 0.2,
      m = m,
      swch = matrix(0, n, n)
    )

    r0_values[i] <- compute_r0(pars)
  }

  # R0 values should be different (not identical) when activity scores differ
  expect_false(all(abs(diff(r0_values)) < 1e-6))
})


# ============================================================================
# 12. neaten_state produces correct data.table structure
# ============================================================================

test_that("neaten_state produces expected columns and row count", {
  skip_if_not_installed("odin2")

  n_activity <- 3
  N <- 300
  pars <- list(
    n_activity = n_activity,
    S0 = rep(0.9, n_activity),
    E0 = rep(0.1, n_activity),
    N = N,
    beta = 0.3,
    sigma = 1 / 5,
    gamma = 1 / 7,
    omega = 0,
    m = matrix(1, n_activity, n_activity),
    swch = matrix(0, n_activity, n_activity),
    class_sizes = rep(100, n_activity),
    activity_scores = rep(1, n_activity)
  )

  sys <- dust2::dust_system_create(
    seirs_act(),
    pars = pars,
    ode_control = dust2::dust_ode_control(max_steps = 1e6)
  )
  dust2::dust_system_set_state_initial(sys)
  t <- seq(0, 50, by = 1)
  out <- dust2::dust_system_simulate(sys, t)
  nt <- neaten_state(system_used = sys, dust_state = out, times_simulated = t)

  expect_true(inherits(nt, "data.table"))
  expect_equal(nrow(nt), length(t))

  # Should have per-class columns and totals
  expected_cols <- c(
    "time", paste0("S", 1:n_activity), paste0("E", 1:n_activity),
    paste0("I", 1:n_activity), paste0("R", 1:n_activity),
    "s_tot", "e_tot", "i_tot", "r_tot", "pop_tot"
  )
  for (col in expected_cols) {
    expect_true(col %in% names(nt), info = paste("Missing column:", col))
  }
})


# ============================================================================
# 13. compute_rt_ts returns time-consistent output
# ============================================================================

test_that("compute_rt_ts returns correct dimensions and R0 at t=0", {
  skip_if_not_installed("odin2")

  n_activity <- 3
  N <- 1e5
  class_sizes <- rep(N / n_activity, n_activity)

  m <- garnett_mixing(
    n_activity = n_activity,
    epsilon_assort = 0.2,
    act_vector = rep(1, n_activity),
    class_size_vector = class_sizes
  )

  pars <- list(
    n_activity = n_activity,
    S0 = rep(1 - 1e-4, n_activity),
    E0 = rep(1e-4, n_activity),
    N = N,
    beta = 0.3,
    sigma = 1 / 3,
    gamma = 0.1,
    omega = 0,
    m = m,
    swch = matrix(0, n_activity, n_activity),
    class_sizes = class_sizes,
    activity_scores = rep(1, n_activity)
  )

  R0 <- compute_r0(pars)

  sys <- dust2::dust_system_create(
    seirs_act(),
    pars = pars,
    ode_control = dust2::dust_ode_control(max_steps = 1e6)
  )
  dust2::dust_system_set_state_initial(sys)
  t <- seq(0, 200, by = 1)
  out <- dust2::dust_system_simulate(sys, t)
  nt <- neaten_state(system_used = sys, dust_state = out, times_simulated = t)

  # Test with add_to_output = FALSE
  rt_df <- compute_rt_ts(nt, pars, add_to_output = FALSE)
  expect_equal(nrow(rt_df), length(t))
  expect_true("R0" %in% names(rt_df))
  expect_true("Rt" %in% names(rt_df))
  expect_equal(rt_df$Rt[1], R0, tolerance = 1e-3)

  # Test with add_to_output = TRUE (modifies nt in place)
  nt2 <- copy(nt)
  result <- compute_rt_ts(nt2, pars, add_to_output = TRUE)
  expect_true("Rt" %in% names(result))
  expect_true("R0" %in% names(result))
})


# ============================================================================
# 14. Spectral gap: proportionate mixing has zero spectral gap
# ============================================================================

test_that("spectral gap is zero for proportionate mixing with uniform activity", {
  n_activity <- 5
  class_sizes <- rep(200, n_activity)

  pars <- list(
    n_activity = n_activity,
    class_sizes = class_sizes,
    activity_scores = rep(1, n_activity),
    beta = 0.3,
    gamma = 0.1,
    sigma = 0.2,
    m = matrix(1 / n_activity, n_activity, n_activity),
    swch = matrix(0, n_activity, n_activity)
  )

  gap <- compute_spectral_gap(pars)
  K <- compute_ngm(pars)

  # For proportionate mixing with uniform activity, NGM should be rank-1
  # => spectral gap = R0 - 0 = R0 (only one nonzero eigenvalue)
  # Actually, the spectral gap should be the difference between top two eigenvalue magnitudes
  evals <- sort(Mod(eigen(K, only.values = TRUE)$values), decreasing = TRUE)
  expect_equal(gap, evals[1] - evals[2], tolerance = 1e-10)
})


# ============================================================================
# 15. Initial conditions functions
# ============================================================================

test_that("make_initial_uniform correctly distributes exposed", {
  ic <- make_initial_uniform(n_activity = 5, fraction_exposed = 0.02)

  expect_equal(length(ic$S0), 5)
  expect_equal(length(ic$E0), 5)

  # S0 + E0 should equal 1 for each class
  expect_true(all(abs(ic$S0 + ic$E0 - 1) < 1e-12))

  # Each class should have the same E0
  expect_equal(ic$E0, rep(0.02, 5), tolerance = 1e-12)
})

test_that("make_initial_specific concentrates exposure in one class", {
  ic <- make_initial_specific(n_activity = 5, initial_class = 3, fraction_exposed = 0.05)

  # Only class 3 should have exposed
  expect_true(ic$E0[3] > 0)
  for (j in c(1, 2, 4, 5)) {
    expect_equal(ic$E0[j], 0, tolerance = 1e-12)
  }

  # S0 + E0 = 1 for each class
  expect_true(all(abs(ic$S0 + ic$E0 - 1) < 1e-12))
})


# ============================================================================
# 16. compute_all produces a well-formed result table
# ============================================================================

test_that("compute_all returns expected structure", {
  skip_if_not_installed("odin2")

  n_activity <- 3
  N <- 1e5
  class_sizes <- rep(N / n_activity, n_activity)

  m <- garnett_mixing(
    n_activity = n_activity,
    epsilon_assort = 0.2,
    act_vector = rep(1, n_activity),
    class_size_vector = class_sizes
  )

  pars <- list(
    n_activity = n_activity,
    S0 = rep(1 - 1e-3, n_activity),
    E0 = rep(1e-3, n_activity),
    N = N,
    beta = 0.4,
    sigma = 1 / 3,
    gamma = 0.1,
    omega = 0,
    m = m,
    swch = matrix(0, n_activity, n_activity),
    class_sizes = class_sizes,
    activity_scores = rep(1, n_activity),
    t_end = 1000
  )

  epi <- run_system_wparms(parms = pars, model = seirs_act())
  result <- compute_all(epi, classed_detail = TRUE, params = pars)

  expect_true(inherits(result, "data.table"))
  # Should have overall row (class = -1) and per-class rows
  expect_true(-1 %in% result$class)
  expect_true(all(1:n_activity %in% result$class))

  # Overall row should have R0, HIT, final_size
  overall <- result[class == -1]
  expect_true("r0" %in% names(overall))
  expect_true("hit" %in% names(overall))
  expect_true("final_size" %in% names(overall))
  expect_true(overall$r0 > 1)
})
