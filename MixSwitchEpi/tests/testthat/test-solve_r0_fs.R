test_that("solve_r0_fs agrees with finalsize::final_size and homogeneous sim", {
  # numeric agreement with finalsize package for a few R0 values
  r0_values <- c(1.5, 2, 3)
  for (r0 in r0_values) {
    # allow a small numerical tolerance between solvers
    expect_equal(
      solve_r0_fs(r0),
      finalsize::final_size(r0 = r0)$p_infected,
      tolerance = 1e-4
    )
  }

  # quick homogeneous simulation check (small, deterministic)
  num_activity <- 5
  has_scores <- c(0.936825, 2.182864, 3.492136, 5.271117, 9.467058)
  mean_act <- mean(has_scores)
  scalars <- basic_scalar_parms("sars-cov-2")

  prelim <- list()
  prelim$beta <- scalars$beta
  prelim$sigma <- scalars$sigma
  prelim$gamma <- scalars$gamma
  prelim$n_activity <- num_activity
  prelim$N <- 1e6
  prelim$t_end <- 1200
  prelim$class_sizes <- rep(prelim$N / num_activity, num_activity)
  prelim$activity_scores <- rep(mean_act, num_activity)
  prelim$epsilon <- 0
  prelim$swch <- simple_switching(switch_rate = 0, num_classes = num_activity)
  prelim$m <- garnett_mixing(
    n_activity = num_activity, epsilon_assort = 0,
    act_vector = prelim$activity_scores,
    class_size_vector = prelim$class_sizes
  )
  ic <- make_initial_conditions(
    n_activity = num_activity,
    initial_condition_rule = "uniform",
    fraction_exposed = 1e-3
  )
  prelim$E0 <- ic$E0
  prelim$S0 <- ic$S0
  prelim$omega <- 0

  target_r0 <- 2
  cal <- calibrate_parms(prelim, target_r0 = target_r0)

  epi_out <- run_system_wparms(parms = cal, model = seirs_act())
  res <- compute_all(epi_out, classed_detail = FALSE, params = cal)
  overall <- as.data.frame(res)
  overall <- overall[overall$class == -1, , drop = FALSE]

  expect_equal(overall$r0, target_r0, tolerance = 1e-6)
  expect_equal(overall$final_size, solve_r0_fs(overall$r0), tolerance = 5e-4)
})
