test_that("heterogeneous no-switch final size <= homogeneous final size for same R0", {
  num_activity <- 5
  has_scores <- c(0.936825, 2.182864, 3.492136, 5.271117, 9.467058)
  mean_act <- mean(has_scores)
  scalars <- basic_scalar_parms("sars-cov-2")

  # homogeneous params
  hom <- list()
  hom$n_activity <- num_activity
  hom$N <- 1e5
  hom$t_end <- 800
  hom$class_sizes <- rep(hom$N / num_activity, num_activity)
  hom$activity_scores <- rep(mean_act, num_activity)
  hom$epsilon <- 0
  hom$swch <- simple_switching(switch_rate = 0, num_classes = num_activity)
  hom$m <- garnett_mixing(
    n_activity = num_activity, epsilon_assort = 0,
    act_vector = hom$activity_scores,
    class_size_vector = hom$class_sizes
  )
  ic <- make_initial_conditions(
    n_activity = num_activity,
    initial_condition_rule = "uniform",
    fraction_exposed = 1e-3
  )
  hom$E0 <- ic$E0
  hom$S0 <- ic$S0
  hom$omega <- 0
  hom$beta <- scalars$beta
  hom$sigma <- scalars$sigma
  hom$gamma <- scalars$gamma

  # heterogeneous, no switching (activity heterogeneity present)
  het <- hom
  het$activity_scores <- has_scores
  het$m <- garnett_mixing(
    n_activity = num_activity, epsilon_assort = 0,
    act_vector = het$activity_scores,
    class_size_vector = het$class_sizes
  )

  target_r0 <- 2
  hom_cal <- calibrate_parms(hom, target_r0)
  het_cal <- calibrate_parms(het, target_r0)

  epi_hom <- run_system_wparms(parms = hom_cal, model = seirs_act())
  res_hom <- compute_all(epi_hom, classed_detail = FALSE, params = hom_cal)
  overall_hom <- as.data.frame(res_hom)
  overall_hom <- overall_hom[overall_hom$class == -1, , drop = FALSE]

  epi_het <- run_system_wparms(parms = het_cal, model = seirs_act())
  res_het <- compute_all(epi_het, classed_detail = FALSE, params = het_cal)
  overall_het <- as.data.frame(res_het)
  overall_het <- overall_het[overall_het$class == -1, , drop = FALSE]

  # sanity
  expect_equal(overall_hom$r0, target_r0, tolerance = 1e-6)
  expect_equal(overall_het$r0, target_r0, tolerance = 1e-6)

  # heterogeneous final size should not exceed homogeneous final size for same R0
  expect_lte(overall_het$final_size, overall_hom$final_size + 1e-6)
  # both final sizes are between 0 and 1
  expect_true(overall_hom$final_size > 0 && overall_hom$final_size < 1)
  expect_true(overall_het$final_size > 0 && overall_het$final_size < 1)
})
