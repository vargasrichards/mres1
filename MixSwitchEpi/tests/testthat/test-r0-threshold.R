test_that("Epidemic threshold: disease dies out when R0 < 1", {
  # Scenario: R0 = 0.8
  target_r0 <- 0.8
  n_activity <- 3

  params <- make_parms(
    n_activity = n_activity,
    epsilon = 0.5,
    switch_rate = 0.1,
    initial_condition_rule = "uniform",
    infected_fraction = 0.001,
    activity_scheme = "linear",
    mixing_model = "garnett_mixing",
    switching_model = "simple_switching",
    target_r0 = target_r0
  )

  # Ensure the calculated R0 is indeed what we requested
  expect_equal(compute_r0(params), target_r0, tolerance = 1e-6)

  # Run simulation
  output <- run_system_wparms(params, seirs_act())
  metrics <- characterise_end(output, classed_detail = FALSE)

  # Final size should be very close to the initial infected fraction (no significant spread)
  # characterise_end returns final_size as a fraction of total population
  expect_lt(metrics$final_size, 0.005) # Very small spread
})

test_that("Epidemic threshold: epidemic occurs when R0 > 1", {
  # Scenario: R0 = 1.5
  target_r0 <- 1.5
  n_activity <- 3

  params <- make_parms(
    n_activity = n_activity,
    epsilon = 0.5,
    switch_rate = 0.1,
    initial_condition_rule = "uniform",
    infected_fraction = 0.001,
    activity_scheme = "linear",
    mixing_model = "garnett_mixing",
    switching_model = "simple_switching",
    target_r0 = target_r0
  )

  # Ensure the calculated R0 is indeed what we requested
  expect_equal(compute_r0(params), target_r0, tolerance = 1e-6)

  # Run simulation
  output <- run_system_wparms(params, seirs_act())
  metrics <- characterise_end(output, classed_detail = FALSE)

  # Final size should be significantly larger than the initial seed
  expect_gt(metrics$final_size, 0.1) # Significant spread

  # Check for a peak
  peak_metrics <- characterise_peaks(output, get_info(output), classed_detail = FALSE, pars = params)
  expect_gt(peak_metrics$psI, params$E0[1] + 1e-6)
})

test_that("Threshold behavior is consistent across mixing models", {
  models <- c("garnett_mixing", "exponential", "polynomial")

  for (m in models) {
    # R0 < 1
    p_low <- make_parms(
      n_activity = 3, epsilon = 0.2, switch_rate = 0.05,
      initial_condition_rule = "uniform", infected_fraction = 0.001,
      mixing_model = m, switching_model = "simple_switching", target_r0 = 0.7,
      activity_scheme = "linear",
      power_param = if (m == "polynomial") 2 else NULL
    )
    out_low <- run_system_wparms(p_low, seirs_act())
    fs_low <- characterise_end(out_low, classed_detail = FALSE)$final_size
    expect_lt(fs_low, 0.005, label = paste("fs_low for", m))

    # R0 > 1
    p_high <- make_parms(
      n_activity = 3, epsilon = 0.2, switch_rate = 0.05,
      initial_condition_rule = "uniform", infected_fraction = 0.001,
      mixing_model = m, switching_model = "simple_switching", target_r0 = 2.0,
      activity_scheme = "linear",
      power_param = if (m == "polynomial") 2 else NULL
    )
    out_high <- run_system_wparms(p_high, seirs_act())
    fs_high <- characterise_end(out_high, classed_detail = FALSE)$final_size
    expect_gt(fs_high, 0.5, label = paste("fs_high for", m))
  }
})
