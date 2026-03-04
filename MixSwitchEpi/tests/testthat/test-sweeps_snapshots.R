library(testthat)
library(MixSwitchEpi)

test_that("scan_mixswitch_framework produces stable sweep results", {
  skip_if_not_installed("odin2")

  base_pars <- make_parms(
    n_activity = 3,
    epsilon = 0.0,
    switch_rate = 0.0,
    initial_condition_rule = "uniform",
    infected_fraction = 1e-03,
    activity_scheme = "null",
    mixing_model = "garnett_mixing",
    switching_model = "simple_switching",
    pathogen_string = "sars-cov-2",
    target_r0 = 2
  )

  epsilon_vals <- c(0, 0.5, 1)
  switch_vals <- c(0, 0.1)
  
  out <- scan_mixswitch_framework(
    base_pars = base_pars,
    epsilon_vals = epsilon_vals,
    switch_vals = switch_vals,
    classed_detail = FALSE,
    switch_model = "simple_switching",
    mix_model = "garnett_mixing"
  )

  results_dt <- data.table::rbindlist(out, fill = TRUE)
  
  # Snapshot the first few rows but avoiding time-dependent fields if they differ 
  expect_snapshot(head(results_dt))
})

test_that("haslemere_empirical produces stable sweep results", {
  skip_if_not_installed("odin2")

  pkg_root <- tryCatch(
    rprojroot::find_package_root_file(),
    error = function(e) getwd()
  )
  path_to_empirical <- file.path(pkg_root, "data", "Haslemere.csv")
  skip_if_not(
    file.exists(path_to_empirical),
    message = paste("Required empirical contact file not found at", path_to_empirical)
  )

  # Small R0 list for faster test execution
  out <- haslemere_empirical(
    target_r0 = 2,
    calibrate_each = FALSE,
    contact_file = path_to_empirical
  )
  # Look at the metrics shape
  res <- out$results
  # Snapshot the first few structured results
  expect_snapshot(head(res))
})
