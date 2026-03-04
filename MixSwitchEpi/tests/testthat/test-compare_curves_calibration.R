test_that("compare_curves parameter sets share expected calibration properties", {
  skip_if_not_installed("odin2")

  N <- 1e6
  target_r0 <- 2
  num_activity <- 5
  haslemere_activity_scores <- c(
    0.936825,
    2.182864,
    3.492136,
    5.271117,
    9.467058
  )

  # Resolve data file relative to repo root, robust under both devtools::test()
  # and Rscript invocation from the repo root
  pkg_root <- tryCatch(
    rprojroot::find_package_root_file(),
    error = function(e) getwd()
  )
  path_to_empirical <- file.path(pkg_root, "data", "Haslemere.csv")
  skip_if_not(
    file.exists(path_to_empirical),
    message = paste("Required empirical contact file not found at", path_to_empirical)
  )

  parsets <- make_parameter_sets(
    N = N,
    r0 = target_r0,
    num_activity = num_activity,
    haslemere_activity_scores = haslemere_activity_scores,
    path_to_empirical_mix = path_to_empirical
  )

  model_names <- vapply(parsets, function(p) p$name, character(1))

  expect_setequal(
    model_names,
    c(
      "Homogeneous population",
      "Heterogeneous pop (proportionate mix)",
      "Heterogeneous pop (empirical mix)",
      "With switching"
    )
  )

  # Extract calibrated parameter lists
  cal_pars <- lapply(parsets, function(p) p$calibrated_parms)
  names(cal_pars) <- model_names

  # Helper to compute R0 safely
  get_r0 <- function(p) {
    compute_r0(p)
  }

  r0_vals <- vapply(cal_pars, get_r0, numeric(1))

  # 1) Homogeneous case is calibrated to target R0
  expect_equal(
    r0_vals[["Homogeneous population"]],
    target_r0,
    tolerance = 1e-6,
    label = "Homogeneous R0 should equal target_r0"
  )

  # # 2) Heterogeneous proportionate mix is calibrated to the SAME target R0
  # expect_equal(
  #   r0_vals[["Heterogeneous pop (proportionate mix)"]],
  #   target_r0,
  #   tolerance = 1e-6,
  #   label = "Proportionate-mix R0 should equal target_r0"
  # )

  # 3) Empirical mix and switching cases inherit homogeneous beta (fixed-beta family)
  beta_vals <- vapply(cal_pars, function(p) p$beta, numeric(1))
  beta_hom <- beta_vals[["Homogeneous population"]]

  expect_equal(
    beta_vals[["Heterogeneous pop (empirical mix)"]],
    beta_hom,
    tolerance = 1e-12,
    label = "Empirical-mix beta should equal homogeneous beta"
  )
  expect_equal(
    beta_vals[["With switching"]],
    beta_hom,
    tolerance = 1e-12,
    label = "With-switching beta should equal homogeneous beta"
  )

  # 4) Under the shared beta, empirical + switching cases should generally
  #    have R0 values different from the homogeneous case (they encode
  #    structural differences in mixing/switching).
  expect_false(
    isTRUE(all.equal(
      r0_vals[["Heterogeneous pop (empirical mix)"]],
      r0_vals[["Homogeneous population"]]
    )),
    label = "Empirical-mix R0 should differ from homogeneous R0 under shared beta"
  )
  expect_false(
    isTRUE(all.equal(
      r0_vals[["With switching"]],
      r0_vals[["Homogeneous population"]]
    )),
    label = "With-switching R0 should differ from homogeneous R0 under shared beta"
  )

  # 5) Smoke test: compare_curves_main() uses these calibrated parameters
  #    without error and produces a metrics_df consistent with the above.
  out <- compare_curves_main(
    N = N,
    r0 = target_r0,
    num_activity = num_activity,
    haslemere_activity_scores = haslemere_activity_scores,
    path_to_empirical_mix = path_to_empirical,
    output_dir = NULL,
    save = FALSE
  )

  md <- out$metrics_df
  expect_true(all(c("model_type", "r0") %in% names(md)))

  # Ensure R0 column matches direct compute_r0 values for each calibrated set
  for (i in seq_len(nrow(md))) {
    mt <- as.character(md$model_type[i])
    r0_direct <- r0_vals[[mt]]
    expect_equal(
      md$r0[i], r0_direct,
      tolerance = 1e-6,
      label = paste("metrics_df R0 mismatch for", mt)
    )
  }
})
