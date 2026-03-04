#' Sanitise the exposed fraction
#'
#' Makes sure that the exposed fraction is between
#' 0 and 1.
#' @param fraction_exposed the fraction of the pop
#' which is put in the E compartment(s)
#' @export
check_fraction_exposed <- function(fraction_exposed) {
  stopifnot(
    assertthat::assert_that(0 <= fraction_exposed & 1 >= fraction_exposed)
  )
}

#' Initial conditions where each compartment is inoculated simultaenously with
#' an equal fraction of exposed individuals (or, rather, stuff...)
#'
#' This is a simple initial condition,
#' Everything else is in the Susceptible
#' compartment.
#'
#' @param n_activity the number of activity classes in the model
#' @param fraction_exposed the fraction of each activity class which is
#' exposed to the infection.
#' @returns list, whose items are `E0` and `S0`.
#' These arrays are of length (n_activity).
#' @export
make_initial_uniform <- function(n_activity, fraction_exposed) {
  # Return values as fractions of each class. odin initialisers use
  # `initial(S[]) <- S0[i] * class_sizes[i]` so S0 and E0 should be the
  # fraction of that class in S and E respectively. If the requested
  # `fraction_exposed` is a fraction of the total population and exposure
  # is distributed uniformly across classes, then the fraction exposed in
  # each class equals `fraction_exposed`.
  e0 <- rep(fraction_exposed, n_activity)
  s0 <- rep(1 - fraction_exposed, n_activity)
  list(
    E0 = e0,
    S0 = s0
  )
}

#' Specific activity class infected first.
#'
#' initial condition where we only expose one activity class (user-specified)
#' rather than exposing all activity classes initially.
#' This has nontrivial consequences
#' for the ensuing dynamics,
#'
#' @param n_activity The number of activity classes in the model.
#' @param initial_class the initial class which contains the exposeds
#' @param fraction_exposed The fraction of the activity class which
#'  is moved to the Exposed
#' (E) compartment.
#'
#' @export
make_initial_specific <- function(n_activity,
                                  initial_class,
                                  fraction_exposed) {
  check_fraction_exposed(fraction_exposed)
  # For a specific-class exposure, set E0 so that E0 * N = fraction_exposed * class_size
  # E0 should be the fraction of the specified class that is exposed.
  e0 <- rep(0, n_activity)
  e0[initial_class] <- fraction_exposed
  # Susceptible fraction per class: other classes are fully susceptible (1),
  # the infected class has 1 - fraction_exposed susceptible.
  s0 <- rep(1, n_activity)
  s0[initial_class] <- 1 - fraction_exposed
  list(E0 = e0, S0 = s0)
}


#' Report the parameter scan results
#'
#' Writes out the summary statistics for the parameter scan,
#' along with the time at which
#' started and completed. Makes a directory for this.
#'
#' perhaps should be updated to measure the amount of cpu time the scan takes?
#'
#' @param results_tibble a tibble/df containing the results of a parameter scan
#' @param parameter_tibble the tibble containing the parameters used for that
#' scan.
#' @returns NULL
#'
#' @export
report_scan_results <- function(results_tibble,
                                parameter_tibble,
                                create_dir = FALSE) {
  if (create_dir == TRUE) {
    pathtime <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    dir.create(file.path(pathtime))
    setwd(file.path(pathtime))

    write.csv(results_tibble, file = "scan_results.csv", sep = ",")
    write.csv(parameter_tibble, file = "parameters.csv", sep = ",")
    message(glue::glue("wrote out results data to directory {pathtime}"))
  }
}


#' Scan 2D parameter space (epsilon × switch_rate.
#'
#' Basic flexible framework function which can be adapted significantly
#' as we require. Runs a set of workers in parallel.
#' We precompute the mixing and switching matrices to reduce
#' computational time.
#'
#' @param base_pars Calibrated base parameters from default_setup().
#' @param epsilon_vals Vector of epsilon values to scan (default: 0 to 1)
#' @param switch_vals Vector of switching rate values to scan (default: 0 to 1)
#' @param n_epsilon Number of epsilon points (if epsilon_vals not provided)
#' @param n_switch Number of switch rate points (if switch_vals not provided)
#' @param classed_detail Whether to compute class-specific metrics
#' @param switch_model the model of switching to be used.
#' allowed values include:
#' @param mix_model the model for
#' @param initial_condition_rule the rule  to generate the initial cond.
#' @param frac_infected_initial the initial fraction infected of the population
#' This will be put in the E (exposed) compartment.
#' @param calibration_mode mode of calibration. allowed values are 'none', 'switch', 'assort', 'both'.
#' If 'none', no recalibration. If 'switch', calibrate once for each switch value to the target r0.
#' If 'assort', dont calibrate for switch but calibrate for assort. If 'both', calibrate each point in mix-switch space.
#' @param t_end simulation horizon in days (default 1200)
#' @returns data.table with all epidemic metrics for each parameter combination
#' @export
scan_mixswitch_framework <- function(base_pars,
                                     epsilon_vals,
                                     switch_vals,
                                     n_epsilon = 11,
                                     n_switch = 11,
                                     classed_detail = FALSE,
                                     switch_model = "simple_switching",
                                     mix_model, power_param = NULL,
                                     initial_condition_rule = NULL,
                                     frac_infected_initial = 1e-4,
                                     calibration_mode = "none",
                                     t_end = 1200) {
  base_pars$t_end <- t_end
  n_activity <- base_pars$n_activity
  cli::cli_alert_info(glue::glue("Precomputing contact matrices using {mix_model}"))
  cli::cli_alert_info(glue::glue("calibration mode: {calibration_mode}"))
  cli::cli_alert_info(glue::glue("Precomputing matrices {switch_model})"))
  if (switch_model == "simple_switching") {
    switching_mats <- lapply(switch_vals, function(sw) {
      MixSwitchEpi::simple_switching(switch_rate = sw, num_classes = n_activity)
    })
  } else if (switch_model == "adjacent_switching_wrapped") {
    switching_mats <- lapply(switch_vals, function(sw) {
      MixSwitchEpi::adjacent_switching_wrapped(
        switch_rate = sw,
        num_classes = n_activity
      )
    })
  } else {
    cli::cli_alert_danger(glue::glue("unrecognised model {switch_model}.
    allowed: `simple_switching`, `adjacent_switching_wrapped`"))
    stop("unrecognised switching model")
  }
  if (is.null(epsilon_vals)) {
    epsilon_vals <- seq(0, 1, length.out = n_epsilon)
  }
  cli::cli_alert_info(glue::glue("Precomputing contact matrices
   using {mix_model}"))
  if (mix_model == "garnett_mixing") {
    mixing_mats <- lapply(epsilon_vals, function(eps) {
      MixSwitchEpi::garnett_mixing(
        n_activity = n_activity,
        epsilon_assort = eps,
        act_vector = base_pars$activity_scores,
        class_size_vector = base_pars$class_sizes
      )
    })
  } else if (mix_model == "exponential") {
    mixing_mats <- lapply(epsilon_vals, function(eps) {
      MixSwitchEpi::exponential_mixing(
        n_activity = n_activity,
        epsilon_assort = eps,
        act_vector = base_pars$activity_scores,
        class_size_vector = base_pars$class_sizes
      )
    })
  } else if (mix_model == "polynomial") {
    mixing_mats <- lapply(epsilon_vals, function(eps) {
      MixSwitchEpi::polynomial_mixing(
        n_activity = n_activity,
        power_param = power_param,
        epsilon_assort = eps,
        act_vector = base_pars$activity_scores,
        class_size_vector = base_pars$class_sizes
      )
    })
  } else {
    stop(glue::glue("unrecognised mixing model {mix_model}"))
  }
  names(mixing_mats) <- as.character(epsilon_vals)
  names(switching_mats) <- as.character(switch_vals)

  n_activity <- base_pars$n_activity

  param_grid <- expand.grid(
    epsilon = epsilon_vals,
    switch_rate = switch_vals,
    stringsAsFactors = FALSE
  )

  run_single <- function(i) {
    eps <- param_grid$epsilon[i]
    sw <- param_grid$switch_rate[i]
    pars <- base_pars
    pars$m <- mixing_mats[[as.character(eps)]]
    pars$swch <- switching_mats[[as.character(sw)]]
    pars$epsilon <- eps
    pars$switch_rate <- sw
    if (calibration_mode != "none") {
      calib_pars <- base_pars
      calib_eps <- if (calibration_mode %in% c("assort", "both")) eps else base_pars$epsilon
      calib_sw <- if (calibration_mode %in% c("switch", "both")) sw else base_pars$switch_rate

      if (mix_model == "garnett_mixing") {
        calib_pars$m <- MixSwitchEpi::garnett_mixing(n_activity, calib_eps, base_pars$activity_scores, base_pars$class_sizes)
      } else if (mix_model == "exponential") {
        calib_pars$m <- MixSwitchEpi::exponential_mixing(n_activity, calib_eps, base_pars$activity_scores, base_pars$class_sizes)
      } else if (mix_model == "polynomial") {
        calib_pars$m <- MixSwitchEpi::polynomial_mixing(n_activity, power_param, calib_eps, base_pars$activity_scores, base_pars$class_sizes)
      }

      if (switch_model == "simple_switching") {
        calib_pars$swch <- MixSwitchEpi::simple_switching(switch_rate = calib_sw, num_classes = n_activity)
      } else if (switch_model == "adjacent_switching_wrapped") {
        calib_pars$swch <- MixSwitchEpi::adjacent_switching_wrapped(switch_rate = calib_sw, num_classes = n_activity)
      }

      calib_pars$epsilon <- calib_eps
      calib_pars$switch_rate <- calib_sw

      calibrated <- MixSwitchEpi::calibrate_parms(
        preliminary_parms = calib_pars,
        target_r0 = base_pars$target_r0
      )
      pars$beta <- calibrated$beta
    }
    model_output <- MixSwitchEpi::run_system_wparms(
      parms = pars,
      model = seirs_act()
    )
    metrics <- MixSwitchEpi::compute_all(model_output, params = pars)
    metrics$epsilon <- eps
    metrics$switch_rate <- sw
    metrics
  }
  message(glue::glue("Running {nrow(param_grid)} parameter combinations..."))
  results_list <- parallel::mclapply(
    seq_len(nrow(param_grid)),
    run_single,
    mc.cores = 10,
    mc.preschedule = TRUE
  )
  cli::cli_alert_success("Scan complete!")
  results_list
}
