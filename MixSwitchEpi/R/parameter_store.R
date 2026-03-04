#' Get the scalar epi parameters for given pathogen
#'
#' Convenience function which returns the scalar parameters
#' (ie transmission rate (beta),
#' rate of progression to I (sigma) and the recovery rate (gamma).
#'
#' @param pathogen_string string describing the desired pathogen
#' for which the parameters
#' are required. Default value is "sars-cov-2" for which we use the parameterisation
#' from (Britton et al., 2020) Science
#' @returns a list of scalars whose elements are `$beta`, `$sigma`, `$gamma`
#' @family parameterisation
#' @export
basic_scalar_parms <- function(pathogen_string = "sars-cov-2") {
  if (pathogen_string == "sars-cov-2") {
    # cli::cli_alert_info("using Britton et al., 2020 Science
    # values for scalar parameters.")
    scalars <- list(
      beta = 0.1,
      sigma = 1 / 3,
      gamma = 1 / 4
    )
  } else if (pathogen_string == "longer_recovery") {
    scalars <- list(
      beta = 0.1,
      sigma = 1 / 3,
      gamma = 1 / 8
    )
  } else {
    cli::cli_alert_danger(glue::glue("unknown pathogen {pathogen_string}
    specified. please specify e.g., `sars-cov-2` or `longer_recovery`."))
    stop()
  }
  scalars
}

#' Changes the β value until the target R0 is met
#'
#' Uses a very simple approach by dividing the target R0 by the current R0
#' and scaling β by that factor.
#'
#' @param preliminary_parms list of preliminary parameters
#' including mixing matrix,
#' activity scores, class sizes, etc.
#' @param target_r0 desired R0 value
#' @returns updated parameters with β set to achieve target R0
#' @family calibration
#' @export
calibrate_parms <- function(preliminary_parms, target_r0) {
  current_r0 <- compute_r0(preliminary_parms)
  scaling_factor <- target_r0 / current_r0
  preliminary_parms$beta <- preliminary_parms$beta * scaling_factor
  assertthat::assert_that(
    abs(compute_r0(preliminary_parms) - target_r0) < 1e-6,
    msg = "Failed to calibrate parameters to desired R0"
  )
  calibrated_parms <- preliminary_parms
  calibrated_parms
}

#' Make the switching matrix according to supplied model
#'
#' Wrapper - mostly for concision - which sanitises the
#' switch rate input into the function and uses the supplied
#' model to produce the CTMC generator matrix supplied to odin
#'
#' @param n_activity the number of activity classes in the odin model.
#' @param switch_rate the rate of switching. this = 1/Expected residence time
#' in a particular activity class
#' @param switching_model allowed values are `simple_switching` and `adjacent_
#' switching`.
#' @export
#' @family switching_model
#' @examples
#' # Example call
#' swch <- make_switch_matrix(
#'   n_activity = 5,
#'   switch_rate = 0.1,
#'   switching_model = "simple_switching"
#' )
#' print(swch)
#' # Example (adjacent wrap-around switching)
#' swch <- make_switch_matrix(
#'   n_activity = 5,
#'   switch_rate = 0.1,
#'   switching_model = "adjacent_switching"
#' )
#' print(swch)
make_switch_matrix <- function(n_activity,
                               switch_rate,
                               switching_model) {
  stopifnot(assertthat::assert_that(
    0 <= switch_rate & switch_rate <= 1
  ))
  # cli::cli_alert_info((glue::glue("using
  # switching model {switching_model} at rate {switch_rate}")))
  if (switching_model == "simple_switching") {
    sw <- MixSwitchEpi::simple_switching(
      switch_rate = switch_rate,
      num_classes = n_activity
    )
  } else if (switching_model == "adjacent_switching") {
    sw <- MixSwitchEpi::adjacent_switching(
      switch_rate = switch_rate,
      num_classes = n_activity
    )
  } else {
    cli::cli_alert_danger(glue::glue("unrecognised switching model
     {switching_model} supplied. allowed values are `simple_switching`
     and `adjacent_switching`"))
    stop()
  }
  sw
}

#' Construct the contact matrix
#'
#'
#' @param n_activity the number of activity groups in the model
#' @param epsilon the assortativity for the contact/mixing model
#' @param activity_vector the vector of activity scores
#' @param class_size_vector the vector of class sizes in the model
#' @param power_param the power parameter for
#' @param mixing_model the model of contact/mixing
#' @export
#' @family contact_models
make_contact_matrix <- function(n_activity,
                                epsilon,
                                activity_vector,
                                class_size_vector,
                                power_param = NULL, # needed for polynomial model
                                mixing_model) {
  # cli::cli_alert_info(glue::glue("making contact matrix for parametric model
  # {mixing_model}"))
  if (mixing_model == "garnett_mixing") {
    contact_mat <- MixSwitchEpi::garnett_mixing(
      n_activity = n_activity,
      epsilon_assort = epsilon,
      act_vector = activity_vector,
      class_size_vector = class_size_vector
    )
  } else if (mixing_model == "polynomial") {
    contact_mat <- polynomial_mixing(
      n_activity = n_activity,
      power_param = power_param,
      epsilon_assort = epsilon,
      act_vector = activity_vector,
      class_size_vector = class_size_vector
    )
  } else if (mixing_model == "exponential") {
    contact_mat <- exponential_mixing(
      n_activity = n_activity,
      epsilon_assort = epsilon,
      act_vector = activity_vector,
      class_size_vector = class_size_vector
    )
  } else {
    cli::cli_alert_danger(glue::glue("unrecognised mixmat model
    {mixing_model} supplied. allowed values are `garnett_mixing`,`exponential`
    and `polynomial`"))
    stop()
  }
  contact_mat
}

#' Make the initial conditions
#'
#' uses the user-supplied initial condition rule
#' and the fraction infected, to generate initial
#' conditions for forward simulation
#'
#' @param initial_condition_rule the rule used to
#' @param frac_infected fraction infected at first
#' @param infected_class if the initial condition rule
#'
#' @export
make_initial_conditions <- function(n_activity,
                                    initial_condition_rule,
                                    fraction_exposed,
                                    infected_class = NULL) {
  if (initial_condition_rule == "uniform") {
    ic <- make_initial_uniform(
      n_activity = n_activity,
      fraction_exposed = fraction_exposed
    )
  } else if (initial_condition_rule == "specific") {
    ic <- make_initial_specific(
      n_activity = n_activity,
      initial_class = infected_class,
      fraction_exposed = fraction_exposed
    )
  } else {
    cli::cli_alert_danger(glue::glue("Unerecognised initial condition rule
    {initial_condition_rule} supplied."))
    stop()
  }
  ic
}

#' Make a default parameterisaiton
#'
#' This function makes the parameterisation for mix-switch
#' simulation.
#' N (total population) is set to 1 million people by default.
#' Switching rate has default 0.1.
#' Initial conditions are set according to the specified rule
#' Function assumes that the class sizes are of equal size.
#'
#' @param n_activity number of activity classes
#' @param epsilon assortativity parameter for mixing matrix
#' @param switch_rate rate of switching between activity classes
#' @param initial_condition_rule rule for initial conditions. allowed values are
#' `uniformly_infected` and `specific_class_infected`
#' @param infected_class if `initial_condition_rule` is
#' `specific_class_infected`,
#' this parameter specifies which class is initially infected.
#' @param activity_scheme scheme for activity scores. Allowed values are
#' `linear`, `null` and `gamma`. Default is `gamma`.
#' This specifies how activity scores
#' are generated for each activity class.
#' @param gamma_dist if `activity_scheme` is `gamma`, this parameter
#' is a list containing the parameters for the gamma distribution.
#'  It can contain either
#' `mean` and `variance`, or `shape` and `rate`. (But not a mixture! :p)
#' @param power_param optional power parameter for the `polynomial` mixing model.
#'
#' @param pathogen_string string describing the pathogen to be modelled.
#' Currently only `sars-cov-2` is implemented.
#' @returns a list of parameters for the simulation.
#' @export
make_parms <- function(n_activity,
                       epsilon,
                       switch_rate,
                       initial_condition_rule,
                       infected_class = NULL,
                       infected_fraction,
                       activity_scheme = "gamma",
                       gamma_dist = NULL,
                       mixing_model,
                       switching_model,
                       target_r0,
                       power_param = NULL,
                       pathogen_string = "sars-cov-2") {
  pars <- list()
  scalars <- basic_scalar_parms(pathogen_string)
  pars$beta <- scalars$beta
  pars$sigma <- scalars$sigma
  pars$gamma <- scalars$gamma

  pars$n_activity <- n_activity
  pars$epsilon <- epsilon
  pars$switch_rate <- switch_rate
  pars$N <- 1e6
  pars$t_end <- 1000
  pars$class_sizes <- rep(pars$N / n_activity, n_activity)

  pars$swch <- MixSwitchEpi::make_switch_matrix(
    n_activity = n_activity,
    switch_rate = switch_rate,
    switching_model = switching_model
  )

  if (activity_scheme == "linear") {
    activity_vector <- activity_linear(
      min_act = 1,
      max_act = n_activity,
      n_classes = n_activity
    )
  } else if (activity_scheme == "null") {
    activity_vector <- activity_null(n_classes = n_activity)
  } else if (activity_scheme == "gamma") {
    if (is.null(gamma_dist)) {
      stop("please supply gamma_dist parameter when activity_scheme = 'gamma'")
    }
    shape <- NULL
    rate <- NULL
    mean <- NULL
    variance <- NULL
    if (!is.null(gamma_dist$shape)) {
      shape <- gamma_dist$shape
      rate <- gamma_dist$rate
    } else if (!is.null(gamma_dist$mean)) {
      mean <- gamma_dist$mean
      variance <- gamma_dist$variance
    } else {
      stop("gamma_dist must contain either (shape and rate)
      or (mean and variance)")
    }

    cutset <- qgamma(seq(0, 1, length.out = n_activity + 1),
      shape = if (!is.null(shape)) shape else (mean^2) / variance,
      rate = if (!is.null(rate)) rate else mean / variance
    )
    activity_vector <- find_gamma_expectations(
      cutoffs = cutset,
      mean = mean,
      variance = variance,
      shape = shape,
      rate = rate,
      num_bins = n_activity
    )
  } else {
    stop("unknown activity scheme.
    Please choose from 'linear', 'null' or 'gamma'")
  }
  pars$activity_scores <- activity_vector

  initial_conditions <- make_initial_conditions(n_activity,
    initial_condition_rule,
    infected_fraction,
    infected_class = infected_class
  )
  pars$E0 <- initial_conditions$E0
  pars$S0 <- initial_conditions$S0

  pars$m <- make_contact_matrix(
    n_activity = n_activity,
    epsilon = epsilon,
    activity_vector = activity_vector,
    class_size_vector = pars$class_sizes,
    mixing_model = mixing_model,
    power_param = power_param
  )
  pars$omega <- 0
  pars$target_r0 <- target_r0
  calibrated_pars <- calibrate_parms(
    preliminary_parms = pars,
    target_r0 = target_r0
  )
  calibrated_pars
}

#' Find expected values (means) within gamma distribution bins
#'
#' @param cutoffs Vector of bin boundaries
#' @param shape Shape parameter of gamma distribution (optional
#' if mean/variance provided)
#' @param rate Rate parameter (optional if mean/variance provided)
#' @param scale Scale parameter (alternative to rate)
#' @param mean Mean of gamma distribution (alternative to shape/rate)
#' @param variance Variance of gamma distribution (alternative to shape/rate)
#' @param num_bins Number of bins (not used, just for compatibility)
#' @return Vector of conditional expected values for each bin
#' @export
find_gamma_expectations <- function(cutoffs, shape = NULL, rate = NULL,
                                    scale = NULL, mean = NULL,
                                    variance = NULL, num_bins = NULL) {
  if (!is.null(mean) && !is.null(variance)) {
    shape <- (mean^2) / variance
    rate <- mean / variance
  }
  if (!is.null(scale)) {
    rate <- 1 / scale
  }
  if (is.null(shape) || is.null(rate)) {
    stop("Must provide either (shape and rate) or (mean and variance)")
  }
  cutoffs <- sort(cutoffs)
  n_intervals <- length(cutoffs) - 1
  expectations <- numeric(n_intervals)
  for (i in 1:n_intervals) {
    lower <- cutoffs[i]
    upper <- cutoffs[i + 1]
    prob <- pgamma(upper, shape = shape, rate = rate) -
      pgamma(lower, shape = shape, rate = rate)
    integrand <- function(x) {
      x * dgamma(x, shape = shape, rate = rate)
    }
    integral_result <- integrate(
      integrand,
      lower = lower,
      upper = upper
    )
    expectations[i] <- integral_result$value / prob
  }
  expectations
}
