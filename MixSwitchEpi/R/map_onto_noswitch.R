# numerically finds the closest non-switching system for given nonzero
# switching and a mixing/contact matrix
# A. Vargas Richards, Oct 2025. Imperial College London, DIDE
# this implementation is very basic. significant scope for optimisation.
# proof-of-concept script.

#' Make a vector for input into optimisation fun.
#' @param M the matrix from which we get the x for optimisation
#' @export
vec_from_mat <- function(M) {
  M[upper.tri(M, diag = TRUE)]
}

#' Make a matrix suitable for use as a contact matrix from a supplied vector
#' @param vec_opt the vector used by the optimiser
#' @param ndim The number of activity classes used
#' @export
mat_from_vec <- function(vec_opt, ndim) {
  M <- matrix(0, ndim, ndim)
  M[upper.tri(M, diag = TRUE)] <- vec_opt
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  M
}

#' Make an initial guess of an 'effective' mixing matrix
#'
#' Function takes the mixing and switching matrices to be
#' reproduced and applies an operator to them to make
#' the 'effective' mixing matrix. That is, the matrix that if
#' @param mix_mat The contacts between different activity classes (and within)
#' those activity classes.
#' @param switch_mat Describes the transitions between the different activity
#' classes.
#' @param parameterisation The parameterisation (other than the mix/switch)
#' matrices. So, recovery rate, transmission rate, etc.
#' @export
make_initial_guess <- function(mix_mat, switch_mat, parameterisation) {
  # we need the parameterisation to get the death rate
  # this initial guess is based on Markov chain theory.
  v <- vec_from_mat(mix_mat)
  print(glue::glue("initial guess for the effective mixing {v}"))
  v
}


#' Fit non-switching system
#'
#' Produce a fitted time-invariant system which is the 'best' approximation
#' to a time-varying epidemic system.
#'
#' This function produces, for a given parameterisation of the activity class
#' switching matrix, the mixing/contact matrix producing the
#'  most similar results for the
#' Exposed compartment.
#' @param parms a parameterisation for the modelled system
#' which includes switching rate > 0 (can also be zero for testing)
#' @param model odin model used for simulation
#' @param eps_init initial value of epsilon which is used for the reference
#' system.
#' @param metric the way
#' @returns sol solution of the optimisation problem
#' @export
fit_noswitch <- function(parms, model, eps_init = NULL,
                         metric = c("i_tot", "e_tot"), normalize = TRUE,
                         verbose = TRUE, nl_opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 200, xtol_rel = 1e-6)) {
  metric <- match.arg(metric)

  n_classes <- parms$n_activity
  if (is.null(n_classes) || !is.numeric(n_classes)) stop("parms$n_activity must be present and numeric")

  if (verbose) cli::cli_alert_info("Simulating reference system (with switching)...")
  reference <- run_system_wparms(parms = parms, model = model)
  ref_ts <- reference[[metric]]
  if (normalize && !is.null(parms$N)) ref_ts <- ref_ts / parms$N

  eps0 <- if (!is.null(eps_init)) eps_init else if (!is.null(parms$epsilon)) parms$epsilon else 0.5

  obj_fn <- function(x) {
    eps <- x[1]
    cmat <- generate_mixingmat(
      n_activity = parms$n_activity,
      epsilon_assort = eps,
      act_vector = parms$activity_scores,
      class_size_vector = parms$class_sizes
    )
    cand_parms <- parms
    cand_parms$m <- cmat
    cand_parms$swch <- matrix(0, nrow = n_classes, ncol = n_classes)

    sim <- run_system_wparms(
      parms = cand_parms,
      model = model
    )

    sim_ts <- sim[[metric]]
    sim_ts <- sim_ts / parms$N

    error <- sqrt(mean((sim_ts - ref_ts)^2))
    error
  }

  lb <- # we're using bounds here as ε must be between 0 and 1 by definition
    ub <- 1
  nlres <- nloptr::nloptr(
    x0 = c(eps0),
    eval_f = obj_fn,
    lb = c(lb),
    ub = c(ub),
    opts = nl_opts
  )

  best_eps <- nlres$solution[1]
  best_cost <- nlres$objective
  best_matrix <- MixSwitchEpi::generate_mixingmat(
    n_activity = parms$n_activity,
    epsilon_assort = best_eps,
    act_vector = parms$activity_scores,
    class_size_vector = parms$class_sizes
  )
  best_parms <- parms
  best_parms$m <- best_matrix
  best_parms$swch <- matrix(0, nrow = n_classes, ncol = n_classes)

  if (verbose == TRUE) {
    cli::cli_alert_success("Optimisation finished: best epsilon = {best_eps},
                                      cost = {round(best_cost, 6)}")
  }
  list(
    best_eps = best_eps,
    best_matrix = best_matrix,
    best_parms = best_parms,
    best_cost = best_cost,
    nl_result = nlres,
    reference = reference,
    metric = metric
  )
}


#' Score difference between model outputs
#'
#' @description Produces a score of the difference in model result
#' between the original mixing/switching parameterisation (`true_matrices`)
#' and the candidate mixing matrix.
#' @param x the vector used by the optimiser to represent the mixing matrix.
#' @param reference_pattern gives the desired pattern of e.g., E(t) through
#' time, which we're going to try and match as close as poss. with no switching.
#' @param shared_parms the shared parameters (e.g., transmission rate, recovery
#' rate).
#' @param model the odin/monty model
#' @param similar_matrices bool. Indicates whether to
#'  reward effective contact matrices
#' which are similar to the original contact matrix. Experimental feature.
#' @export
score_difference <- function(x,
                             reference_pattern,
                             shared_parms,
                             model,
                             similar_matrices) {
  mixmat_forscore <- mat_from_vec(x, ndim = nrow(shared_parms$m))
  # no switching for the candidate
  smat <- matrix(0, nrow = nrow(mixmat_forscore), ncol = nrow(mixmat_forscore))
  shared_parms$swch <- smat
  shared_parms$m <- mixmat_forscore

  # run the system with candidate parameters and compute difference
  cand <- run_system_wparms(parms = shared_parms, model = model)
  # compute difference on total exposed (or total infectious) if available
  if (!("e_tot" %in% names(cand)) || !("e_tot" %in% names(reference_pattern))) {
    # fall back to i_tot
    diff_ts <- cand$i_tot - reference_pattern$i_tot
  } else {
    diff_ts <- cand$e_tot - reference_pattern$e_tot
  }
  sc <- sum(abs(diff_ts))
  sc
}

#' Make a generator matrix
#'
#' Makes a generator for the switching process.
#'
#' @param prelim_matrix A preliminary matrix.
#' @param diag_element The diagonal element value.
#' @export
make_generator <- function(prelim_matrix, diag_element) {
  for (i in seq_len(ncol(prelim_matrix))) {
    for (j in seq_len(ncol(prelim_matrix))) {
      if (i == j) {
        prelim_matrix[i, j] <- diag_element
      }
    }
  }
  prelim_matrix
}


#' Function to plot the approximated system against the original system
#'
#' Allows for immediate visual comparison of the trajectory of the
#' possible epidemic, both with and without movemement between activity
#' classes. This function plots the trajectories next to each other.
#'
#' It also computes the
#'
#' @param shared_parms the parameters which are shared by both the approxi-
#' mate system and the exact system. Usually this is just the  scalar parameters
#' @param original_parms the original parameters: usually nonzero switching matrix
#' and the contact matrix.
#' @param approximate_parms the optimised parameters giving similar epidemic
#' trajectories as the original parameters
#'
#' @export
compare_approx <- function(ref_parms, approx_parms, model) {
  ref_run <- run_system_wparms(parms = ref_parms, model = model)
  approx_run <- run_system_wparms(parms = approx_parms, model = model)

  # align and compute difference on total infectious
  if (!("i_tot" %in% names(ref_run)) || !("i_tot" %in% names(approx_run))) {
    stop("Expected 'i_tot' column in simulation outputs")
  }
  diff_runs <- data.frame(
    time = ref_run$time,
    ref_i = ref_run$i_tot,
    approx_i = approx_run$i_tot,
    diff = ref_run$i_tot - approx_run$i_tot
  )
  diff_runs
}
