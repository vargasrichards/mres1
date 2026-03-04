#' Convienience function to compute the spectral radius of a matrix
#'
#' The spectral radius is the maximum absolute value of the eigenvalues of the matrix.
#' @param a matrix
#' @returns the spectral radius (max abs eigenvalue)
#' @export
#' @seealso `compute_r0`, `compute_rt`
spectral_radius <- function(a) {
  eigenvals <- eigen(a, only.values = TRUE)$values
  max(abs(eigenvals))
}


#' Generates the rho matrix from parameters
#'
#' The matrix of mixing probabilities between activity classes. This is the matrix of
#' probabilities that an individual in activity class i has contact with an
#' individual in activity class j.
#' See the documentation for `rho_entry` for details on how this is constructed.
#'
#' @param params the model parameters
#' @export
rho_matrix <- function(params) {
  # If an explicit mixing matrix `m` is supplied in params, use it.
  if (!is.null(params$m)) {
    return(params$m)
  }

  # Otherwise construct a default mixing matrix from parameters.
  # Default to the garnett mixing model when epsilon/activity/class_sizes are provided.
  n_activity <- params$n_activity
  epsilon <- params$epsilon %||% 0
  act_vector <- params$activity_scores
  class_size_vector <- params$class_sizes

  if (is.null(act_vector) || is.null(class_size_vector) || is.null(n_activity)) {
    stop("rho_matrix: insufficient parameters to construct mixing matrix (need n_activity, activity_scores, class_sizes)")
  }

  mmat <- garnett_mixing(
    n_activity = n_activity,
    epsilon_assort = epsilon,
    act_vector = act_vector,
    class_size_vector = class_size_vector
  )
  mmat
}

#' Generates a single entry of the rho matrix
#'
#' @param i the row index
#' @param j the column index
#' @param params the model parameters
#' @export
rho_entry <- function(i, j, params) {
  m <- rho_matrix(params)
  m[i, j]
}

#' Generates a single entry of the T matrix
#'
#' The rate of transmission from class j to class i
#' @seealso `t_matrix`
#' @export
t_entry <- function(i, j, params) {
  # t[i,j] = beta * activity_scores[j] * m[j, i]
  m <- rho_matrix(params)
  params$beta * params$activity_scores[j] * m[j, i]
}

#' Generates the T matrix from parameters
#'
#' This T matrix describes the transmission component of the
#' NGM. See Mathematical Tools for Understanding Infectious Disease Dynamics,
#' Diekmann et al  for details.
#' @param params the model parameters
#' @seealso `t_entry`
#' @export
t_matrix <- function(params) {
  n_activity <- params$n_activity
  m <- rho_matrix(params)
  # t_entry logic: beta * activity_scores[j] * m[j, i]
  # T[i,j] = beta * t(M)[i,j] * act[j]
  tmat <- params$beta * t(m)
  tmat <- sweep(tmat, 2, params$activity_scores, "*")
  tmat
}

#' Makes the sigma matrix
#'
#' Constructs the matrix of the state-transitions.
#'
#'
#' @param params the model parameters
#' @family r0_calc
#' @export
sigma_matrix <- function(params) {
  n_activity <- params$n_activity
  # the precomputed switch matrix must be provided by callers
  if (is.null(params$swch)) stop("sigma_matrix: 'swch' (switch matrix) must be provided in params")
  sigma_mat <- params$swch
  gamma <- params$gamma
  for (i in 1:n_activity) {
    sigma_mat[i, i] <- sigma_mat[i, i] - gamma
  }
  sigma_mat
}

#' Constructs the next-generation matrix
#'
#' Make the next-generation matrix (NGM) for the specified model. This
#' uses the usual next-generation approach for compartmental models with a
#' large domain.
#' @param pars The model parameters (list)
#' @family r0_calc
#' @export
compute_ngm <- function(pars) {
  t_mat <- t_matrix(pars)
  sigma_mat <- sigma_matrix(pars)
  k <- -t_mat %*% solve(sigma_mat)
  k
}

#' Computes the next-generation matrix and R0
#'
#' @param pars the parameters for the model run
#' @family r0_calc
#' @export
compute_r0 <- function(pars) {
  k <- compute_ngm(pars)
  r0 <- spectral_radius(k)
  r0
}

#' spectral gap calculation
#'
#' compute the difference in absolute value between the
#' dominant and second-most dominant eigenvals.
#' Should provide some measure of the time of duration
#' transient behaviour.
#'
#' @param pars parameters for model run
#' @family epidemic_metrics
#' @export
compute_spectral_gap <- function(pars) {
  k <- compute_ngm(pars)
  eigvals <- eigen(k, only.values = TRUE)$values
  magnitudes <- sort(Mod(eigvals), decreasing = TRUE)
  sgap <- magnitudes[1] - magnitudes[2]
  sgap
}


#' Computes the Rt using the population state at a particular time t
#'
#' @param pars the model parameters
#' @param popstate the state of the population at the time for Rt.
#' Must have a list element $s_sizes giving the number of people
#' who are susceptible to infection in each activity class. This
#' will then be divided through by the vector of activity class
#' sizes.
#' @export
#' @family epidemic_metrics
#' @seealso `compute_r0`
compute_rt <- function(pars, popstate) {
  class_sizes <- pars$class_sizes
  n_activity <- pars$n_activity
  sigma_mat <- sigma_matrix(pars)
  t_mat <- t_matrix(pars)
  k <- -t_mat %*% solve(sigma_mat)
  # popstate$s_sizes may either be absolute counts or already scaled
  # (fractions in [0,1]). Detect which case to avoid double-scaling.
  s_sizes_vec <- popstate$s_sizes
  # coerce list-like or data.frame columns to numeric vector
  if (is.list(s_sizes_vec) || is.data.frame(s_sizes_vec)) {
    s_sizes_vec <- as.numeric(unlist(s_sizes_vec))
  } else {
    s_sizes_vec <- as.numeric(s_sizes_vec)
  }

  if (max(s_sizes_vec, na.rm = TRUE) <= 1 + 1e-8) {
    sfrac <- diag(s_sizes_vec,
      nrow = n_activity,
      ncol = n_activity
    )
  } else {
    sfrac <- diag(s_sizes_vec / class_sizes,
      nrow = n_activity,
      ncol = n_activity
    )
  }

  k_eff <- sfrac %*% k
  rt <- spectral_radius(k_eff)
  rt
}

#' Convenience function to get population info from model output
#'
#' Gets number of activity classes and total population size.
#' @param model_output Model output data.table
#' @returns info, a list containing relevant information
#' @family epidemic_info
#' @export
get_info <- function(model_output) {
  pop_total <- max(model_output[["pop_tot"]], na.rm = TRUE)
  stopifnot(length(pop_total) == 1)

  num_classes <- max(as.numeric(gsub(
    "[^0-9]", "",
    names(model_output)
  )), na.rm = TRUE)

  class_sizes <- integer(num_classes)
  for (i in seq_len(num_classes)) {
    class_label <- paste0("class", i)
    class_sizes[i] <- unique(model_output[[class_label]])[1]
  }
  stopifnot(isTRUE(all.equal(sum(class_sizes), pop_total, tolerance = 1e-8)))

  list(
    pop_total = pop_total,
    num_classes = num_classes,
    class_size_all = class_sizes
  )
}

#' Characterise epidemic peaks for total and by class
#'
#' This takes in scaled epidemic output and computes peak sizes.
#'
#' @param model_output Model output (scaled)
#' @param info Info list from get_info()
#' @param classed_detail If TRUE, compute for every activity class
#' @family epidemic_metrics
#' @export
characterise_peaks <- function(model_output,
                               info,
                               classed_detail = TRUE,
                               pars) {
  if (!inherits(model_output, "data.table")) {
    model_output <- data.table::as.data.table(model_output)
  }
  pop_size <- info$pop_total
  num_classes <- info$num_classes
  class_sizes <- info$class_size_all

  overall_peak <- MixSwitchEpi::find_peaks(
    epi_data = model_output,
    e_label = "e_tot",
    i_label = "i_tot",
    pop_size = pop_size,
    class = -1, pars = pars
  )

  if (classed_detail) {
    class_results <- vector("list", num_classes)

    for (i in seq_len(num_classes)) {
      class_size <- class_sizes[i]
      susceptibles <- paste0("S", i)
      exposeds <- paste0("E", i)
      infecteds <- paste0("I", i)

      class_results[[i]] <- find_peaks(
        epi_data = model_output,
        s_label = susceptibles,
        e_label = exposeds,
        i_label = infecteds,
        pop_size = class_size,
        class = i, pars = pars
      )
    }
    do.call(rbind, c(list(overall_peak), class_results))
  } else {
    overall_peak
  }
}

#' Compute how synchronised epidemic peaks are
#'
#'
#' @param peak_data Data frame of peak times per class
#' @export
peak_synchronicity <- function(peak_data) {}

#' Scale results by occupancy fraction for each activity class
#'
#' @param model_output Output from dust/odin run
#' @family epidemic_metrics
#' @export
scale_results <- function(model_output) {
  if (!inherits(model_output, "data.table")) {
    model_output <- data.table::copy(data.table::as.data.table(model_output))
  } else {
    model_output <- data.table::copy(model_output)
  }

  epi_info <- MixSwitchEpi::get_info(model_output)
  class_sizes <- epi_info$class_size_all
  num_classes <- epi_info$num_classes
  pop_tot <- epi_info$pop_total

  model_output[, `:=`(
    s_tot = s_tot / pop_tot,
    e_tot = e_tot / pop_tot,
    i_tot = i_tot / pop_tot,
    r_tot = r_tot / pop_tot
  )]

  for (i in seq_len(num_classes)) {
    susceptibles <- paste0("S", i)
    exposeds <- paste0("E", i)
    infecteds <- paste0("I", i)
    removeds <- paste0("R", i)
    this_class_size <- class_sizes[i]

    model_output[, (susceptibles) := get(susceptibles) / this_class_size]
    model_output[, (exposeds) := get(exposeds) / this_class_size]
    model_output[, (infecteds) := get(infecteds) / this_class_size]
    model_output[, (removeds) := get(removeds) / this_class_size]
  }

  model_output
}

#' Compute time to threshold prevalence
#'
#' @param model_output Output from dust/odin
#' @param prevalence_threshold Numeric threshold for (E + I) / N_j
#' @param num_classes Number of activity classes
#' @param classed_detail Whether to compute per class
#' @family epidemic_metrics
#' @export
compute_ttprev <- function(model_output, prevalence_threshold, num_classes,
                           classed_detail = TRUE) {
  stopifnot(prevalence_threshold >= 0 && prevalence_threshold <= 1)

  if (!inherits(model_output, "data.table")) {
    model_output <- data.table::as.data.table(model_output)
  }

  # Overall time to prevalence
  model_output[, ei_tot := i_tot + e_tot]
  overall_ttp_rows <- model_output[ei_tot >= prevalence_threshold]
  ttp_overall <- min(overall_ttp_rows[["time"]], na.rm = TRUE)

  overall_ttp_results <- data.frame(
    time_to_prev = ttp_overall,
    class = -1,
    prev_thresh = prevalence_threshold
  )

  if (classed_detail) {
    class_results <- vector("list", num_classes)

    for (i in seq_len(num_classes)) {
      exposeds <- paste0("E", i)
      infecteds <- paste0("I", i)

      model_output[, ei := get(exposeds) + get(infecteds)]
      ttp_rows <- model_output[ei >= prevalence_threshold]
      ttp_class <- min(ttp_rows[["time"]], na.rm = TRUE)

      class_results[[i]] <- data.frame(
        time_to_prev = ttp_class,
        class = i,
        prev_thresh = prevalence_threshold
      )
    }

    do.call(rbind, c(class_results, list(overall_ttp_results)))
  } else {
    overall_ttp_results
  }
}

#' Compute epidemic duration and final size
#'
#'
#'
#' @param model_output Model output
#' @param classed_detail Whether to compute per-class results
#' @family epidemic_metrics
#' @export
characterise_end <- function(model_output, classed_detail = TRUE) {
  if (!inherits(model_output, "data.table")) {
    model_output <- data.table::as.data.table(model_output)
  }

  info <- get_info(model_output)
  num_classes <- info$num_classes
  pop_total <- info$pop_total
  overall_start <- 0
  # Determine epidemic end as the first time after the epidemic peak
  # when the number of infecteds falls below 1 (i_tot < 1).
  # Avoid picking up the trivial early-time case when i_tot starts at 0.
  if (nrow(model_output) == 0) {
    overall_end <- NA
    whole_epiduration <- NA
    overall_fs <- NA
  } else {
    peak_idx <- which.max(model_output[["i_tot"]])
    peak_time <- model_output[peak_idx, ][["time"]]
    # look for first time after peak where i_tot < 1
    post_peak_rows <- model_output[model_output[["time"]] > peak_time & model_output[["i_tot"]] < 1, ]
    if (nrow(post_peak_rows) == 0) {
      overall_end <- NA
      whole_epiduration <- NA
      overall_fs <- NA
    } else {
      overall_end <- min(post_peak_rows[["time"]], na.rm = TRUE)
      whole_epiduration <- overall_end - overall_start
      min_idx <- which.min(post_peak_rows[["time"]])
      overall_fs <- post_peak_rows[min_idx, ][["r_tot"]] / pop_total
    }
  }
  overall_dets <- data.frame(
    class = -1,
    start_time = overall_start,
    end_time = overall_end,
    duration = whole_epiduration,
    final_size = overall_fs
  )

  if (classed_detail) {
    class_results <- vector("list", num_classes)

    for (i in seq_len(num_classes)) {
      class_label <- paste0("class", i)
      exposeds <- paste0("E", i)
      infecteds <- paste0("I", i)
      removeds <- paste0("R", i)

      start_rows <- model_output[model_output[[exposeds]] > 0, ]
      if (nrow(start_rows) == 0) {
        start_time <- NA
      } else {
        start_time <- min(start_rows[["time"]], na.rm = TRUE)
      }

      end_rows <- model_output[model_output[[infecteds]] < 1 &
        model_output[[removeds]] > 0, ]
      if (nrow(end_rows) == 0) {
        end_time <- NA
        final_size <- NA
        duration <- NA
      } else {
        end_time <- min(end_rows[["time"]], na.rm = TRUE)
        min_idx <- which.min(end_rows[["time"]])
        end_row <- end_rows[min_idx, ]
        final_unscaled <- end_row[[removeds]]
        in_class <- end_row[[class_label]]
        final_size <- final_unscaled / in_class
        duration <- end_time - start_time
      }
      class_results[[i]] <- data.frame(
        class = i,
        start_time = start_time,
        end_time = end_time,
        duration = duration,
        final_size = final_size
      )
    }
    do.call(rbind, c(class_results, list(overall_dets)))
  } else {
    overall_dets
  }
}

#' Compute empirical herd immunity threshold (HIT)
#'
#' @description Computes the empirical herd immunity threshold (HIT). This is
#' defined as the proportion of the population that has been infected when the
#' effective reproduction number Rt crosses below 1 for the first time.
#'
#' @param model_output_rt Model output data.table
#' with Rt column already calculated
#' @returns A list with elements:
#' - HIT: the empirical herd immunity threshold
#' - HIT_time: the time when Rt crosses below 1
#' - Rt_at_HIT: the value of Rt at the HIT time (should be close to 1)
#' @family epidemic_metrics
#' @export
compute_hit <- function(model_output_rt) {
  # Ensure data.table and copy to avoid modifying user input
  if (!inherits(model_output_rt, "data.table")) {
    model_output_rt <- data.table::as.data.table(model_output_rt)
  } else {
    model_output_rt <- data.table::copy(model_output_rt)
  }
  if (!"Rt" %in% names(model_output_rt)) {
    stop("compute_hit: model_output_rt must contain an 'Rt' column")
  }
  idx <- which(model_output_rt$Rt < 1)
  if (length(idx) == 0) {
    # Rt never crossed below 1 (can occur at fast switching / short t_end).
    # Report only when debug mode is active; otherwise return NA silently so
    # the caller (compute_all) can substitute the analytic approximation.
    if (Sys.getenv("MIXSWITCH_DEBUG") == "1") {
      cli::cli_alert_danger(
        "compute_hit: Rt never drops below 1; HIT set to NA (analytic fallback will be used). \\
        Rt range: [{round(min(model_output_rt$Rt, na.rm=TRUE),4)}, \\
        {round(max(model_output_rt$Rt, na.rm=TRUE),4)}]"
      )
    }
    return(list(HIT = NA_real_, HIT_time = NA_real_, Rt_at_HIT = NA_real_))
  }
  crit <- model_output_rt[idx[1], ]
  if (!"s_tot" %in% names(crit)) {
    stop("compute_hit: model_output_rt must contain an 's_tot'
    column (susceptible total)")
  }
  s_tot_val <- crit$s_tot

  # If s_tot is in absolute counts (>1), scale using pop_tot if available
  if (!is.na(s_tot_val) && s_tot_val > 1 + 1e-8) {
    if ("pop_tot" %in% names(crit)) {
      s_frac <- s_tot_val / crit$pop_tot
    } else {
      stop("compute_hit: s_tot appears to be in counts
      but 'pop_tot' is missing;
      provide scaled s_tot or include pop_tot in model_output_rt")
    }
  } else {
    s_frac <- s_tot_val
  }

  hit <- 1 - s_frac
  stopifnot(
    "HIT must be between 0 and 1 - check for
    inappropriate scaling" = (hit >= 0 && hit <= 1)
  )
  hit_time <- crit$time
  rt_at_hit <- crit$Rt
  list(
    HIT = hit,
    HIT_time = hit_time,
    Rt_at_HIT = rt_at_hit
  )
}

#' Find epidemic peaks (for I or E + I)
#'
#' This function requires that the data are scaled.
#'
#' @param epi_data data.table with epidemic trajectories
#' @param s_label Column name for susceptibles
#' @param e_label Column name for exposed
#' @param i_label Column name for infected
#' @param pop_size Total (sub)population size
#' @param class Activity class identifier
#' @family epidemic_metrics
#' @export
find_peaks <- function(epi_data,
                       s_label,
                       e_label,
                       i_label,
                       pop_size,
                       class,
                       pars) {
  stopifnot(length(pop_size) == 1)

  if (!inherits(epi_data, "data.table")) {
    epi_data <- data.table::as.data.table(epi_data)
  }

  epi_data[, ei := get(e_label) + get(i_label)]
  peak_ei_idx <- which.max(epi_data[["ei"]])
  peak_ei <- epi_data[peak_ei_idx]

  peak_size_ei <- peak_ei[["ei"]]
  peak_time_ei <- peak_ei[["time"]]

  peak_i_idx <- which.max(epi_data[[i_label]])
  peak_i <- epi_data[peak_i_idx]
  peak_size_i <- peak_i[[i_label]]
  peak_time_i <- peak_i[["time"]]

  d <- data.frame(
    ptEI = peak_time_ei,
    ptI = peak_time_i,
    psEI = peak_size_ei,
    psI = peak_size_i,
    class = class
  )
  d
}

#' Wrapper function to compute all metrics
#'
#' @param model_output Model output (must contain Rt column or
#' params must be provided)
#' @param prevalence_threshold Prevalence threshold for ttprev
#' @param params Model parameters used to generate model output
#' @family epidemic_metrics
#' @export
compute_all <- function(model_output,
                        classed_detail = TRUE,
                        prevalence_threshold = NULL,
                        params) {
  if (!inherits(model_output, "data.table")) {
    model_output <- data.table::as.data.table(model_output)
  }

  # Ensure Rt column exists; if params provided we can compute Rt time series
  if (!"Rt" %in% names(model_output)) {
    if (is.null(params)) {
      stop("compute_all: either 'Rt' must be present in model_output or 'params' must be supplied to compute it")
    }
    model_output <- compute_rt_ts(model_output, params, add_to_output = TRUE)
  }

  # compute hit info (may return NA values if Rt never < 1)
  comphit <- compute_hit(model_output)
  hit_data <- comphit$HIT
  hit_time <- comphit$HIT_time
  rt_at_hit <- comphit$Rt_at_HIT

  epi_info <- get_info(model_output)

  end_details <- characterise_end(model_output, classed_detail)
  rescaled_results <- scale_results(model_output)

  peak_details <- characterise_peaks(
    rescaled_results,
    epi_info,
    classed_detail,
    params
  )

  result_dt <- data.table::as.data.table(end_details)
  peak_dt <- data.table::as.data.table(peak_details)

  data.table::setkey(result_dt, class)
  data.table::setkey(peak_dt, class)

  result <- result_dt[peak_dt]
  data.table::setorder(result, class)

  # r0 and spectral gap (require params)
  if (!is.null(params)) {
    r0 <- compute_r0(params)
    sgap <- compute_spectral_gap(params)
  } else {
    r0 <- NA_real_
    sgap <- NA_real_
  }

  # attach hit / r0 / spectral gap to overall row (class == -1)
  result[class == -1, `:=`(hit = hit_data, hit_time = hit_time, rt_at_hit = rt_at_hit, r0 = r0, spectral_gap = sgap)]

  # optionally compute time-to-prevalence if requested
  if (!is.null(prevalence_threshold)) {
    ttp <- compute_ttprev(model_output, prevalence_threshold, epi_info$num_classes, classed_detail = classed_detail)
    ttp_dt <- data.table::as.data.table(ttp)
    data.table::setkey(ttp_dt, class)
    # merge ttp into result by class
    result <- ttp_dt[result]
    data.table::setorder(result, class)
  }

  result
}


#' Compute R0 and Rt time series from model output
#'
#' Uses an effective NGM approach, accounting for the depletion of susceptibles
#'
#' @param model_output Model output from dust/odin
#' @param pars Model parameters
#' @param add_to_output If TRUE, adds R0 and Rt columns to model_output
#' @returns If add_to_output is TRUE, returns model_output with R0 and Rt
#' columns added. Otherwise, returns a data.frame with time, R0, Rt columns.
#' @export
compute_rt_ts <- function(model_output,
                          pars,
                          add_to_output = TRUE) {
  if (!inherits(model_output, "data.table")) {
    model_output <- data.table::as.data.table(model_output)
  } else {
    model_output <- data.table::copy(model_output)
  }

  n_activity <- pars$n_activity
  class_sizes <- pars$class_sizes

  t_mat <- t_matrix(pars)
  sigma_mat <- sigma_matrix(pars)
  k_mat <- -t_mat %*% solve(sigma_mat)

  r0 <- spectral_radius(k_mat)

  s_cols <- paste0("S", seq_len(n_activity))
  s_matrix <- as.matrix(model_output[, s_cols, with = FALSE])
  # Ensure numeric matrix (handle cases where columns are list-columns)
  s_matrix <- matrix(as.numeric(s_matrix), nrow = nrow(model_output), ncol = n_activity)
  # If S columns are already fractions (max <= 1)
  # we don't divide by class_sizes;
  # otherwise assume S columns are absolute counts and convert to fractions.
  if (max(s_matrix, na.rm = TRUE) <= 1 + 1e-8) {
    s_frac_matrix <- s_matrix
  } else {
    s_frac_matrix <- sweep(s_matrix, 2, class_sizes, "/")
  }

  n_times <- nrow(s_frac_matrix)
  rt_values <- numeric(n_times)
  for (t in seq_len(n_times)) {
    s_frac_t <- as.vector(s_frac_matrix[t, ])
    sfrac_diag <- diag(s_frac_t, nrow = n_activity, ncol = n_activity)
    k_eff <- sfrac_diag %*% k_mat
    rt_values[t] <- spectral_radius(k_eff)
  }

  # Optional debug output when MIXSWITCH_DEBUG env var is set
  if (Sys.getenv("MIXSWITCH_DEBUG") == "1") {
    cat("DEBUG: t_mat:\n")
    print(t_mat)
    cat("DEBUG: sigma_mat:\n")
    print(sigma_mat)
    cat("DEBUG: k_mat:\n")
    print(k_mat)
    cat("DEBUG: rt_values:\n")
    print(rt_values)
  }

  if (add_to_output) {
    # DEBUG: print rt_values for troubleshooting
    # cat("DEBUG rt_values:", paste(rt_values, collapse=","), "\n")
    data.table::set(model_output, j = "R0", value = r0)
    data.table::set(model_output, j = "Rt", value = rt_values)
    return(model_output)
  } else {
    return(data.frame(
      time = model_output$time,
      R0 = r0,
      Rt = rt_values
    ))
  }
}

#' Computes the outcome metrics for a homogeneous population
#'
#' That is, the outcome metrics for a population where
#' each individual has the same activity. This is equivalent
#' to a mass action single pop.
#'
#' @param base_params the basic parameters: beta, sigma, gamma
#' @export
compute_homogeneous_metrics <- function(base_params) {
  # first lets calculate the population mean activity
  # If callers haven't supplied a switch matrix ('swch'), construct a
  # trivial zero-switch matrix here for the homogeneous convenience wrapper.
  if (is.null(base_params$swch)) {
    n_activity <- base_params$n_activity %||% length(base_params$activity_scores)
    base_params$swch <- matrix(0, n_activity, n_activity)
  }
  K <- compute_ngm(base_params)

  # ensure activity_scores are constant (homogeneous assumption)
  assertthat::assert_that(
    is.numeric(base_params$activity_scores) &&
      all(base_params$activity_scores == base_params$activity_scores[1]),
    msg = "compute_homogeneous_metrics: activity_scores must be the same for all classes to compute homogeneous metrics"
  )

  homogeneous_r0 <- spectral_radius(K)
  homogeneous_hit <- 1 - (1 / homogeneous_r0)

  # finalsize::final_size expects a scalar R0
  homogeneous_fs <- finalsize::final_size(r0 = homogeneous_r0)$p_infected

  hom_metrics <- list(
    hom_r0 = homogeneous_r0,
    hom_hit = homogeneous_hit,
    hom_fs = homogeneous_fs
  )
  hom_metrics
}
