#'

# approximations of interest:
# RMSE as well as the best fit when the


#' Apperoximate a switching system with a no-switch system
#' using a single scalar epsilon parameter for assortativity/proportionate
#' mixing.
#'
#' @param reference_parameters A list of reference parameters.
#' @export
approx_epsilon <- function(reference_parameters,
                           initial_condition_rule,
                           pathogen_string) {
  n_activity <- reference_parameters$n_activity
  switch_rate <- reference_parameters$switch_rate


  model <- MixSwitchEpi::seirs_act()

  reference <- MixSwitchEpi::run_system_wparms(
    parms = reference_parameters,
    model = model
  )

  fit_res <- fit_noswitch(parms, model)
  approx_parms <- fit_res$best_parms

  cat(sprintf(
    "Best epsilon: %0.4f (cost=%0.6g)\n",
    fit_res$best_eps, fit_res$best_cost
  ))

  approx_sim <- run_system_wparms(parms = approx_parms, model = model)

  # # prepare data for plotting (normalize by N)
  # Ntot <- parms$N
  # df_ref <- data.frame(time = reference$time, i = reference$i_tot / Ntot, type = "reference")
  # df_approx <- data.frame(time = approx_sim$time, i = approx_sim$i_tot / Ntot, type = "approximation")
  # df_plot <- rbind(df_ref, df_approx)

  # p <- ggplot2::ggplot(df_plot, aes(x = time, y = i, color = type)) +
  #   ggplot2::geom_line(size = 0.8) +
  #   ggplot2::theme_minimal() +
  #   ggplot2::labs(
  #     title = "Reference (switching) vs No-switch approximation",
  #     x = "Time", y = "Prevalence (I / N)"
  #   ) +
  #   ggplot2::scale_color_manual(values = c("reference" = "#1f78b4", "approximation" = "#e31a1c"))

  # ggplot2::ggsave(
  #   filename = "noswitch_approx.svg",
  #   plot = p,
  #   width = 9,
  #   height = 4
  # )
  message("Saved plot to map_onto_noswitch_demo_R.png")
}

#' Approximate a switching system with a minimally constrained contact matrix
#' with no switching
#'
#'
#'
#' @export
#' @seealso approx_epsilon which approximates
#' with a single scalar epsilon parameterising
#' the distribution.
approx_minimally_constrained <- function() {}

#'
#'
#' Approximate a switching system with a no-switch system
#' using both minimally constrained contact matrix and
#' epsilon-constrained approximations.
#'
#' @export
approx_both <- function() {}

#' Characterise error in approximations across
#' parameter space
#'
#' @param gamma_range range of recovery rates to explore
#' @param switch_rate_range range of switching rates to explore
#'
#' @export
approximation_error <- function(gamma_range, switch_rate_range) {
  errors <- lapply(gamma_range, function(gamma_val) {
    lapply(switch_rate_range, function(switch_rate_val) {
      parms <- MixSwitchEpi::make_parms(
        n_activity = 5,
        switch_rate = switch_rate_val,
        initial_condition_rule = "uniform",
        pathogen_string = "sars-cov-2"
      )
      parms$gamma <- gamma_val

      model <- MixSwitchEpi::seirs_act()

      reference <- MixSwitchEpi::run_system_wparms(
        parms = parms,
        model = model
      )

      fit_res <- fit_noswitch(parms, model)

      approx_parms <- fit_res$best_parms

      approx_sim <- run_system_wparms(parms = approx_parms, model = model)

      # compute error metric (e.g., RMSE between reference and approx)
      error_metric <- compute_error(reference, approx_sim)

      list(
        gamma = gamma_val,
        switch_rate = switch_rate_val,
        error = error_metric
      )
    })
  })
}
