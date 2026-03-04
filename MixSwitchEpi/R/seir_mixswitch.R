#  this file contains code for the demonstration figure
# in main_f()

#' Run the system and produce neat output
#'
#' Bare-bones function to simulate the system using dust.
#' can call this many times if needed, for eg optimisation
#'
#' @param parms parameters for the model
#' @param model odin model
#' @family simulation
#' @export
run_system_wparms <- function(parms, model) {
  t <- seq(0,
    parms$t_end,
    length.out = as.integer(parms$t_end * 5)
  ) # can be adjusted for finer res.
  seirs_activity <- dust2::dust_system_create(model, parms,
    ode_control = dust2::dust_ode_control(max_steps = 1e7)
  )
  dust2::dust_system_set_state_initial(seirs_activity)
  sys <- dust2::dust_system_simulate(seirs_activity, t)
  neatsys <- MixSwitchEpi::neaten_state(
    system_used = seirs_activity,
    dust_state = sys,
    times_simulated = t
  )
  neatsys
}

#' Demo function for plotting
#'
#'
#' @export
main_f <- function(target_r0) {
  n_activity <- 5
  N <- 1e5
  epsilon_assort <- 0.5
  cmat <- MixSwitchEpi::garnett_mixing(
    n_activity = n_activity,
    epsilon_assort = epsilon_assort,
    act_vector = 1:n_activity,
    class_size_vector = rep(N / n_activity, n_activity)
  )
  some_switching <- MixSwitchEpi::simple_switching(
    switch_rate = 0.2,
    num_classes = n_activity
  )
  default_parms <- list(
    m = cmat,
    epsilon = epsilon_assort,
    n_activity = n_activity,
    S0 = rep(1 / n_activity, n_activity),
    E0 = c(0.01, 0.01, 0.01, 0.01, 0.01),
    beta = 0.05,
    sigma = 0.1,
    gamma = 0.05,
    N = N,
    t_end = 500,
    omega = 0,
    swch = some_switching,
    class_sizes = rep(N / n_activity, n_activity),
    activity_scores = 2 * (1:n_activity)
  )

  cal_parms <- calibrate_parms(
    preliminary_parms = default_parms,
    target_r0 = target_r0
  )
  t <- seq(0, cal_parms$t_end, length.out = cal_parms$t_end * 6)
  seirs <- dust2::dust_system_create(seirs_act(), pars = cal_parms)
  dust2::dust_system_set_state_initial(sys = seirs)
  sys <- dust2::dust_system_simulate(seirs, t)
  dust2::dust_unpack_state(seirs, state = sys)

  model_output <- MixSwitchEpi::neaten_state(seirs, sys, t)
  model_output_rt <- MixSwitchEpi::compute_rt_ts(model_output, cal_parms)
  metrics <- MixSwitchEpi::compute_all(model_output,
    classed_detail = FALSE,
    prevalence_threshold = 0.01,
    params = cal_parms
  )
  # p1 <- MixSwitchEpi::plot_epidemic_annotated(model_output, metrics)
  p3 <- MixSwitchEpi::plot_epidemic_annotated(model_output,
    metrics,
    compartments = c(
      "S",
      "E",
      "I",
      "R"
    )
  )
  # p4 <- MixSwitchEpi::plot_epidemic_faceted(model_output, metrics)
  rtp <- MixSwitchEpi::plot_epidemic_with_rt(model_output_rt,
    metrics,
    compartments = c("I", "E"), rt_only = T
  )
  p3 <- p3 + ggplot2::scale_x_continuous(expand = expansion(mult = c(0, 0)))
  rtp <- rtp + ggplot2::scale_x_continuous(expand = expansion(mult = c(0, 0)))
  ppanel <- patchwork::wrap_plots(p3, rtp, nrow = 2, ncol = 1, axes = "collect")
  ppanel <- ppanel + patchwork::plot_annotation(tag_levels = "A")

  # ggplot2::ggsave(ppanel, "panel_example.pdf")
  list(
    model_output = model_output,
    model_output_rt = model_output_rt,
    metrics = metrics,
    plots = list(
      p3 = p3,
      rtp = rtp,
      combi = ppanel
    )
  )
}
