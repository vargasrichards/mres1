##' @useDynLib MixSwitchEpi, .registration = TRUE
##' @importFrom odin2 odin
##' @importFrom dust2 dust_system_create dust_system_set_state_initial dust_system_simulate
##' @importFrom data.table :=
##' @importFrom data.table .SD
##' @importFrom data.table .N
NULL

cache <- new.env()

# Avoid R CMD check notes for data.table/dplyr columns
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      "s_tot", "e_tot", "i_tot", "r_tot", "ei_tot", "ei", "pop_tot",
      "S", "E", "I", "R", "time", "class", "final_size", "peak_time", "peak_size",
      "residence_time", "calib_at_r0", "metric_label", "metric_code",
      "hom_value", "x_label_pos", "infec_state", "model_name",
      "res_time_plot", "res_fac", "value", "infection_state", "metric",
      "S0", "E0", "I0", "R0", "n_activity", "class_sizes", "activity_scores",
      "model_type", "hit_hom_pred", "fs_hom_pred", "hom_pred",
      "inf_value", "legend_label",
      "prevalence", "scenario", "res_label"
    )
  )
}
