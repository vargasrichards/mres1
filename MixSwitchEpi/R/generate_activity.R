#' Retrieve fitted activity scores for a given setting
#'
#' This function retrieves the fitted activity scores for a specific setting
#' based on empirical data from (Pung et al., 2024) Interface
#' @param setting_string A string indicating the setting
#' @param n_classes The number of activity classes
#' @returns A numeric vector of activity scores
#' @examples
#' activity_empirical("haslemere", n_classes = 5)
#' @export
activity_empirical <- function(setting_string, n_classes) {
  if (setting_string == "haslemere") {
    # Taken from haslemere_M.csv scores or hardcoded in compare_curves.R
    scores <- c(0.3146078, 1.1479220, 2.3361036, 4.2211382, 9.3466695)
    if (n_classes != 5) {
      warning("Haslemere data usually expects 5 classes, returning 5.")
    }
    return(scores)
  } else if (setting_string == "cruise") {
    cli::cli_alert_danger("Cruise setting not yet implemented.")
    stop()
  } else {
    stop("Unknown setting string")
  }
}
