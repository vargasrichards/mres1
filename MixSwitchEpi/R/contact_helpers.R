## Helpers for generating mixing matrices

#' Generate a mixing/contact matrix with selectable model
#'
#' Wrapper to produce a mixing matrix. Defaults to the Garnett model.
#'
#' @param n_activity integer number of activity classes
#' @param epsilon_assort assortativity parameter
#' @param act_vector vector of activity scores
#' @param class_size_vector vector of class sizes
#' @param mixing_model character; one of 'garnett', 'polynomial', 'exponential'
#' @export
generate_mixingmat <- function(n_activity,
                               epsilon_assort,
                               act_vector,
                               class_size_vector,
                               mixing_model = c(
                                 "garnett",
                                 "polynomial",
                                 "exponential"
                               ),
                               power_param = NULL) {
  mixing_model <- match.arg(mixing_model)
  if (mixing_model == "garnett") {
    return(garnett_mixing(
      n_activity = n_activity,
      epsilon_assort = epsilon_assort,
      act_vector = act_vector,
      class_size_vector = class_size_vector
    ))
  } else if (mixing_model == "polynomial") {
    if (is.null(power_param)) stop("polynomial mixing requires power_param")
    return(polynomial_mixing(
      n_activity = n_activity,
      power_param = power_param,
      epsilon_assort = epsilon_assort,
      act_vector = act_vector,
      class_size_vector = class_size_vector
    ))
  } else if (mixing_model == "exponential") {
    return(exponential_mixing(
      n_activity = n_activity,
      epsilon_assort = epsilon_assort,
      act_vector = act_vector,
      class_size_vector = class_size_vector
    ))
  } else {
    stop("unrecognised mixing model")
  }
}
