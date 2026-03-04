# generative models for contact (mixing) matrices and switching matrices.
#
# Using single or dual parameters, the matrix of required dimensionality
# (dictated by the number of activity classes) can be generated.
#

#' parameterically generate mixing matrix for the desired
#' number ofactivity classes
#'
#' Produces a symmetric matrix of dimension n_activity x n_activity. The matrix
#' needs to be multiplied by the state vectors S(t) and I(t) in the
#' transmission term. Since the each entry in matrix m: m i,j is per capita
#' susceptible class per capita infectious class, then
#' the rows of m sum to 1, making it a stochastic matrix.
#'
#' @param n_activity integer: the number of activity classes to account for
#' @param assortativity float: the degree to which classes mix with themselves
#' preferentially. when assortativity = 1 the function will
#' return the identity matrix
#' @param act_vector a vector containing the scale of the activity classes
#' @family mixing_models
#' @export
garnett_mixing <- function(n_activity,
                           epsilon_assort,
                           act_vector,
                           class_size_vector) {
  tot_act_size <- sum(act_vector * class_size_vector)

  # Every row i initially has the same base: (1 - epsilon) * (AjNj / sum(AkNk))
  row_base <- (1 - epsilon_assort) * (act_vector * class_size_vector) / tot_act_size

  # Fill the matrix row-by-row
  m <- matrix(row_base, nrow = n_activity, ncol = n_activity, byrow = TRUE)

  # Add epsilon assortativity to the diagonal entries where i == j
  diag(m) <- diag(m) + epsilon_assort

  m
}
#' Parametrically generate mixing matrix for the desired
#'
#' Note that this has an extra parameter compared to the other mixing
#' models. Thus the call signature needs to be adjusted accordingly
#'
#' @param n_activity the number of activity classes in the model
#' @param power_param the power parameter. if = 1
#' @param epsilon_assort the assortativity.
#' @param act_vector the vector of activity class scores.
#' Gives a measure of the activity of each class.
#' @param class_size_vector the vector of activity class sizes
#' @export
polynomial_mixing <- function(n_activity,
                              power_param,
                              epsilon_assort,
                              act_vector,
                              class_size_vector) {
  summed_prod <- sum(act_vector * class_size_vector)
  unnorm <- outer(
    1:n_activity,
    1:n_activity,
    function(i, j) {
      MixSwitchEpi::kronecker_delta(i, j) * epsilon_assort +
        (1 - epsilon_assort
        ) * ((abs(i - j) /
          (n_activity - 1))^power_param) * (act_vector[j]
        * class_size_vector[j]
          / summed_prod)
    }
  )
  unnorm <- as.matrix(unnorm)
  # row-wise normalisation: sum_j rho_{ij} = 1
  rho <- sweep(unnorm, 1, rowSums(unnorm), FUN = "/")
  rho
}


#' Parameterically generate exponential mixing model
#'
#' Automatically computes the normalisation
#' constants needed
#' @param n_activity the number of activity classes
#' @param epsilon_assort the assortativity of the contact pattern
#' @param act_vector the vector of activity class scores
#' @param class_size_vector the vector of activity class sizes
#' @export
exponential_mixing <- function(n_activity,
                               epsilon_assort,
                               act_vector,
                               class_size_vector) {
  summed_prod <- sum(act_vector * class_size_vector)
  unnorm <- outer(
    1:n_activity,
    1:n_activity,
    function(i, j) {
      (act_vector[j] * class_size_vector[j] / summed_prod) *
        exp(-epsilon_assort * abs(i - j))
    }
  )
  unnorm <- as.matrix(unnorm)
  # row-wise normalisation: sum_j rho_{ij} = 1
  rho <- sweep(unnorm, 1, rowSums(unnorm), FUN = "/")
  rho
}


#' Generate the activity-class scaler for each activity class
#' Linear activity scores
#'
#' Simple function which linearly interpolates between the lowest-activity and
#' the highest-activity class. Hence, for each class i we have a positive number
#' A_i which gives a measure of activity.
#'
#' (More complex models of how activity
#' changes with the activity class are possible.)
#' @export
activity_linear <- function(min_act, max_act, n_classes) {
  activity_vector <- seq(from = min_act, to = max_act, length.out = n_classes)
  activity_vector
}

#' Generate a null activity vector
#'
#' Function of use when no activity structure is desired.
#' Especially relevant for testing purposes.
#'
#' @param n_classes Number of activity classes
#' @returns A numeric vector of activity scores (all equal to 1)
#' @export
activity_null <- function(n_classes) {
  activity_vector <- rep(1, n_classes)
  activity_vector
}
