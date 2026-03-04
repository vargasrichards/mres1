#' Make the simplest possible switching matrix.
#'
#' User specifies the number of activity classes to be used and the
#' switching rate = 1 - rate(remain) in class i, conditional on not
#' recovering.
#' @param switch_rate rate of switching between activity classes
#' @param num_classes number of activity classes
#' @returns A square matrix of dimension num_classes x num_classes
#' @export
simple_switching <- function(switch_rate, num_classes) {
  w <- matrix(nrow = num_classes, ncol = num_classes)
  for (i in 1:num_classes) {
    for (j in 1:num_classes) {
      if (i == j) {
        w[i, j] <- -switch_rate
      } else {
        w[i, j] <- switch_rate / (num_classes - 1)
      }
    }
  }
  w
}

#' Switching matrix where you can only move between adjacent classes
#'
#' Hence there is twice expected residence time in the extreme classes.
#'
#'
#' @param switch_rate switcjing rate. in general = 1/expected residence time
#' for the non-extreme classes (ie, 2,3,4) for a 5-class model.
#' @param num_classes number of activity classes
#' @seealso [`simple_switching`] which uses uniform movements
#' where you can move more than distance of 1 class at once.
#' @export
adjacent_switching <- function(switch_rate, num_classes) {
  w <- matrix(0, nrow = num_classes, ncol = num_classes)
  diag(w) <- -switch_rate
  w[row(w) == col(w) - 1] <- switch_rate / 2
  w[row(w) == col(w) + 1] <- switch_rate / 2
  w[1, 1] <- -switch_rate / 2
  w[num_classes, num_classes] <- -switch_rate / 2
  w
}

#' Switching matrix where you can only move between adjacent classes but wraps lowest and highest activities.
#'
#' Note that we can have min activity -> max activity in a single step.
#'
#' @param switch_rate switching rate. in general = 1/expected residence time
#' for the non-extreme classes (ie, 2,3,4) for a 5-class model.
#' @param num_classes number of activity classes
#' @export
adjacent_switching_wrapped <- function(switch_rate, num_classes) {
  w <- matrix(0, nrow = num_classes, ncol = num_classes)
  diag(w) <- -switch_rate
  for (i in seq_len(num_classes)) {
    # target j = i-1 and i+1 with wrapping
    j_low <- if (i == 1) num_classes else i - 1
    j_high <- if (i == num_classes) 1 else i + 1

    w[i, j_low] <- switch_rate / 2
    w[i, j_high] <- switch_rate / 2
  }
  w
}

#' Elaborated switching
#'
#' Here we assume that w(i,j) is the flow from i -> j. Note this is different
#' to some other conventions.
#' Particularly, the ngm calculation using Sigma matrix.
#'
#'
#' @export
elaborated_switching <- function(switch_rate, num_classes, increment_scalar) {
  w <- matrix(nrow = num_classes, ncol = num_classes)

  for (i in 1:num_classes) {
    for (j in 1:num_classes) {
      switch_rate <- switch_rate * increment_scalar^(i - 1)
      if (i == j) {
        w[i, j] <- -switch_rate
      } else {
        w[i, j] <- switch_rate / (num_classes - 1)
      }
    }
  }
}


# HELPER #
#' Find the steady-state distribution of a continuous-time Markov chain
#' given its transition rate matrix W.
#'
#'
#' @param W A square matrix where W(i,j) is the rate of transition
#'
#' @export
find_steadystate <- function(W) {
  n <- nrow(W)

  # The steady state π satisfies: π^T W = 0 and sum(π) = 1
  # Equivalently: W^T π = 0 and sum(π) = 1

  # Build system of equations:
  # [W^T    ] π = [0]
  # [1 1 1...] π = [1]

  A <- rbind(t(W), rep(1, n))
  b <- c(rep(0, n), 1)

  # Solve the overdetermined system using least squares
  # (W^T has rank n-1, so we need the sum constraint)
  pi <- as.vector(solve(t(A) %*% A) %*% t(A) %*% b)

  pi <- pmax(pi, 0)
  pi <- pi / sum(pi)

  pi
}
