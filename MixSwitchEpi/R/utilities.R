#' Neaten up the dust2 output into a data.table
#'
#' Efficient helper function.
#'
#' @export
neaten_state <- function(system_used, dust_state, times_simulated) {
  unpacked_sys <- dust2::dust_unpack_state(system_used, dust_state)

  nclasses <- nrow(unpacked_sys$S)
  ntimes <- ncol(unpacked_sys$S)

  s_mat <- t(unpacked_sys$S)
  e_mat <- t(unpacked_sys$E)
  i_mat <- t(unpacked_sys$I)
  r_mat <- t(unpacked_sys$R)

  ncols_total <- nclasses * 5 + 5 + 1
  result_mat <- matrix(0, nrow = ntimes, ncol = ncols_total)

  s_cols <- nclasses + seq_len(nclasses)
  e_cols <- 2 * nclasses + seq_len(nclasses)
  i_cols <- 3 * nclasses + seq_len(nclasses)
  r_cols <- 4 * nclasses + seq_len(nclasses)

  result_mat[, s_cols] <- s_mat
  result_mat[, e_cols] <- e_mat
  result_mat[, i_cols] <- i_mat
  result_mat[, r_cols] <- r_mat

  for (j in seq_len(nclasses)) {
    result_mat[, j] <- s_mat[, j] + e_mat[, j] + i_mat[, j] + r_mat[, j]
  }

  tot_idx <- 5 * nclasses + 1
  result_mat[, tot_idx] <- rowSums(s_mat) # s_tot
  result_mat[, tot_idx + 1] <- rowSums(e_mat) # e_tot
  result_mat[, tot_idx + 2] <- rowSums(i_mat) # i_tot
  result_mat[, tot_idx + 3] <- rowSums(r_mat) # r_tot

  # Compute total population as the sum of per-class totals (columns 1:nclasses).
  result_mat[, tot_idx + 4] <- rowSums(result_mat[, 1:nclasses, drop = FALSE])
  result_mat[, ncols_total] <- times_simulated

  result <- data.table::as.data.table(result_mat)

  col_names <- c(
    paste0("class", seq_len(nclasses)),
    paste0("S", seq_len(nclasses)),
    paste0("E", seq_len(nclasses)),
    paste0("I", seq_len(nclasses)),
    paste0("R", seq_len(nclasses)),
    "s_tot", "e_tot", "i_tot", "r_tot", "pop_tot", "time"
  )
  data.table::setnames(result, col_names)
  result
}

#' Simplifies epidemic output
#'
#' Especially useful when comparing epidemic trajectories across
#' homogeneous-inhomogeneous parameter sets
#'
#' @param classed_epidf epidemic with class structure
#' @return epidemic with S, E, I, R and time columns
#' @export
declass <- function(classed_epidf) {
  declassed <- classed_epidf |>
    dplyr::select(s_tot:time) |>
    dplyr::mutate(
      S = s_tot / pop_tot,
      E = e_tot / pop_tot,
      I = i_tot / pop_tot,
      R = r_tot / pop_tot
    ) |>
    dplyr::select(c(S, E, I, R, time))
  declassed
}


#' Writes out the epidemic to a csv
#'
#'
write_epidemic <- function(neatened_epidemic) {
  write.csv(neatened_epidemic,
    file = "epidemic_output.csv",
    row.names = FALSE
  )
}

#' kronecker_delta
#'
#' @param i first index
#' @param j second index
#' @return 1 if i == j, 0 otherwise
#' @export
kronecker_delta <- function(i, j) {
  ifelse(i == j, 1, 0)
}

#' Pretty-print contact and switching matrices
#'
#' @description Pretty-print the mixing and switching matrix
#' @param mix_mat the contact or mixing matrix describing intra-
#' and inter-activity class contatc frequencies
#' @param swtch_mat the matrix describing the switching of hosts between
#' @export
display_mix_switch <- function(mix_mat, swtch_mat) {
  cli::cli_h1("Mixing / Contact Matrix")
  print(mix_mat)
  cli::cli_h1("Activity Class Switching Matrix")
  print(swtch_mat)
}


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
