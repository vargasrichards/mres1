#' Deterministic SEIR(S) epidemic model w/ activity classes and switching
#' @description
#' The main deterministic model framework used. Comprises a
#' Susceptible-Exposed-Infectious-Recovered-Suscetpible (SEIRS model) with
#' a user-specified nunber of activity classes which inidividuals can move
#' between at rates specified by a 'switching matrix'. This evolves determini-
#' stically according to a system of ODEs
#' @export

# the switching
# matrix is applied to
# the activity class state occupancy vector.
switch_s_tmp[, ] <- S[i] * swch[i, j] # need to check whether i or j
switch_s[] <- sum(switch_s_tmp[, i])
switch_e_tmp[, ] <- E[i] * swch[i, j]
switch_e[] <- sum(switch_e_tmp[, i])
switch_i_tmp[, ] <- I[i] * swch[i, j]
switch_i[] <- sum(switch_i_tmp[, i])
switch_r_tmp[, ] <- R[i] * swch[i, j]
switch_r[] <- sum(switch_r_tmp[, i])
deriv(S[]) <- -lambda[i] * S[i] + omega * R[i] + switch_s[i]
deriv(E[]) <- lambda[i] * S[i] - sigma * E[i] + switch_e[i]
deriv(I[]) <- sigma * E[i] - gamma * I[i] + switch_i[i]
deriv(R[]) <- gamma * I[i] - omega * R[i] + switch_r[i]

# s_ij calculation:
# i is source class (from I[i]), j is target class.
# i is target class, j is source class.
s_ij[, ] <- m[i, j] * I[j] / class_sizes[j]
lambda[] <- beta * activity_scores[i] * sum(s_ij[i, ])

# Initialise compartments using per-class sizes (class_sizes), not total N.
# S0 and E0 are fractions per class, so multiply by class_sizes[i] to get counts.
initial(S[]) <- S0[i] * class_sizes[i]
initial(E[]) <- E0[i] * class_sizes[i]
initial(I[]) <- 0
initial(R[]) <- 0


S0 <- parameter()
E0 <- parameter()
# N <- parameter()
m <- parameter() # mixing matrix: defining the
swch <- parameter() # switching matrix
beta <- parameter() # scalar transmission parameter.
sigma <- parameter() # going from E to I (progression to symptomaticity)
gamma <- parameter() # recovery rate (I -> R)
omega <- parameter() # parameterises our immune waning rate

n_activity <- parameter() # the number of activity classes
class_sizes <- parameter()
activity_scores <- parameter()

dim(S, S0, E0, lambda, R, I, E, class_sizes, activity_scores) <- n_activity
dim(swch, m, s_ij, switch_s_tmp, switch_e_tmp, switch_i_tmp, switch_r_tmp) <- c(n_activity, n_activity)
dim(switch_s, switch_e, switch_i, switch_r) <- n_activity
