library(testthat)
library(MixSwitchEpi)
library(finalsize)

# test_that("compute_homogeneous_metrics works", {
#   params <- list(beta = 0.4, gamma = 0.2)
#   res <- MixSwitchEpi::compute_homogeneous_metrics(params)

#   expect_equal(res$hom_r0, 2)
#   expect_equal(res$hom_hit, 0.5)
#   expect_true(res$hom_fs > 0.5)
# })

test_that("final size matches finalsize::final_size for homogeneous no-switch model", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("finalsize")

  n_activity <- 3
  N <- 1e6
  class_sizes <- rep(N / n_activity, n_activity)
  beta <- 0.5
  gamma <- 0.1
  R0 <- beta / gamma # 5

  # Analytical final size from finalsize package
  fs_analytical <- finalsize::final_size(r0 = R0)$p_infected

  # Simulate with proportionate mixing (homogeneity assumption) and no switching
  m <- matrix(1 / n_activity, n_activity, n_activity)
  pars <- list(
    n_activity = n_activity,
    S0 = rep(1 - 1e-4, n_activity),
    E0 = rep(1e-4, n_activity),
    N = N,
    beta = beta,
    sigma = 1 / 3,
    gamma = gamma,
    omega = 0,
    m = m,
    swch = matrix(0, n_activity, n_activity),
    class_sizes = class_sizes,
    activity_scores = rep(1, n_activity),
    t_end = 1000
  )

  sys <- dust2::dust_system_create(
    seirs_act(),
    pars = pars,
    ode_control = dust2::dust_ode_control(max_steps = 1e6)
  )
  dust2::dust_system_set_state_initial(sys)
  t <- seq(0, 1000, by = 1)
  out <- dust2::dust_system_simulate(sys, t)
  nt <- neaten_state(system_used = sys, dust_state = out, times_simulated = t)

  final_R <- sum(nt$R1[nrow(nt)], nt$R2[nrow(nt)], nt$R3[nrow(nt)]) # Assuming 3 classes
  simulated_fs <- final_R / N

  # Tolerance: 5% relative error is acceptable for ODE vs analytical
  expect_equal(simulated_fs, fs_analytical, tolerance = 0.05)
})

test_that("final size matches finalsize::final_size for heterogeneous mixing across assortativities", {
  skip_if_not_installed("odin2")
  skip_if_not_installed("finalsize")

  check_scenario <- function(n_act, eps) {
    N <- 1000 * n_act # Scale N with classes
    class_sizes <- rep(N / n_act, n_act)
    # Generate activity scores: linearly spaced from 1 to 10
    act_scores <- seq(1, 10, length.out = n_act)

    m <- garnett_mixing(
      n_activity = n_act,
      epsilon_assort = eps,
      act_vector = act_scores,
      class_size_vector = class_sizes
    )

    # Use scalar parameter helper
    scalars <- basic_scalar_parms()

    # R0 depends on configuration
    pars <- list(
      n_activity = n_act,
      class_sizes = class_sizes,
      activity_scores = act_scores,
      beta = scalars$beta,
      gamma = scalars$gamma,
      sigma = scalars$sigma,
      m = m,
      swch = matrix(0, n_act, n_act),
      S0 = rep(1 - 1e-4, n_act),
      E0 = rep(1e-4, n_act),
      N = N,
      omega = 0,
      t_end = 2000
    )

    # 1. Compute R0 and NGM
    r0_mix <- compute_r0(pars)
    ngm <- compute_ngm(pars)

    # Normalize NGM so its spectral radius is 1.0, to satisfy finalsize check
    # But pass r0_mix explicitly to final_size, which will multiply it back
    if (abs(r0_mix) < 1e-10) {
      ngm_norm <- ngm
    } else {
      ngm_norm <- ngm / r0_mix
    }

    # 2. Predicted (Analytical)
    fs_analytical <- finalsize::final_size(
      r0 = r0_mix,
      contact_matrix = ngm_norm,
      demography_vector = rep(1, n_act), # Use uniform as NGM handles structure
      susceptibility = matrix(1, nrow = n_act, ncol = 1),
      p_susceptibility = matrix(1, nrow = n_act, ncol = 1)
    )$p_infected

    # 3. Simulated
    sys <- dust2::dust_system_create(
      seirs_act(),
      pars,
      ode_control = dust2::dust_ode_control(max_steps = 1e6)
    )
    dust2::dust_system_set_state_initial(sys)
    t <- seq(0, 1000, by = 1)
    out <- dust2::dust_system_simulate(sys, t)
    nt <- neaten_state(system_used = sys, dust_state = out, times_simulated = t)

    sim_p_inf <- numeric(n_act)
    for (i in 1:n_act) {
      sim_p_inf[i] <- nt[[paste0("R", i)]][nrow(nt)] / class_sizes[i]
    }

    # Comparison
    # We verify each group's final size is close
    expect_equal(sim_p_inf, fs_analytical,
      tolerance = 0.05,
      label = paste("Mismatch:: n_act =", n_act, ", eps =", eps)
    )
  }

  # Test across a range of classes and assortativities
  # We use a mix of class numbers to cover the requested 2-15 range efficiently
  # without running excessive simulations.
  class_configs <- c(2, 5, 10, 15)
  eps_configs <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1)

  for (n in class_configs) {
    for (eps in eps_configs) {
      check_scenario(n, eps)
    }
  }
})
