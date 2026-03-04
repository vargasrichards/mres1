library(flodia)
library(ggplot2)
library(patchwork)
library(MixSwitchEpi)

seir_class <- function(x = 1, y = 1,
                       class_num,
                       col = list()) {
  r <- 0.1
  xgap <- 0.2

  S <- node(
    x = x, y = y, r = r,
    label = bquote(S[.(class_num)]),
    node_col = col$S
  )

  E <- node(
    x = S$x1 + xgap, y = y, r = r,
    label = bquote(E[.(class_num)]),
    node_col = col$E
  )

  I <- node(
    x = E$x1 + xgap, y = y, r = r,
    label = bquote(I[.(class_num)]),
    node_col = col$I
  )

  R <- node(
    x = I$x1 + xgap, y = y, r = r,
    label = bquote(R[.(class_num)]),
    node_col = col$R
  )

  ## Disease progression
  flowx(from = S, to = E, label = bquote(lambda[.(class_num)]))
  flowx(from = E, to = I, label = expression(sigma))
  flowx(from = I, to = R, label = expression(gamma))

  list(
    x0 = S$x0, x1 = R$x1,
    y0 = S$y0, y1 = S$y1
  )
}


seir_activity <- function(K = 3) {
  ## vertical spacing between activity classes
  ygap <- 0.5

  col <- list(
    S = light_palette("gnbu"),
    E = light_palette("rdor"),
    I = light_palette("ylgn"),
    R = light_palette("bupu")
  )

  classes <- vector("list", K)

  ## draw and group each activity class
  for (i in seq_len(K)) {
    classes[[i]] <- group(
      f = seir_class,
      args = list(
        x = 1,
        y = 1 - (i - 1) * ygap,
        class_num = i,
        col = col
      ),
      group_col = grey(0.95),
      oma = c(0.05, 0.05, 0.05, 0.05)
    )
  }

  for (i in seq_len(K - 1)) {
    for (j in (i + 1):K) {
      ij <- as.integer(paste(i, j, sep = ""))
      ji <- as.integer(paste(j, i, sep = ""))

      ## jump distance (how many classes we are jumping over)
      jump <- abs(i - j)
      ## offset flows to the sides (left for i->j, right for j->i)
      ## vary offset by jump distance to avoid overlapping lines
      ## use smaller offsets to keep labels within bounds
      offset_ij <- -0.05 * jump
      offset_ji <- 1.0 + 0.05 * jump

      ## i -> j (left side)
      flowy(
        from = classes[[i]],
        to = classes[[j]],
        label = bquote(W[.(ji)]),
        pos = offset_ij
      )

      ## j -> i (right side)
      flowy(
        from = classes[[j]],
        to = classes[[i]],
        label = bquote(W[.(ij)]),
        pos = offset_ji
      )
    }
  }

  list(
    x0 = min(vapply(classes, `[[`, numeric(1), "x0")),
    x1 = max(vapply(classes, `[[`, numeric(1), "x1")),
    y0 = classes[[K]]$y0,
    y1 = classes[[1]]$y1
  )
}

#' Make Figure 1: flodia SEIR activity diagram + Rt/epidemic metrics
#'
#' Draws a 2-class flodia SEIR diagram (LHS) and a homogeneous SEIR
#' epidemic (RHS) with an Rt panel aligned below the epidemic plot. The
#' combined figure is saved as `figure1.pdf` in the project root.
make_figure1 <- function(
  K = 2,
  width = 10,
  height = 5,
  target_r0 = 2,
  x_truncation = 300
) {
  # Build homogeneous (single-population) calibrated parameters
  base_pars <- MixSwitchEpi::make_parms(
    n_activity = 1,
    epsilon = 0,
    switch_rate = 0,
    initial_condition_rule = "uniform",
    infected_fraction = 1e-3,
    activity_scheme = "null",
    mixing_model = "garnett_mixing",
    switching_model = "simple_switching",
    pathogen_string = "sars-cov-2",
    target_r0 = target_r0
  )

  # simulate the homogeneous epidemic using existing helper
  epi_out <- MixSwitchEpi::run_system_wparms(
    parms = base_pars,
    model = MixSwitchEpi:::seirs_act()
  )

  # Compute metrics and add Rt timeseries (adds R0 and Rt columns)
  metrics <- MixSwitchEpi::compute_all(epi_out, classed_detail = FALSE, params = base_pars)
  epi_with_rt <- MixSwitchEpi::compute_rt_ts(epi_out, base_pars, add_to_output = TRUE)

  p_epi <- MixSwitchEpi::plot_epidemic_annotated(
    model_output = epi_with_rt,
    metrics = metrics,
    compartments = c("S", "E", "I", "R"),
    show_total = TRUE,
    annotate_peak = TRUE,
    annotate_final = TRUE,
    title = sprintf("Homogeneous SEIR (R0 = %.2f)", target_r0),
    x_truncation = x_truncation,
    axis_text_size = 13,
    axis_title_size = 13
  )
  p_rt <- MixSwitchEpi::plot_rt_timeseries(epi_with_rt,
    highlight_threshold = TRUE,
    x_truncation = x_truncation,
    axis_text_size = 13,
    axis_title_size = 13,
    metrics = metrics
  )

  rt_epi_plot <- p_epi / p_rt + patchwork::plot_annotation(tag_levels = "A")
  ggplot2::ggsave(plot = rt_epi_plot, filename = "homogeneous_epidemic_w_rt.pdf")
  # Compose final figure and save using grid viewports so flodia (grid)
  # output can be drawn directly alongside ggplot objects.
  outfile_pdf <- "figure1.pdf"
  outfile_png <- "figure1.png"
  # make RHS text smaller and remove lower panel title
  p_epi <- p_epi + ggplot2::theme(
    plot.title = ggplot2::element_text(size = 12),
    axis.text = ggplot2::element_text(size = 10),
    axis.title = ggplot2::element_text(size = 10)
  )
  p_rt <- p_rt + ggplot2::labs(title = NULL) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 10)
    )
  grid::grid.newpage()
  # Make left-hand (flodia) panel slightly larger than before
  layout <- grid::grid.layout(
    nrow = 1, ncol = 2,
    widths = grid::unit.c(grid::unit(0.40, "npc"), grid::unit(0.60, "npc"))
  )
  grid::pushViewport(grid::viewport(layout = layout))

  # Left: draw flodia directly
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  tryCatch(
    {
      flodia::flodia(function() seir_activity(K), oma = rep(0.2, 4))
    },
    error = function(e) {
      grid::grid.text("(flodia diagram failed to render)", gp = grid::gpar(col = "red"))
    }
  )
  grid::popViewport()

  # Right: print the stacked epidemic + Rt patchwork
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  tryCatch(
    {
      print(rt_epi_plot, newpage = FALSE)
    },
    error = function(e) {
      grid::grid.text("(epidemic plot failed to render)", gp = grid::gpar(col = "red"))
    }
  )

  grDevices::dev.off()

  # Also save as PNG (high-res)
  png_dpi <- 300
  grDevices::png(outfile_png,
    width = width * png_dpi, height = height * png_dpi,
    res = png_dpi, type = "cairo"
  )
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = layout))
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  tryCatch(
    {
      flodia::flodia(function() seir_activity(K), oma = rep(0.2, 4))
    },
    error = function(e) {
      grid::grid.text("(flodia diagram failed to render)",
        gp = grid::gpar(col = "red")
      )
    }
  )
  grid::popViewport()
  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  tryCatch(
    {
      print(rt_epi_plot, newpage = FALSE)
    },
    error = function(e) {
      grid::grid.text("(epidemic plot failed to render)",
        gp = grid::gpar(col = "red")
      )
    }
  )
  grid::popViewport()
  grDevices::dev.off()

  # Return filenames and the RHS patchwork object
  # (as the combined plot can't be a single ggplot)
  invisible(list(pdf = outfile_pdf, png = outfile_png, rhs_plot = rt_epi_plot))
  rt_epi_plot
}
