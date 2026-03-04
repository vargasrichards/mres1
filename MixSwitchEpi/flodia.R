library(flodia)

seir_class <- function(x = 1, y = 1,
                       class_num,
                       col = list()) {
  r <- 0.1
  xgap <- 0.3

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

## ============================================================
## Full SEIR model with all-to-all activity switching
## ============================================================

seir_activity <- function(K) {
  ## vertical spacing between activity classes
  ygap <- 0.7

  ## colour scheme (matching vignette examples)
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

  ## ------------------------------------------------------------
  ## All-to-all activity switching (group-level)
  ## ------------------------------------------------------------

  for (i in seq_len(K - 1)) {
    for (j in (i + 1):K) {
      ij <- as.integer(paste(i, j, sep = ""))
      ji <- as.integer(paste(j, i, sep = ""))

      ## i -> j
      flowy(
        from = classes[[i]],
        to = classes[[j]],
        label = bquote(W[.(ji)]),
        pos = 0.45
      )

      ## j -> i
      flowy(
        from = classes[[j]],
        to = classes[[i]],
        label = bquote(W[.(ij)]),
        pos = 0.55
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

## ============================================================
## Render diagram
## ============================================================

# flodia_png(
#   f = seir_activity,
#   filepath = "fig",
#   res = 200
# )
