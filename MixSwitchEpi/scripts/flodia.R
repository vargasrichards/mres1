# flodia for mixswitch
library("flodia")
# note that i didnt write this code, L Whittles did.


sedar <- function(S_x = 1, S_y = 1,
                  label = list(S = "S_1", E = "E", D = "D", A = "A", R = "R"),
                  RS_pos = 0.5) {
  # define the radius of the nodes
  r <- 0.1
  xgap <- 0.3
  ygap <- 0.2

  col <- list(
    S = light_palette("gnbu"),
    E = light_palette("rdor"),
    I = light_palette("ylgn"),
    R = light_palette("bupu")
  )


  # draw S at (S_x, S_y)
  S <- node(x = S_x, y = S_y, r = r, label = label$S, node_col = col$S)
  # draw E xgap units to the right of
  E <- node(x = S$x1 + xgap, y = S$y, r = r, label = label$E, node_col = col$E)
  # draw D 1.5 * xgap units to the right, and ygap units above E
  D <- node(
    x = E$x1 + xgap * 1.8, y = S$y + ygap, r = r, label = label$D,
    node_col = col$I
  )
  # draw A at the same x-coordinate as d, and ygap units below E
  A <- node(x = D$x, y = S$y - ygap, r = r, label = label$A, node_col = col$I)
  # draw R xgap units to the right of A and D
  R <- node(x = D$x1 + xgap, y = S$y, r = r, label = label$R, node_col = col$R)

  # add flows
  # connect S to E
  flowx(from = S, to = E, label = expression(lambda))
  # connect S to A and D
  forkx(
    from = E, to0 = A, to1 = D, label_from = expression(gamma[E]),
    label_to0 = expression(1 - p[D]),
    label_to1 = expression(p[D])
  )
  # connect D to R
  bendy(
    from = D, to = R, pos_to = 0.7,
    label_to = expression(gamma[D])
  )
  # connect A to R
  bendy(
    from = A, to = R, pos_to = 0.3,
    label_to = expression(gamma[A]),
    label_to_gap = -0.05
  )
  # connect R to S - record co-ordinates in object rs as this turn extends the
  # top of the flodia (y1) so must be output
  rs <- turny(
    from = R, mid_y = D$y1 + 0.05, to = S,
    label = expression(gamma[R]), pos_to = RS_pos
  )


  # return a list of the left, right, bottom and top co-ordinates
  # additionally return node S co-ordinates
  list(x0 = S$x0, x1 = R$x1, y0 = A$y0, y1 = rs$y1, S = S)
}


sedar_group <- function(S_x = 1, S_y = 1,
                        label = list(S = "S", E = "E", D = "D", A = "A", R = "R"),
                        RS_pos = 0.5) {
  # group the nodes together
  g <- group(
    f = sedar, args = list(
      S_x = S_x, S_y = S_y, label = label,
      RS_pos = RS_pos
    ),
    group_col = grey(0.95), oma = c(0.1, 0.1, 0.15, 0.1)
  )

  # specify length of inflow and outflow
  l <- 0.3

  # births into S
  births <- flowx(to = g$S, length = l, label = expression(alpha), label_pos = 0.4)
  # deaths out of group
  deaths <- flowx(from = g, length = l, label = expression(mu))

  # return a list of the left, right, bottom and top co-ordinates, and node S
  list(x0 = births$x0, x1 = deaths$x1, y0 = g$y0, y1 = g$y1, S = g$S)
}


vax <- function() {
  labels_U <- c(S = "S", E = "E", D = "D", A = "A", R = "R")
  labels_V <- sprintf("%s*", labels_U)
  names(labels_V) <- labels_U

  U <- sedar_group(label = as.list(labels_U))
  V <- sedar_group(S_y = U$y0 - 0.8, label = as.list(labels_V), RS_pos = 0.8)

  # add vaccination flow - aligning label between U and V
  flowy(
    from = U$S, to = V$S, pos = 0.5, label = expression(nu),
    label_y = calc_pos(U$y0, V$y1)
  )
  # add waning flow
  flowy(from = V, to = U, label = expression(omega))

  # return a list of the left, right, bottom and top co-ordinates
  list(x0 = U$x0, x1 = U$x1, y0 = V$y0, y1 = U$y1)
}

flodia::flodia_pdf(f = vax, filepath = "flodia.pdf")
#
