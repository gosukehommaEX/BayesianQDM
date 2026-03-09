#' Plot Method for pbayesdecisionprob2cont Objects
#'
#' Displays operating characteristics for two-continuous-endpoint results
#' returned by \code{\link{pbayesdecisionprob2cont}}.
#'
#' When the input scenarios form a regular grid over
#' \code{(mu_t1, mu_t2)} (i.e., every combination of the unique values of
#' \code{mu_t1} and \code{mu_t2} is present), the function produces a
#' \strong{filled tile plot}: each panel (Go, Gray, NoGo) is coloured by its
#' own probability on a continuous gradient (white to the panel colour), so
#' intensity directly reflects the probability magnitude.
#' Otherwise the function falls back to a \strong{scatter plot} in which point
#' colour encodes the decision probability on a continuous scale.
#'
#' When \code{which = "all"}, the three panels are arranged side-by-side using
#' \code{gridExtra::grid.arrange}, so each panel retains its own independent
#' colour scale.  This requires the \pkg{gridExtra} package.
#'
#' For \code{design = 'controlled'} or \code{design = 'external'}, both axes
#' are expressed as treatment-minus-control differences:
#' \eqn{\theta_1 = \mu_{t1} - \bar{\mu}_{c1}} and
#' \eqn{\theta_2 = \mu_{t2} - \bar{\mu}_{c2}},
#' where \eqn{\bar{\mu}_{c1}} and \eqn{\bar{\mu}_{c2}} are the means of the
#' supplied \code{mu_c1} and \code{mu_c2} vectors.
#' For \code{design = 'uncontrolled'}, the axes represent \eqn{\mu_{t1}} and
#' \eqn{\mu_{t2}} directly.
#'
#' Vertical and horizontal reference lines are drawn at the decision thresholds:
#' \itemize{
#'   \item When \code{prob = 'posterior'}: vertical lines at \eqn{\theta_{TV1}}
#'         and \eqn{\theta_{MAV1}} (x-axis) and horizontal lines at
#'         \eqn{\theta_{TV2}} and \eqn{\theta_{MAV2}} (y-axis).
#'   \item When \code{prob = 'predictive'}: a single vertical line at
#'         \eqn{\theta_{NULL1}} and a single horizontal line at
#'         \eqn{\theta_{NULL2}}.
#' }
#'
#' @param x An object of class \code{pbayesdecisionprob2cont}.
#' @param which A character string specifying which decision probability to
#'        plot.  Must be one of \code{"Go"}, \code{"Gray"}, \code{"NoGo"}, or
#'        \code{"all"}.  When \code{"all"}, all three panels are arranged
#'        side-by-side via \code{gridExtra::grid.arrange}.
#'        Default is \code{"Go"}.
#' @param title A character string for the plot title.  Defaults to
#'        \code{NULL} (no title displayed).
#' @param xlab A character string or expression for the x-axis label.
#'        Defaults to \code{NULL}, which auto-generates a label based on
#'        \code{design}.
#' @param ylab A character string or expression for the y-axis label.
#'        Defaults to \code{NULL}, which auto-generates a label based on
#'        \code{design}.
#' @param col_go A character string specifying the high-end fill colour for the
#'        Go probability gradient.  Default is \code{"#004C97"}.
#' @param col_nogo A character string specifying the high-end fill colour for
#'        the NoGo probability gradient.  Default is \code{"#F0B323"}.
#' @param col_gray A character string specifying the high-end fill colour for
#'        the Gray probability gradient.  Default is \code{"gray60"}.
#' @param base_size A positive numeric scalar specifying the base font size
#'        (in points) passed to \code{theme_bw()}.  Default is \code{28}.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns a \code{ggplot} object (single panel) or a
#'         \code{gtable} object (\code{which = "all"}).
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_vline geom_hline annotate
#'   scale_fill_gradient scale_x_continuous scale_y_continuous coord_cartesian
#'   labs theme_bw theme element_text element_line element_blank margin unit
#'   geom_point scale_color_gradient
#' @importFrom gridExtra grid.arrange
#' @export
plot.pbayesdecisionprob2cont <- function(x,
                                         which     = "Go",
                                         title     = NULL,
                                         xlab      = NULL,
                                         ylab      = NULL,
                                         col_go    = "#004C97",
                                         col_nogo  = "#F0B323",
                                         col_gray  = "gray60",
                                         base_size = 28,
                                         ...) {

  # --- Input validation ---
  if (!which %in% c("Go", "Gray", "NoGo", "all")) {
    stop("'which' must be one of \"Go\", \"Gray\", \"NoGo\", or \"all\"")
  }
  if (nrow(x) < 2L) {
    stop("'x' must contain at least 2 scenarios (rows) to produce a plot")
  }

  # --- Extract attributes ---
  design      <- attr(x, "design")
  prob        <- attr(x, "prob")
  theta_TV1   <- attr(x, "theta_TV1")
  theta_MAV1  <- attr(x, "theta_MAV1")
  theta_TV2   <- attr(x, "theta_TV2")
  theta_MAV2  <- attr(x, "theta_MAV2")
  theta_NULL1 <- attr(x, "theta_NULL1")
  theta_NULL2 <- attr(x, "theta_NULL2")
  gamma_go    <- attr(x, "gamma_go")
  gamma_nogo  <- attr(x, "gamma_nogo")

  # --- Build x/y axis variables ---
  if (design %in% c("controlled", "external")) {
    ax1     <- x$mu_t1 - mean(x$mu_c1)
    ax2     <- x$mu_t2 - mean(x$mu_c2)
    x_label <- if (is.null(xlab)) expression(theta[1] == mu[t1] - mu[c1]) else xlab
    y_label <- if (is.null(ylab)) expression(theta[2] == mu[t2] - mu[c2]) else ylab
  } else {
    ax1     <- x$mu_t1
    ax2     <- x$mu_t2
    x_label <- if (is.null(xlab)) expression(mu[t1]) else xlab
    y_label <- if (is.null(ylab)) expression(mu[t2]) else ylab
  }

  # --- Threshold values used directly as axis coordinates (no transformation) ---
  if (prob == "posterior") {
    vx_TV  <- theta_TV1;  vx_MAV <- theta_MAV1
    vy_TV  <- theta_TV2;  vy_MAV <- theta_MAV2
  } else {
    vx_NULL <- theta_NULL1
    vy_NULL <- theta_NULL2
  }

  # --- Detect grid layout ---
  u1      <- sort(unique(ax1))
  u2      <- sort(unique(ax2))
  is_grid <- (length(u1) * length(u2) == nrow(x))

  # --- Helper: axis breaks ---
  axis_breaks <- function(vals) {
    b <- pretty(vals, n = 6L)
    b[b >= min(vals) & b <= max(vals)]
  }

  # --- Helper: add threshold reference lines and labels ---
  add_thresholds <- function(p, x_rng, y_rng, bs) {
    off_x <- diff(x_rng) * 0.02
    off_y <- diff(y_rng) * 0.02
    if (prob == "posterior") {
      p <- p +
        ggplot2::geom_vline(xintercept = vx_MAV, color = "gray30",
                            linetype = "dotted", linewidth = 1.0) +
        ggplot2::geom_vline(xintercept = vx_TV,  color = "gray30",
                            linetype = "dotted", linewidth = 1.0) +
        ggplot2::geom_hline(yintercept = vy_MAV, color = "gray30",
                            linetype = "dotted", linewidth = 1.0) +
        ggplot2::geom_hline(yintercept = vy_TV,  color = "gray30",
                            linetype = "dotted", linewidth = 1.0) +
        ggplot2::annotate("text",
                          x = vx_MAV - off_x, y = y_rng[2],
                          label = "theta[MAV1]", parse = TRUE,
                          color = "gray30", hjust = 1, vjust = 1,
                          size = bs / 3) +
        ggplot2::annotate("text",
                          x = vx_TV  + off_x, y = y_rng[2],
                          label = "theta[TV1]",  parse = TRUE,
                          color = "gray30", hjust = 0, vjust = 1,
                          size = bs / 3) +
        ggplot2::annotate("text",
                          x = x_rng[1], y = vy_MAV - off_y,
                          label = "theta[MAV2]", parse = TRUE,
                          color = "gray30", hjust = 0, vjust = 1,
                          size = bs / 3) +
        ggplot2::annotate("text",
                          x = x_rng[1], y = vy_TV  + off_y,
                          label = "theta[TV2]",  parse = TRUE,
                          color = "gray30", hjust = 0, vjust = 0,
                          size = bs / 3)
    } else {
      p <- p +
        ggplot2::geom_vline(xintercept = vx_NULL, color = "gray30",
                            linetype = "dotted", linewidth = 1.0) +
        ggplot2::geom_hline(yintercept = vy_NULL, color = "gray30",
                            linetype = "dotted", linewidth = 1.0) +
        ggplot2::annotate("text",
                          x = vx_NULL + off_x, y = y_rng[2],
                          label = "theta[NULL1]", parse = TRUE,
                          color = "gray30", hjust = 0, vjust = 1,
                          size = bs / 3) +
        ggplot2::annotate("text",
                          x = x_rng[1], y = vy_NULL + off_y,
                          label = "theta[NULL2]", parse = TRUE,
                          color = "gray30", hjust = 0, vjust = 0,
                          size = bs / 3)
    }
    p
  }

  # --- Helper: common theme ---
  common_theme <- function(bs) {
    ggplot2::theme_bw(base_size = bs) +
      ggplot2::theme(
        legend.position   = "bottom",
        legend.margin     = ggplot2::margin(t = -10, r = 0, b = 0, l = 0),
        legend.box.margin = ggplot2::margin(t = -10, r = 0, b = 0, l = 0),
        legend.text       = ggplot2::element_text(size = bs * 0.54),
        legend.title      = ggplot2::element_text(size = bs * 0.54),
        legend.key.width  = ggplot2::unit(1.5, "cm"),
        panel.grid.minor  = ggplot2::element_blank(),
        panel.border      = ggplot2::element_blank(),
        axis.line         = ggplot2::element_line(color = "black",
                                                  linewidth = 0.8),
        strip.text        = ggplot2::element_text(size = bs * 0.64,
                                                  face = "bold")
      )
  }

  # ---------------------------------------------------------------------------
  # TILE MODE
  # ---------------------------------------------------------------------------
  if (is_grid) {

    x_rng <- range(ax1)
    y_rng <- range(ax2)

    # Helper: build a single gradient tile panel
    make_tile_panel <- function(prob_col, high_col, panel_title, bs) {
      df <- data.frame(ax1      = ax1,
                       ax2      = ax2,
                       prob_val = x[[prob_col]])

      p <- ggplot2::ggplot(df, ggplot2::aes(x = ax1, y = ax2)) +
        ggplot2::geom_tile(ggplot2::aes(fill = prob_val)) +
        ggplot2::scale_fill_gradient(
          name   = paste0("P(", prob_col, ")"),
          low    = "white",
          high   = high_col,
          limits = c(0, 1)
        ) +
        ggplot2::scale_x_continuous(breaks = axis_breaks(ax1),
                                    expand = c(0, 0)) +
        ggplot2::scale_y_continuous(breaks = axis_breaks(ax2),
                                    expand = c(0, 0)) +
        ggplot2::coord_cartesian(xlim = x_rng, ylim = y_rng) +
        ggplot2::labs(title = panel_title, x = x_label, y = y_label) +
        common_theme(bs)

      p <- add_thresholds(p, x_rng, y_rng, bs)
      p
    }

    if (which == "all") {
      if (!requireNamespace("gridExtra", quietly = TRUE)) {
        stop("Package 'gridExtra' is required for which = \"all\". ",
             "Please install it with install.packages(\"gridExtra\").")
      }
      p_go   <- make_tile_panel("Go",   col_go,   "Go",   base_size)
      p_gray <- make_tile_panel("Gray", col_gray, "Gray", base_size)
      p_nogo <- make_tile_panel("NoGo", col_nogo, "NoGo", base_size)

      out <- gridExtra::grid.arrange(
        p_go, p_gray, p_nogo,
        nrow = 1L,
        top  = if (!is.null(title)) title else ""
      )
      return(invisible(out))

    } else {
      high_col <- switch(which, Go = col_go, Gray = col_gray, NoGo = col_nogo)
      p <- make_tile_panel(which, high_col, title, base_size)
    }

    # ---------------------------------------------------------------------------
    # SCATTER MODE (non-grid)
    # ---------------------------------------------------------------------------
  } else {

    prob_col <- if (which == "all") {
      message("Non-grid input: 'which = \"all\"' shows Go probability only.")
      "Go"
    } else {
      which
    }

    high_col  <- switch(prob_col, Go = col_go, NoGo = col_nogo, Gray = col_gray)
    x_rng     <- range(ax1)
    y_rng     <- range(ax2)

    df_scatter <- data.frame(ax1      = ax1,
                             ax2      = ax2,
                             prob_val = x[[prob_col]])

    p <- ggplot2::ggplot(df_scatter,
                         ggplot2::aes(x = ax1, y = ax2, color = prob_val)) +
      ggplot2::geom_point(size = 4) +
      ggplot2::scale_color_gradient(
        name   = paste0("P(", prob_col, ")"),
        low    = "white",
        high   = high_col,
        limits = c(0, 1)
      ) +
      ggplot2::scale_x_continuous(breaks = axis_breaks(ax1),
                                  expand = c(0.05, 0)) +
      ggplot2::scale_y_continuous(breaks = axis_breaks(ax2),
                                  expand = c(0.05, 0)) +
      ggplot2::labs(title = title, x = x_label, y = y_label) +
      common_theme(base_size)

    p <- add_thresholds(p, x_rng, y_rng, base_size)
  }

  print(p)
  invisible(p)
}
