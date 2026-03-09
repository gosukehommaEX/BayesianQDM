#' Plot Method for pbayesdecisionprob1bin Objects
#'
#' Displays an operating characteristics curve of Go/NoGo/Gray decision
#' probabilities against the true treatment effect for binary endpoint results
#' returned by \code{\link{pbayesdecisionprob1bin}}.
#'
#' For \code{design = 'controlled'} or \code{design = 'external'}, the
#' x-axis represents the treatment-minus-control difference
#' \eqn{\theta = \pi_t - \bar{\pi}_c}, where \eqn{\bar{\pi}_c} is the mean
#' of the supplied \code{pi_c} values. For \code{design = 'uncontrolled'},
#' the x-axis represents \eqn{\pi_t} directly.
#'
#' Vertical reference lines are drawn at the decision thresholds:
#' \itemize{
#'   \item When \code{prob = 'posterior'}: lines at \eqn{\theta_{TV}} and
#'         \eqn{\theta_{MAV}} (converted to the x-axis scale).
#'   \item When \code{prob = 'predictive'}: a single line at \eqn{\theta_{NULL}}.
#' }
#'
#' @param x An object of class \code{pbayesdecisionprob1bin}.
#' @param title A character string for the plot title. Defaults to
#'        \code{NULL} (no title displayed).
#' @param xlab A character string or expression for the x-axis label.
#'        Defaults to \code{NULL}, which auto-generates a label based on
#'        \code{design}.
#' @param col_go A character string specifying the colour for the Go curve.
#'        Default is \code{"#004C97"}.
#' @param col_nogo A character string specifying the colour for the NoGo curve.
#'        Default is \code{"#F0B323"}.
#' @param col_gray A character string specifying the colour for the Gray curve.
#'        Default is \code{"gray60"}.
#' @param base_size A positive numeric scalar specifying the base font size
#'        (in points) passed to \code{theme_bw()}. Default is \code{28}.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns a \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_vline annotate
#'   scale_color_manual scale_linetype_manual scale_x_continuous
#'   scale_y_continuous labs theme_bw theme element_text element_line
#'   element_blank margin unit
#' @export
plot.pbayesdecisionprob1bin <- function(x,
                                        title    = NULL,
                                        xlab     = NULL,
                                        col_go   = "#004C97",
                                        col_nogo = "#F0B323",
                                        col_gray = "gray60",
                                        base_size = 28,
                                        ...) {

  # --- Validate minimum number of scenarios ---
  if (nrow(x) < 2L) {
    stop("'x' must contain at least 2 scenarios (rows) to produce a line plot")
  }

  # --- Extract attributes ---
  design        <- attr(x, "design")
  prob          <- attr(x, "prob")
  theta_TV      <- attr(x, "theta_TV")
  theta_MAV     <- attr(x, "theta_MAV")
  theta_NULL    <- attr(x, "theta_NULL")

  # --- Build x-axis variable ---
  if (design %in% c("controlled", "external")) {
    # Use the mean of pi_c as the reference for the difference scale
    pi_c_mean <- mean(x$pi_c)
    theta      <- x$pi_t - pi_c_mean
    x_label    <- if (is.null(xlab)) {
      expression(theta == pi[t] - pi[c])
    } else {
      xlab
    }
  } else {
    # Uncontrolled: x-axis is pi_t directly
    theta   <- x$pi_t
    x_label <- if (is.null(xlab)) expression(pi[t]) else xlab
  }

  # --- Reshape to long format ---
  oc_long <- data.frame(
    theta    = rep(theta, 3L),
    prob     = c(x$Go, x$NoGo, x$Gray),
    Decision = rep(c("Go", "NoGo", "Gray"), each = nrow(x))
  )
  oc_long$Decision <- factor(oc_long$Decision, levels = c("Go", "NoGo", "Gray"))

  # --- Compute axis limits with a small margin ---
  x_min  <- min(theta)
  x_max  <- max(theta)
  x_pad  <- (x_max - x_min) * 0.02
  x_rng  <- c(x_min - x_pad, x_max + x_pad)

  x_breaks <- pretty(theta, n = 6L)
  x_breaks <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]

  # --- Compute threshold positions on the x-axis ---
  # theta_TV, theta_MAV, and theta_NULL are already on the difference scale
  # (pi_t - pi_c), so they map directly to the x-axis without adjustment.
  if (prob == "posterior") {
    vline_TV  <- theta_TV
    vline_MAV <- theta_MAV
  } else {
    vline_NULL <- theta_NULL
  }

  # --- Annotation y position: slightly below the top of the plot ---
  annot_y <- 0.95

  # --- Compute horizontal offset for annotation text ---
  # Use 2% of the x range so the label clears the vline regardless of scale
  x_offset <- (x_max - x_min) * 0.02

  # --- Base plot ---
  p <- ggplot2::ggplot(oc_long,
                       ggplot2::aes(x = theta, y = prob,
                                    color    = Decision,
                                    linetype = Decision)) +
    ggplot2::geom_line(linewidth = 1.5)

  # --- Add threshold vertical lines and labels ---
  if (prob == "posterior") {
    p <- p +
      ggplot2::geom_vline(xintercept = vline_MAV, color = "gray30",
                          linetype = "dotted", linewidth = 1.2) +
      ggplot2::geom_vline(xintercept = vline_TV,  color = "gray30",
                          linetype = "dotted", linewidth = 1.2) +
      ggplot2::annotate("text",
                        x     = vline_MAV - x_offset,
                        y     = annot_y,
                        label = "theta[MAV]",
                        parse = TRUE,
                        color = "gray30",
                        hjust = 1,
                        size  = base_size / 3) +
      ggplot2::annotate("text",
                        x     = vline_TV + x_offset,
                        y     = annot_y,
                        label = "theta[TV]",
                        parse = TRUE,
                        color = "gray30",
                        hjust = 0,
                        size  = base_size / 3)
  } else {
    p <- p +
      ggplot2::geom_vline(xintercept = vline_NULL, color = "gray30",
                          linetype = "dotted", linewidth = 1.2) +
      ggplot2::annotate("text",
                        x     = vline_NULL + x_offset,
                        y     = annot_y,
                        label = "theta[NULL]",
                        parse = TRUE,
                        color = "gray30",
                        hjust = 0,
                        size  = base_size / 3)
  }

  # --- Scales, labels, and theme ---
  p <- p +
    ggplot2::scale_color_manual(
      values = c("Go" = col_go, "NoGo" = col_nogo, "Gray" = col_gray)
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Go" = "solid", "NoGo" = "twodash", "Gray" = "dashed")
    ) +
    ggplot2::scale_x_continuous(
      limits = x_rng,
      breaks = x_breaks,
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      title    = title,
      x        = x_label,
      y        = "Probability",
      color    = "Decision",
      linetype = "Decision"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position   = "bottom",
      legend.margin     = ggplot2::margin(t = -10, r = 0, b = 0, l = 0),
      legend.box.margin = ggplot2::margin(t = -10, r = 0, b = 0, l = 0),
      legend.spacing.y  = ggplot2::unit(0, "cm"),
      legend.text       = ggplot2::element_text(size = base_size * 0.54),
      legend.title      = ggplot2::element_blank(),
      legend.key.width  = ggplot2::unit(2, "cm"),
      panel.grid.minor  = ggplot2::element_blank(),
      panel.border      = ggplot2::element_blank(),
      axis.line         = ggplot2::element_line(color = "black", linewidth = 0.8)
    )

  print(p)
  invisible(p)
}
