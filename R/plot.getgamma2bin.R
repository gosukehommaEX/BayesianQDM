#' Plot Method for getgamma2bin Objects
#'
#' Displays calibration curves of marginal Go and NoGo probabilities against
#' the threshold grid \eqn{\gamma} for results returned by
#' \code{\link{getgamma2bin}}.
#'
#' The x-axis represents candidate threshold values \eqn{\gamma \in (0, 1)}.
#' The y-axis represents the marginal probability:
#' \itemize{
#'   \item \strong{Go curve}: \eqn{\Pr(g_{\mathrm{Go}} \ge \gamma)} evaluated
#'         under the Go-calibration scenario (typically the Null scenario).
#'   \item \strong{NoGo curve}: \eqn{\Pr(g_{\mathrm{NoGo}} \ge \gamma)}
#'         evaluated under the NoGo-calibration scenario (typically the
#'         Alternative scenario).
#' }
#'
#' Horizontal reference lines are drawn at \code{target_go} and
#' \code{target_nogo}.  Filled circles (\code{geom_point}) mark the optimal
#' thresholds \eqn{(\gamma_{\mathrm{go}},\, \Pr(\mathrm{Go})_{\mathrm{opt}})}
#' and \eqn{(\gamma_{\mathrm{nogo}},\, \Pr(\mathrm{NoGo})_{\mathrm{opt}})},
#' with their values shown in the legend.
#' If either optimal threshold is \code{NA}, the corresponding point is omitted.
#'
#' @param x An object of class \code{getgamma2bin}.
#' @param title A character string for the plot title.  Defaults to
#'        \code{NULL} (no title displayed).
#' @param col_go A character string specifying the colour for the Go curve.
#'        Default is \code{"#658D1B"}.
#' @param col_nogo A character string specifying the colour for the NoGo curve.
#'        Default is \code{"#D91E49"}.
#' @param base_size A positive numeric scalar specifying the base font size
#'        (in points) passed to \code{theme_bw()}.  Default is \code{28}.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns a \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline annotate
#'   scale_color_manual scale_linetype_manual scale_fill_manual
#'   scale_x_continuous scale_y_continuous guides guide_legend
#'   labs theme_bw theme element_text element_line
#'   element_blank margin unit
#' @importFrom stats setNames
#' @export
plot.getgamma2bin <- function(x,
                              title     = NULL,
                              col_go    = "#658D1B",
                              col_nogo  = "#D91E49",
                              base_size = 28,
                              ...) {

  # --- Input validation ---
  if (!inherits(x, "getgamma2bin")) {
    stop("'x' must be an object of class 'getgamma2bin'")
  }

  # --- Extract target values stored in x ---
  target_go   <- x$target_go
  target_nogo <- x$target_nogo

  # --- Extract grid results ---
  df         <- x$grid_results
  gamma_go   <- x$gamma_go
  gamma_nogo <- x$gamma_nogo

  # --- Reshape to long format ---
  plot_long <- data.frame(
    gamma    = rep(df$gamma_grid, 2L),
    prob     = c(df$PrGo_grid, df$PrNoGo_grid),
    Scenario = rep(c("Go", "NoGo"), each = nrow(df))
  )
  plot_long$Scenario <- factor(plot_long$Scenario, levels = c("Go", "NoGo"))

  # --- Axis limits ---
  x_min <- min(df$gamma_grid)
  x_max <- max(df$gamma_grid)
  x_pad <- (x_max - x_min) * 0.02
  x_rng <- c(x_min - x_pad, x_max + x_pad)

  x_breaks <- pretty(df$gamma_grid, n = 6L)
  x_breaks <- x_breaks[x_breaks >= x_min & x_breaks <= x_max]

  # --- Build data frame for optimal threshold points ---
  # Points are placed at (gamma_go, PrGo_opt) and (gamma_nogo, PrNoGo_opt).
  # The 'fill' aesthetic drives a separate legend with expression labels.
  go_val   <- if (is.na(gamma_go))   "NA" else sprintf("%.2f", gamma_go)
  nogo_val <- if (is.na(gamma_nogo)) "NA" else sprintf("%.2f", gamma_nogo)

  pt_labels <- c(
    paste0("gamma[go]==", go_val),
    paste0("gamma[nogo]==", nogo_val)
  )

  pt_df <- data.frame(
    gamma   = c(gamma_go,       gamma_nogo),
    prob    = c(x$PrGo_opt,     x$PrNoGo_opt),
    pt_col  = c(col_go,         col_nogo),
    pt_grp  = factor(pt_labels, levels = pt_labels)
  )
  # Remove rows where the optimal threshold is NA
  pt_df <- pt_df[!is.na(pt_df$gamma), ]

  # Named colour vector for scale_fill_manual
  pt_col_vals <- stats::setNames(pt_df$pt_col, as.character(pt_df$pt_grp))

  # --- Base plot ---
  p <- ggplot2::ggplot(plot_long,
                       ggplot2::aes(x = gamma, y = prob,
                                    color    = Scenario,
                                    linetype = Scenario)) +
    ggplot2::geom_line(linewidth = 1.5)

  # --- Horizontal reference lines ---
  p <- p +
    ggplot2::geom_hline(yintercept = target_go,
                        color = col_go, linetype = "dotted",
                        linewidth = 1.5) +
    ggplot2::geom_hline(yintercept = target_nogo,
                        color = col_nogo, linetype = "dotted",
                        linewidth = 1.5)

  # --- Optimal threshold points ---
  if (nrow(pt_df) > 0L) {
    p <- p +
      ggplot2::geom_point(
        data        = pt_df,
        mapping     = ggplot2::aes(x = gamma, y = prob, fill = pt_grp),
        color       = pt_df$pt_col,
        shape       = 21,
        size        = 5,
        stroke      = 1.5,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = pt_col_vals,
        labels = lapply(names(pt_col_vals), function(lb) parse(text = lb)[[1L]])
      )
  }

  # --- Scales, labels, and theme ---
  p <- p +
    ggplot2::scale_color_manual(
      values = c("Go" = col_go, "NoGo" = col_nogo)
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Go" = "solid", "NoGo" = "twodash")
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
      x        = expression(gamma),
      y        = "Probability",
      color    = "Scenario",
      linetype = "Scenario",
      fill     = "Optimal threshold"
    ) +
    ggplot2::guides(
      color    = ggplot2::guide_legend(order = 1),
      linetype = ggplot2::guide_legend(order = 1),
      fill     = ggplot2::guide_legend(order = 2,
                                       override.aes = list(size = 5))
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
