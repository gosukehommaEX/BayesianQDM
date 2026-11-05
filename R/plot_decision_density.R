#' Posterior density of the effect of one endpoint with the Go and NoGo
#' regions
#'
#' @param x An object of class \code{"qdm_decision"} for one endpoint and a
#'   criterion based on posterior probability.
#' @param cols Named colors of the decisions (\code{decision_colors()}).
#' @param style The style of the plot (\code{plot_style()}).
#'
#' @return A ggplot object. The density is the central difference of the
#'   distribution function of the effect (\code{decision_cdf()}) on a grid
#'   of about 400 points that contains the thresholds.
#' @noRd
plot_decision_density <- function(x, cols, style) {
  design <- x$design
  crit <- design$criterion
  cdf <- decision_cdf(design, x$posterior_t, x$posterior_c)
  support <- effect_support(design)
  tv <- crit$theta_tv
  mav <- crit$theta_mav
  limits <- range(effect_quantile(cdf, 5e-4, support),
                  effect_quantile(cdf, 1 - 5e-4, support), tv, mav)
  pad <- 0.04 * diff(limits)
  lower <- max(limits[1L] - pad, support[1L])
  upper <- min(limits[2L] + pad, support[2L])
  grid <- sort(unique(c(seq(lower, upper, length.out = 401L), tv, mav)))
  h <- 1e-4 * (upper - lower)
  dens <- pmax((cdf(grid + h) - cdf(grid - h)) / (2 * h), 0)
  curve <- data.frame(x = grid, y = dens)
  less <- crit$direction == "less"
  in_go <- if (less) grid <= tv else grid >= tv
  in_nogo <- if (less) grid >= mav else grid <= mav
  labels <- c("Go region", "NoGo region")
  regions <- rbind(
    data.frame(x = grid[in_go], y = dens[in_go], region = labels[1L]),
    data.frame(x = grid[in_nogo], y = dens[in_nogo], region = labels[2L])
  )
  regions$region <- factor(regions$region, levels = labels)
  caption <- NULL
  if (!is.null(design$fixed_c)) {
    caption <- if (identical(design$effect, "hazard_ratio")) {
      sprintf("Fixed control median survival time: %s",
              format_num(design$fixed_c, 4L))
    } else {
      sprintf("Fixed control value: %s", format_num(design$fixed_c, 4L))
    }
  }
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = regions,
                         ggplot2::aes(x = .data$x, ymin = 0, ymax = .data$y,
                                      fill = .data$region),
                         alpha = 0.5) +
    ggplot2::geom_line(data = curve, ggplot2::aes(x = .data$x, y = .data$y),
                       linewidth = style$linewidth) +
    ggplot2::geom_vline(xintercept = c(mav, tv), colour = "gray40",
                        linetype = "dotted") +
    ggplot2::annotate("text", x = c(mav, tv), y = Inf,
                      label = c("theta[MAV]", "theta[TV]"), parse = TRUE,
                      hjust = -0.15, vjust = 1.5, colour = "gray30",
                      size = style$label_size) +
    ggplot2::scale_fill_manual(values = stats::setNames(
      unname(cols[c("Go", "NoGo")]), labels
    ), drop = FALSE) +
    # Room above the density for the labels of the thresholds
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.15))
    ) +
    ggplot2::labs(x = effect_label(design), y = "Posterior density",
                  fill = NULL, title = sprintf("Decision: %s", x$decision),
                  subtitle = decision_subtitle(x), caption = caption) +
    style$theme +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
