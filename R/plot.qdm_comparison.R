#' Plot the operating characteristics of several designs
#'
#' For one endpoint, draws one panel per decision with the probability of
#' that decision against the true treatment effect (or the parameter of the
#' treatment group without a concurrent control group, as in
#' \code{\link{plot.qdm_oc}}), with one curve per design. For two endpoints,
#' draws heat maps of the probabilities, with one row of panels per design
#' and one column per decision.
#'
#' The horizontal axis must be the same for all designs: the scenarios are
#' the same, but a design without a concurrent control group has the
#' parameter of the treatment group on the axis, and the thresholds are
#' shifted by the location of the control group (see
#' \code{\link{plot.qdm_oc}}). The thresholds of the criteria are drawn when
#' they are the same for all designs.
#'
#' @param x An object of class \code{"qdm_comparison"}.
#' @param decisions The decisions to draw, some of \code{"Go"},
#'   \code{"NoGo"}, \code{"Gray"} and \code{"Miss"}, or \code{NULL} (the
#'   default) for Go, NoGo and Gray, and Miss if it has positive probability
#'   in some scenario.
#' @param col_design Colors of the curves of the designs (one endpoint), one
#'   per design, or \code{NULL} for the palette \code{"Dark 3"} of
#'   \code{grDevices::hcl.colors()}.
#' @param col_go,col_nogo,col_gray,col_miss Colors of the heat maps of the
#'   decisions (two endpoints).
#' @param linewidth,point_size,label_size,base_size,legend_position As in
#'   \code{\link{plot.qdm_oc}}.
#' @param cell_values Two endpoints only: \code{TRUE} (the default) to write
#'   the probability in each cell of the heat maps.
#' @param ... Not used; arguments given here are ignored with a warning.
#'
#' @return A ggplot object.
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' cmp <- qdm_compare(list(n12 = des, n24 = qdm_sample_size(des, 24, 24)),
#'                    data.frame(pi_t = seq(0.2, 0.8, by = 0.1), pi_c = 0.3))
#' plot(cmp)
#' plot(cmp, decisions = "Go")
#' @export
plot.qdm_comparison <- function(x, decisions = NULL, col_design = NULL,
                                col_go = "#658D1B", col_nogo = "#D91E49",
                                col_gray = "#939597", col_miss = "#E69F00",
                                linewidth = 1, point_size = 1.5,
                                label_size = NULL, base_size = 11,
                                legend_position = "bottom",
                                cell_values = TRUE, ...) {
  designs <- attr(x, "designs")
  if (is.null(designs)) {
    stop("'x' must be an unmodified result of qdm_compare().", call. = FALSE)
  }
  warn_unused_dots(list(...), "plot.qdm_comparison")
  cols <- decision_colors(col_go = col_go, col_nogo = col_nogo,
                          col_gray = col_gray, col_miss = col_miss)
  style <- plot_style(linewidth = linewidth, point_size = point_size,
                      label_size = label_size, base_size = base_size,
                      legend_position = legend_position)
  check_flag(cell_values, "cell_values")
  all_decisions <- c("Go", "NoGo", "Gray", "Miss")
  if (is.null(decisions)) {
    decisions <- if (all(x$Miss <= 0, na.rm = TRUE)) {
      all_decisions[1:3]
    } else {
      all_decisions
    }
  } else if (!is.character(decisions) || length(decisions) == 0L ||
             anyNA(decisions) || !all(decisions %in% all_decisions) ||
             anyDuplicated(decisions) > 0L) {
    stop(paste("'decisions' must contain some of \"Go\", \"NoGo\", \"Gray\"",
               "and \"Miss\"."), call. = FALSE)
  } else {
    decisions <- intersect(all_decisions, decisions)
  }
  labels <- levels(x$design)
  rows <- lapply(labels, function(lab) x[x$design == lab, , drop = FALSE])
  if (designs[[1L]]$criterion$n_endpoints == 2L) {
    return(plot_compare_two(rows, designs, decisions, cols, style,
                            cell_values))
  }
  if (!cell_values) {
    stop("'cell_values' applies only to the heat maps of two endpoints.",
         call. = FALSE)
  }
  if (is.null(col_design)) {
    col_design <- grDevices::hcl.colors(length(labels), "Dark 3")
  } else {
    if (!is.character(col_design) || length(col_design) != length(labels)) {
      stop(sprintf("'col_design' must give one color for each of the %d %s.",
                   length(labels),
                   if (length(labels) == 1L) "design" else "designs"),
           call. = FALSE)
    }
    for (k in seq_along(col_design)) {
      check_color(col_design[k], "col_design")
    }
  }
  axes <- lapply(seq_along(labels), function(k) {
    oc_axis(designs[[k]], rows[[k]])
  })
  same <- vapply(axes, function(a) {
    isTRUE(all.equal(a$x, axes[[1L]]$x)) &&
      identical(deparse(a$xlab), deparse(axes[[1L]]$xlab))
  }, logical(1))
  if (!all(same)) {
    stop(paste("The designs do not have the same horizontal axis (for",
               "example, one has a concurrent control group and another",
               "not); plot them separately."), call. = FALSE)
  }
  axis <- axes[[1L]]
  if (anyDuplicated(round(axis$x, 10)) > 0L) {
    stop(paste("The scenarios must have distinct values on the horizontal",
               "axis; plot scenarios with the same value separately."),
         call. = FALSE)
  }
  same_lines <- all(vapply(axes, function(a) {
    isTRUE(all.equal(a$lines, axis$lines)) && identical(a$labels, axis$labels)
  }, logical(1)))
  ord <- order(axis$x)
  long <- do.call(rbind, lapply(seq_along(labels), function(k) {
    data.frame(
      x = rep(axis$x[ord], length(decisions)),
      prob = unlist(lapply(decisions, function(d) rows[[k]][[d]][ord]),
                    use.names = FALSE),
      design = labels[k],
      decision = rep(decisions, each = length(ord))
    )
  }))
  long$design <- factor(long$design, levels = labels)
  long$decision <- factor(long$decision, levels = decisions)
  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$x, y = .data$prob,
                                          colour = .data$design))
  if (same_lines) {
    p <- p +
      ggplot2::geom_vline(xintercept = axis$lines, colour = "gray40",
                          linetype = "dotted") +
      # A layer with its own data: annotate() with several labels fails when
      # the plot has facets. The labels are at different heights, so that
      # they do not overlap in narrow panels.
      ggplot2::geom_text(data = data.frame(
                           x = axis$lines, label = axis$labels,
                           y = 1 - 0.07 * (seq_along(axis$lines) - 1)
                         ),
                         ggplot2::aes(x = .data$x, y = .data$y,
                                      label = .data$label),
                         parse = TRUE, hjust = -0.15, vjust = 1,
                         colour = "gray30", size = style$label_size,
                         inherit.aes = FALSE)
  }
  p +
    ggplot2::geom_line(linewidth = style$linewidth) +
    ggplot2::geom_point(size = style$point_size) +
    ggplot2::facet_wrap(ggplot2::vars(.data$decision)) +
    ggplot2::scale_colour_manual(values = stats::setNames(col_design,
                                                          labels)) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                breaks = seq(0, 1, by = 0.2)) +
    ggplot2::labs(x = axis$xlab, y = "Probability", colour = "Design",
                  caption = if (same_lines) axis$caption else NULL) +
    style$theme +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
