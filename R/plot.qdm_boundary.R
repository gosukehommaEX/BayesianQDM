#' Plot a decision boundary
#'
#' Draws the decision (Go, NoGo, Gray, or Miss if it occurs) at each point of
#' the grid computed by \code{\link{qdm_boundary}} as a tile filled by the
#' decision. The edges between regions of different decisions are the
#' decision boundaries of the design.
#'
#' The plot is a \code{ggplot} object, so it can be changed further by adding
#' \pkg{ggplot2} functions, for example \code{ggplot2::labs()},
#' \code{ggplot2::theme()} and \code{ggplot2::scale_fill_manual()}.
#'
#' @inheritParams plot.qdm_oc
#' @param x An object of class \code{"qdm_boundary"}.
#' @param ... Not used; arguments given here are ignored with a warning.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' b <- qdm_boundary(des)
#' plot(b)
#' plot(b, col_gray = "gray80", base_size = 14, legend_position = "right")
#' @export
plot.qdm_boundary <- function(x, col_go = "#658D1B", col_nogo = "#D91E49",
                              col_gray = "#939597", col_miss = "#E69F00",
                              base_size = 11, legend_position = "bottom",
                              ...) {
  axes <- attr(x, "axes")
  if (is.null(axes)) {
    stop("'x' must be an unmodified result of qdm_boundary().",
         call. = FALSE)
  }
  warn_unused_dots(list(...), "plot.qdm_boundary")
  cols <- decision_colors(col_go = col_go, col_nogo = col_nogo,
                          col_gray = col_gray, col_miss = col_miss)
  style <- plot_style(base_size = base_size,
                      legend_position = legend_position)
  d <- data.frame(x = x[[axes$x]], y = x[[axes$y]], decision = x$decision)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y,
                                  fill = .data$decision)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = cols, drop = TRUE) +
    ggplot2::labs(x = axes$xlab, y = axes$ylab, fill = "Decision",
                  caption = axes$caption) +
    style$theme +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}
