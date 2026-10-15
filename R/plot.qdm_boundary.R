#' Plot a decision boundary
#'
#' Draws the decision (Go, NoGo, Gray, or Miss if it occurs) at each point of
#' the grid computed by \code{\link{qdm_boundary}} as a tile filled by the
#' decision. The edges between regions of different decisions are the
#' decision boundaries of the design.
#'
#' @param x An object of class \code{"qdm_boundary"}.
#' @param ... Not used.
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
#' plot(qdm_boundary(des))
#' @export
plot.qdm_boundary <- function(x, ...) {
  axes <- attr(x, "axes")
  if (is.null(axes)) {
    stop("'x' must be an unmodified result of qdm_boundary().",
         call. = FALSE)
  }
  d <- data.frame(x = x[[axes$x]], y = x[[axes$y]], decision = x$decision)
  cols <- c(Go = "#658D1B", NoGo = "#D91E49", Gray = "#939597",
            Miss = "#E69F00")
  ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y,
                                  fill = .data$decision)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = cols, drop = TRUE) +
    ggplot2::labs(x = axes$xlab, y = axes$ylab, fill = "Decision",
                  caption = axes$caption) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid = ggplot2::element_blank())
}
