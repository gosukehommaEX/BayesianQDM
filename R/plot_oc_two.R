#' Heat maps of operating characteristics of two endpoints
#'
#' Draws one panel per decision with the scenarios on a grid of the true
#' treatment effects of the two endpoints (or the true parameters of the
#' treatment group for a design without a concurrent control group). Each
#' cell is shaded from white (probability 0) to the color of the decision
#' (probability 1) and shows the probability of the decision if
#' \code{cell_values} is \code{TRUE}; scenarios that are not evaluable are
#' white and marked "NE".
#'
#' @param x An object of class \code{"qdm_oc"}.
#' @param design The design of \code{x}, with two endpoints.
#' @param cols Colors of the decisions, named by decision.
#' @param style Result of \code{plot_style()}.
#' @param cell_values Logical; write the probabilities in the cells.
#'
#' @return A \code{ggplot} object.
#' @noRd
plot_oc_two <- function(x, design, cols, style, cell_values = TRUE) {
  grid <- oc_two_grid(x, design)
  decisions <- grid$decisions
  long <- data.frame(
    x = rep(grid$x, length(decisions)),
    y = rep(grid$y, length(decisions)),
    decision = factor(rep(decisions, each = nrow(x)), levels = decisions),
    prob = unlist(lapply(decisions, function(d) x[[d]]), use.names = FALSE)
  )
  long$label <- ifelse(is.na(long$prob), "NE",
                       if (cell_values) sprintf("%.2f", long$prob) else "")
  long$fill <- shade_color(long$prob, cols[as.character(long$decision)])
  ggplot2::ggplot(long, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$fill), colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$label),
                       size = style$label_size) +
    ggplot2::facet_wrap(ggplot2::vars(.data$decision)) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(x = grid$labs[1L], y = grid$labs[2L]) +
    style$theme +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}
