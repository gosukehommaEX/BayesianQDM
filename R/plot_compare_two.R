#' Heat maps of the operating characteristics of several designs for two
#' endpoints
#'
#' @param rows A list with the rows of each design (data frames).
#' @param designs The list of the designs.
#' @param decisions The decisions to draw.
#' @param cols Named colors of the decisions (\code{decision_colors()}).
#' @param style The style of the plot (\code{plot_style()}).
#' @param cell_values Whether to write the probabilities in the cells.
#'
#' @return A ggplot object with one row of panels per design and one column
#'   per decision.
#' @noRd
plot_compare_two <- function(rows, designs, decisions, cols, style,
                             cell_values) {
  labels <- names(designs)
  grids <- lapply(seq_along(rows), function(k) {
    oc_two_grid(rows[[k]], designs[[k]])
  })
  same <- vapply(grids, function(g) {
    identical(g$labs, grids[[1L]]$labs) && identical(g$x, grids[[1L]]$x) &&
      identical(g$y, grids[[1L]]$y)
  }, logical(1))
  if (!all(same)) {
    stop(paste("The designs do not have the same axes (for example, one has",
               "a concurrent control group and another not); plot them",
               "separately."), call. = FALSE)
  }
  grid <- grids[[1L]]
  n <- nrow(rows[[1L]])
  long <- do.call(rbind, lapply(seq_along(rows), function(k) {
    data.frame(
      x = rep(grid$x, length(decisions)),
      y = rep(grid$y, length(decisions)),
      design = labels[k],
      decision = rep(decisions, each = n),
      prob = unlist(lapply(decisions, function(d) rows[[k]][[d]]),
                    use.names = FALSE)
    )
  }))
  long$design <- factor(long$design, levels = labels)
  long$decision <- factor(long$decision, levels = decisions)
  long$label <- ifelse(is.na(long$prob), "NE",
                       if (cell_values) sprintf("%.2f", long$prob) else "")
  long$fill <- shade_color(long$prob, cols[as.character(long$decision)])
  ggplot2::ggplot(long, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$fill), colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$label),
                       size = style$label_size) +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$design),
                        cols = ggplot2::vars(.data$decision)) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(x = grid$labs[1L], y = grid$labs[2L]) +
    style$theme +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}
