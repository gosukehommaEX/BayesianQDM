#' Overlay heat map of operating characteristics of two endpoints
#'
#' Draws one heat map with the scenarios on the grid of
#' \code{oc_two_grid()}. Each cell has the color of the decision with the
#' highest probability (\code{dominant_decision()}), shaded in five steps of
#' that probability (below 0.5, 0.5 to 0.6, 0.6 to 0.7, 0.7 to 0.8, and
#' 0.8 or more) from a light shade (weight 0.2 of the color against white)
#' to the color itself, and shows the probability if \code{cell_values} is
#' \code{TRUE}. Scenarios that are not evaluable are white and marked "NE".
#' The legends show the colors of the decisions and the five steps (in
#' gray). The tiles are drawn with \code{show.legend = c(fill = TRUE)}, so
#' that the legend has a key for every decision even if no cell reaches the
#' highest step (from ggplot2 3.5.0, keys are otherwise drawn only for
#' values present in the data), and the keys of the steps use the polygon
#' glyph (the rectangle glyph of ggplot2 4.0.0 ignores the fill).
#'
#' @param x An object of class \code{"qdm_oc"}.
#' @param design The design of \code{x}, with two endpoints.
#' @param cols Colors of the decisions, named by decision.
#' @param style Result of \code{plot_style()}.
#' @param cell_values Logical; write the probabilities in the cells.
#'
#' @return A \code{ggplot} object.
#' @noRd
plot_oc_overlay <- function(x, design, cols, style, cell_values = TRUE) {
  grid <- oc_two_grid(x, design)
  decisions <- intersect(c("Go", "Gray", "NoGo", "Miss"), grid$decisions)
  dom <- dominant_decision(x, decisions)
  weights <- c(0.2, 0.4, 0.6, 0.8, 1)
  steps <- c("< 0.50", "0.50-0.60", "0.60-0.70", "0.70-0.80", ">= 0.80")
  # One fill key per decision and step, and a white key for scenarios that
  # are not evaluable
  keys <- c(paste(rep(decisions, each = 5L), 1:5), "NE")
  key_cols <- c(shade_color(rep(weights, length(decisions)),
                            rep(cols[decisions], each = 5L)), "#FFFFFF")
  names(key_cols) <- keys
  ne <- is.na(dom$decision)
  cells <- data.frame(
    x = grid$x, y = grid$y,
    key = factor(ifelse(ne, "NE", paste(dom$decision, dom$step)),
                 levels = keys),
    label = ifelse(ne, "NE",
                   if (cell_values) sprintf("%.2f", dom$prob) else "")
  )
  # Invisible points that carry the legend of the steps
  legend_steps <- data.frame(x = grid$x[1L], y = grid$y[1L],
                             step = factor(steps, levels = steps))
  gray_steps <- shade_color(weights, "#808080")
  ggplot2::ggplot(cells, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$key), colour = "white",
                       show.legend = c(fill = TRUE, colour = FALSE)) +
    ggplot2::geom_point(data = legend_steps,
                        ggplot2::aes(colour = .data$step), alpha = 0,
                        key_glyph = "polygon") +
    ggplot2::geom_text(ggplot2::aes(label = .data$label),
                       size = style$label_size) +
    ggplot2::scale_fill_manual(values = key_cols, limits = keys,
                               breaks = paste(decisions, 5L),
                               labels = decisions, name = "Decision") +
    ggplot2::scale_colour_manual(values = stats::setNames(gray_steps, steps),
                                 limits = steps, name = "Probability") +
    ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1L),
      colour = ggplot2::guide_legend(
        order = 2L, override.aes = list(fill = gray_steps, alpha = 1)
      )
    ) +
    ggplot2::labs(x = grid$labs[1L], y = grid$labs[2L]) +
    style$theme +
    # The two legends one above the other, so that they fit below the plot
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.box = "vertical")
}
