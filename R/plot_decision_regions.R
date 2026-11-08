#' Probabilities of the regions of a decision on two endpoints
#'
#' @param x An object of class \code{"qdm_decision"} for two endpoints.
#' @param cols Named colors of the decisions (\code{decision_colors()}).
#' @param style The style of the plot (\code{plot_style()}).
#'
#' @return A ggplot object with one tile per region, placed as in
#'   \code{qdm_regions()} but with the favorable side of each endpoint at the
#'   right and at the top, labeled with the number and the probability of the
#'   region. The labels of the axes give the conditions of the regions,
#'   preceded by "future" for a predictive criterion.
#' @noRd
plot_decision_regions <- function(x, cols, style) {
  crit <- x$design$criterion
  n <- if (crit$type == "posterior") 3L else 2L
  reg <- qdm_regions(crit$type, crit$direction)
  # Conditions of the axes with the values of the thresholds
  values <- function(text) {
    for (k in 1:2) {
      if (crit$type == "posterior") {
        text <- gsub(paste0("MAV", k), format_num(crit$theta_mav[k], 4L),
                     text, fixed = TRUE)
        text <- gsub(paste0("TV", k), format_num(crit$theta_tv[k], 4L), text,
                     fixed = TRUE)
      } else {
        text <- gsub(paste0("NULL", k), format_num(crit$theta_null[k], 4L),
                     text, fixed = TRUE)
      }
    }
    # A predictive criterion is on the effects observed in the future study
    if (crit$type == "predictive") {
      text <- paste("future", text)
    }
    text
  }
  # Region k is in column i (endpoint 1) and row j (endpoint 2) of
  # qdm_regions(), k = (i - 1) n + j; column 1 and row 1 are the favorable
  # sides
  tiles <- expand.grid(j = seq_len(n), i = seq_len(n))
  tiles$region <- (tiles$i - 1L) * n + tiles$j
  tiles$x <- n + 1L - tiles$i
  tiles$y <- n + 1L - tiles$j
  tiles$prob <- unname(x$prob_regions[paste0("R", tiles$region)])
  role <- rep("Gray", nrow(tiles))
  if (crit$rule == "joint") {
    role[tiles$region %in% crit$go_regions] <- "Go"
    role[tiles$region %in% crit$nogo_regions] <- "NoGo"
  }
  tiles$fill <- shade_color(0.15 + 0.85 * tiles$prob, cols[role])
  tiles$label <- sprintf("R%d\n%s", tiles$region,
                         formatC(tiles$prob, format = "f", digits = 3L))
  ggplot2::ggplot(tiles, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$fill), colour = "white",
                       linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label),
                       size = style$label_size) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_continuous(breaks = seq_len(n),
                                labels = rev(values(colnames(reg))),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_len(n),
                                labels = rev(values(rownames(reg))),
                                expand = c(0, 0)) +
    ggplot2::labs(x = "Endpoint 1", y = "Endpoint 2",
                  title = sprintf("Decision: %s", x$decision),
                  subtitle = decision_subtitle(x)) +
    style$theme +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}
