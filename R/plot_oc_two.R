#' Heat maps of operating characteristics of two endpoints
#'
#' Draws one panel per decision with the scenarios on a grid of the true
#' treatment effects of the two endpoints (or the true parameters of the
#' treatment group for a design without a concurrent control group). Each
#' cell shows the probability of the decision and is shaded from white
#' (probability 0) to the color of the decision (probability 1); scenarios
#' that are not evaluable are white and marked "NE".
#'
#' @param x An object of class \code{"qdm_oc"}.
#' @param design The design of \code{x}, with two endpoints.
#' @param cols Colors of the decisions, named by decision.
#' @param style Result of \code{plot_style()}.
#'
#' @return A \code{ggplot} object.
#' @noRd
plot_oc_two <- function(x, design, cols, style) {
  par <- switch(design$endpoint,
                binary2 = c("pi_t1", "pi_t2", "pi_c1", "pi_c2"),
                continuous2 = c("mu_t1", "mu_t2", "mu_c1", "mu_c2"))
  concurrent <- design$control == "concurrent"
  if (concurrent) {
    ax1 <- x[[par[1L]]] - x[[par[3L]]]
    ax2 <- x[[par[2L]]] - x[[par[4L]]]
    labs <- c("Treatment effect on endpoint 1",
              "Treatment effect on endpoint 2")
  } else {
    ax1 <- x[[par[1L]]]
    ax2 <- x[[par[2L]]]
    labs <- c("Treatment group, endpoint 1", "Treatment group, endpoint 2")
  }
  ax1 <- round(ax1, 10)
  ax2 <- round(ax2, 10)
  if (anyDuplicated(data.frame(ax1, ax2)) > 0L) {
    stop(paste("The scenarios must have distinct values on the axes; plot",
               "the scenarios with different correlations or control",
               "parameters separately."), call. = FALSE)
  }
  decisions <- c("Go", "NoGo", "Gray", "Miss")
  if (all(x$Miss <= 0, na.rm = TRUE)) {
    decisions <- decisions[1:3]
  }
  lev1 <- sort(unique(ax1))
  lev2 <- sort(unique(ax2))
  long <- data.frame(
    x = factor(rep(ax1, length(decisions)), levels = lev1,
               labels = format(lev1, trim = TRUE)),
    y = factor(rep(ax2, length(decisions)), levels = lev2,
               labels = format(lev2, trim = TRUE)),
    decision = factor(rep(decisions, each = nrow(x)), levels = decisions),
    prob = unlist(lapply(decisions, function(d) x[[d]]), use.names = FALSE)
  )
  long$label <- ifelse(is.na(long$prob), "NE", sprintf("%.2f", long$prob))
  long$fill <- shade_color(long$prob, cols[as.character(long$decision)])
  ggplot2::ggplot(long, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$fill), colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = .data$label),
                       size = style$label_size) +
    ggplot2::facet_wrap(ggplot2::vars(.data$decision)) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(x = labs[1L], y = labs[2L]) +
    style$theme +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}
