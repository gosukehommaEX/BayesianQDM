#' Grid of the scenarios of two endpoints
#'
#' Places the scenarios of operating characteristics for two endpoints on a
#' grid of the true treatment effects of the two endpoints (or the true
#' parameters of the treatment group without a concurrent control group),
#' as used by the heat maps of \code{plot.qdm_oc()}.
#'
#' @param x An object of class \code{"qdm_oc"}.
#' @param design The design of \code{x}, with two endpoints.
#'
#' @return A list with \code{x} and \code{y}, factors with one element per
#'   scenario (the values on the axes, with their sorted unique values as
#'   levels), \code{labs}, the labels of the axes, and \code{decisions},
#'   the decisions shown (Miss only if it has positive probability).
#' @noRd
oc_two_grid <- function(x, design) {
  par <- switch(design$endpoint,
                binary2 = c("pi_t1", "pi_t2", "pi_c1", "pi_c2"),
                continuous2 = c("mu_t1", "mu_t2", "mu_c1", "mu_c2"))
  if (design$control == "concurrent") {
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
  # Labels of the levels, with more digits if 7 significant digits do not
  # distinguish them (equal labels would merge the levels)
  level_labels <- function(lev) {
    lab <- format(lev, trim = TRUE)
    if (anyDuplicated(lab) > 0L) {
      lab <- format(lev, trim = TRUE, digits = 15L)
    }
    lab
  }
  lev1 <- sort(unique(ax1))
  lev2 <- sort(unique(ax2))
  list(x = factor(ax1, levels = lev1, labels = level_labels(lev1)),
       y = factor(ax2, levels = lev2, labels = level_labels(lev2)),
       labs = labs, decisions = decisions)
}
