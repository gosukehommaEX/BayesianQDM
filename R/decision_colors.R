#' Colors of the decisions in a plot
#'
#' Validates the color arguments of a \code{plot()} method and names them by
#' decision.
#'
#' @param ... Colors given as \code{col_go}, \code{col_nogo},
#'   \code{col_gray} and \code{col_miss} (any subset, in any order).
#'
#' @return A named character vector with the names \code{"Go"},
#'   \code{"NoGo"}, \code{"Gray"} and \code{"Miss"} for the arguments given.
#' @noRd
decision_colors <- function(...) {
  cols <- list(...)
  decisions <- c(col_go = "Go", col_nogo = "NoGo", col_gray = "Gray",
                 col_miss = "Miss")
  for (nm in names(cols)) {
    check_color(cols[[nm]], nm)
  }
  stats::setNames(unlist(cols, use.names = FALSE),
                  unname(decisions[names(cols)]))
}
