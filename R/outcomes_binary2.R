#' All outcomes of one group for two binary endpoints
#'
#' Lists all vectors of counts of the four joint response patterns
#' (\code{"00"}, \code{"01"}, \code{"10"}, \code{"11"}) with total \code{n}.
#' There are \code{choose(n + 3, 3)} of them.
#'
#' @param n Number of subjects (a whole number of at least 0).
#'
#' @return An integer matrix with one row per outcome and the columns
#'   \code{"00"}, \code{"01"}, \code{"10"} and \code{"11"}.
#' @noRd
outcomes_binary2 <- function(n) {
  g <- expand.grid(x01 = 0:n, x10 = 0:n, x11 = 0:n)
  g <- g[rowSums(g) <= n, , drop = FALSE]
  out <- cbind(n - rowSums(g), g$x01, g$x10, g$x11)
  storage.mode(out) <- "integer"
  dimnames(out) <- list(NULL, c("00", "01", "10", "11"))
  out
}
