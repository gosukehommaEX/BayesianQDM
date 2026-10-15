#' Simulated counts of two binary endpoints with common random numbers
#'
#' Assigns each simulated subject to one of the four cells by comparing a
#' uniform random number with the cumulative cell probabilities, so that the
#' same random numbers give similar data under similar scenarios.
#'
#' @param u Numeric matrix of uniform random numbers, one row per simulated
#'   study and one column per subject.
#' @param p Cell probabilities (\code{"00"}, \code{"01"}, \code{"10"},
#'   \code{"11"}).
#'
#' @return An integer matrix of counts with one row per simulated study and
#'   the columns \code{"00"}, \code{"01"}, \code{"10"} and \code{"11"}.
#' @noRd
simulate_counts_binary2 <- function(u, p) {
  cum <- cumsum(p)[1:3]
  cell <- 1L + (u > cum[1]) + (u > cum[2]) + (u > cum[3])
  out <- vapply(1:4, function(k) rowSums(cell == k), numeric(nrow(u)))
  out <- matrix(as.integer(out), nrow = nrow(u))
  dimnames(out) <- list(NULL, c("00", "01", "10", "11"))
  out
}
