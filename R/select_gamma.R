#' Smallest threshold that meets a target
#'
#' @param grid Increasing numeric vector of thresholds.
#' @param prob Probabilities of meeting the criterion at the thresholds
#'   (non-increasing in the threshold).
#' @param target Upper bound for the probability.
#'
#' @return The smallest element of \code{grid} whose probability is at most
#'   \code{target} (up to a tolerance of 1e-12), or \code{NA} if there is
#'   none.
#' @noRd
select_gamma <- function(grid, prob, target) {
  ok <- which(prob <= target + 1e-12)
  if (length(ok) == 0L) NA_real_ else grid[min(ok)]
}
