#' Compare probabilities with a decision threshold
#'
#' A probability meets a threshold if it is at least the threshold. Values
#' below the threshold by less than \code{tol} are treated as equal to it, so
#' that ties that arise in symmetric settings (for example a posterior
#' probability of exactly 0.5) are not decided by rounding errors.
#'
#' @param p Numeric vector or matrix of probabilities.
#' @param gamma A single threshold.
#' @param tol Tolerance.
#'
#' @return A logical vector or matrix of the same shape as \code{p}.
#' @noRd
meets_threshold <- function(p, gamma, tol = 1e-10) {
  p >= gamma - tol
}
