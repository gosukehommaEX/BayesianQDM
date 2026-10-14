#' Whether a Go or NoGo criterion is met
#'
#' @param p Probabilities: a vector or matrix for one probability, or a list
#'   of two (one per endpoint) for the marginal rules. For a single decision
#'   under a marginal rule, a numeric vector of length 2 is also accepted.
#' @param gamma Threshold(s): one number, or one per endpoint.
#' @param rule One of \code{"joint"}, \code{"and"} or \code{"or"}.
#'
#' @return A logical vector or matrix (see \code{meets_threshold()}).
#' @noRd
meets_criterion <- function(p, gamma, rule = "joint") {
  if (rule == "joint") {
    return(meets_threshold(p, gamma))
  }
  if (!is.list(p)) {
    p <- as.list(p)
  }
  a <- meets_threshold(p[[1L]], gamma[1L])
  b <- meets_threshold(p[[2L]], gamma[2L])
  if (rule == "and") a & b else a | b
}
