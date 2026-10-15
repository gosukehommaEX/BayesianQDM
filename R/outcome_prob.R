#' Probability of a set of outcomes under each scenario
#'
#' @param ind Logical vector or matrix indicating the outcomes: over the
#'   outcomes (with \code{sample$weights}), one column per scenario (with
#'   \code{sample$weights} of the same size), or a matrix over pairs of
#'   treatment and control outcomes (with \code{sample$weights_t} and
#'   \code{sample$weights_c}).
#' @param sample A list returned by \code{oc_sample()}.
#'
#' @return A numeric vector with one probability per scenario.
#' @noRd
outcome_prob <- function(ind, sample) {
  if (is.null(sample[["weights"]])) {
    return(colSums(sample$weights_t * (ind %*% sample$weights_c)))
  }
  colSums(sample[["weights"]] * ind)
}
