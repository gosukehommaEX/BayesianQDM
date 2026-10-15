#' Probabilities of the outcomes under one scenario
#'
#' @param sample A list returned by \code{oc_sample()}.
#' @param j Index of the scenario.
#'
#' @return A numeric vector of probabilities, or for enumerated outcomes of
#'   two groups a list with elements \code{t} and \code{c} (probabilities of
#'   the treatment and control outcomes).
#' @noRd
sample_weights <- function(sample, j) {
  if (is.null(sample[["weights"]])) {
    return(list(t = sample$weights_t[, j], c = sample$weights_c[, j]))
  }
  sample[["weights"]][, j]
}
