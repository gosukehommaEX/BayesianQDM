#' Decision probabilities from outcomes and their probabilities
#'
#' @param sample A list returned by \code{oc_sample()}.
#' @param gamma_go,gamma_nogo Thresholds of the Go and NoGo criteria.
#'
#' @return A matrix with one row per scenario and columns \code{Go},
#'   \code{NoGo}, \code{Gray} and \code{Miss}.
#' @noRd
oc_from_sample <- function(sample, gamma_go, gamma_nogo) {
  go <- meets_threshold(sample$p_go, gamma_go)
  nogo <- meets_threshold(sample$p_nogo, gamma_nogo)
  w <- sample$weights
  out <- cbind(Go = colSums(w * (go & !nogo)),
               NoGo = colSums(w * (!go & nogo)),
               Gray = colSums(w * (!go & !nogo)),
               Miss = colSums(w * (go & nogo)))
  pmin(pmax(out, 0), 1)
}
