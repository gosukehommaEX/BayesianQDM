#' Decision probabilities from outcomes and their probabilities
#'
#' @param sample A list returned by \code{oc_sample()}.
#' @param gamma_go,gamma_nogo Thresholds of the Go and NoGo criteria.
#' @param rule One of \code{"joint"}, \code{"and"} or \code{"or"}.
#'
#' @return A matrix with one row per scenario and columns \code{Go},
#'   \code{NoGo}, \code{Gray} and \code{Miss} (\code{NA} for scenarios that
#'   are not evaluable).
#' @noRd
oc_from_sample <- function(sample, gamma_go, gamma_nogo, rule = "joint") {
  go <- meets_criterion(sample$p_go, gamma_go, rule)
  nogo <- meets_criterion(sample$p_nogo, gamma_nogo, rule)
  out <- cbind(Go = outcome_prob(go & !nogo, sample),
               NoGo = outcome_prob(!go & nogo, sample),
               Gray = outcome_prob(!go & !nogo, sample),
               Miss = outcome_prob(go & nogo, sample))
  pmin(pmax(out, 0), 1)
}
