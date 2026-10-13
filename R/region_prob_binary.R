#' Probabilities of the Go and NoGo regions for a binary endpoint
#'
#' Computes the posterior (or posterior predictive) probabilities of the Go
#' and NoGo regions of a single binary endpoint from the parameters of the
#' Beta posterior distributions of the two groups.
#'
#' @param criterion An object of class \code{"qdm_criterion"} for one
#'   endpoint.
#' @param a_t,b_t,a_c,b_c Numeric vectors of posterior Beta parameters.
#'
#' @return A list with numeric vectors \code{go} and \code{nogo}.
#' @noRd
region_prob_binary <- function(criterion, a_t, b_t, a_c, b_c) {
  if (criterion$type == "posterior") {
    go <- pbetadiff(criterion$theta_tv, a_t, b_t, a_c, b_c,
                    lower.tail = FALSE)
    nogo <- pbetadiff(criterion$theta_mav, a_t, b_t, a_c, b_c)
  } else {
    go <- pbetabinomdiff(criterion$theta_null, criterion$m_t, a_t, b_t,
                         criterion$m_c, a_c, b_c, lower.tail = FALSE)
    nogo <- pbetabinomdiff(criterion$theta_null, criterion$m_t, a_t, b_t,
                           criterion$m_c, a_c, b_c)
  }
  list(go = go, nogo = nogo)
}
