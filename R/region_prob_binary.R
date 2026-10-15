#' Probabilities of the Go and NoGo regions for a binary endpoint
#'
#' Computes the posterior (or posterior predictive) probabilities of the Go
#' and NoGo regions of a single binary endpoint from the parameters of the
#' Beta posterior distributions of the two groups, or of the treatment group
#' and a fixed control response rate.
#'
#' @param criterion An object of class \code{"qdm_criterion"} for one
#'   endpoint.
#' @param a_t,b_t,a_c,b_c Numeric vectors of posterior Beta parameters
#'   (\code{a_c} and \code{b_c} are not used for a fixed control).
#' @param fixed Fixed control response rate, or \code{NULL}.
#'
#' @return A list with numeric vectors \code{go} and \code{nogo}.
#' @noRd
region_prob_binary <- function(criterion, a_t, b_t, a_c, b_c, fixed = NULL) {
  if (!is.null(fixed)) {
    if (criterion$type == "posterior") {
      # pi_t - fixed > TV if and only if pi_t > TV + fixed
      go <- stats::pbeta(criterion$theta_tv + fixed, a_t, b_t,
                         lower.tail = FALSE)
      nogo <- stats::pbeta(criterion$theta_mav + fixed, a_t, b_t)
    } else {
      go <- ppred_fixed_binary(criterion$theta_null, criterion$m_t, a_t, b_t,
                               fixed, lower.tail = FALSE)
      nogo <- ppred_fixed_binary(criterion$theta_null, criterion$m_t, a_t,
                                 b_t, fixed)
    }
    return(list(go = go, nogo = nogo))
  }
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
