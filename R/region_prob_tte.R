#' Probabilities of the Go and NoGo regions for a time-to-event endpoint
#'
#' Posterior probabilities are computed exactly from the Gamma distribution
#' of the treatment hazard for a fixed control and by \code{psurvdiff()}
#' otherwise. Posterior predictive probabilities are estimated by
#' simulating \code{design$nMC} future studies for each posterior
#' distribution (\code{cpp_pred_tte()}); the Go and NoGo probabilities are
#' estimated from the same simulations and add up to 1.
#'
#' @param design An object of class \code{"qdm_design_tte"}.
#' @param post_t List with numeric vectors \code{shape} and \code{rate} of the
#'   posterior distributions of the treatment hazard.
#' @param post_c The same for the control hazard (vectors of the same length
#'   as those of \code{post_t}, or of length 1), or \code{NULL} for a fixed
#'   control.
#'
#' @return A list with numeric vectors \code{go} and \code{nogo}.
#' @noRd
region_prob_tte <- function(design, post_t, post_c) {
  crit <- design$criterion
  t0 <- design$milestone
  fixed <- design$fixed_c
  if (crit$type == "posterior") {
    if (!is.null(fixed)) {
      # S_t - fixed > q if and only if S_t > q + fixed
      go <- psurv_fixed(crit$theta_tv + fixed, t0, post_t$shape,
                        post_t$rate, lower.tail = FALSE)
      nogo <- psurv_fixed(crit$theta_mav + fixed, t0, post_t$shape,
                          post_t$rate)
    } else {
      go <- psurvdiff(crit$theta_tv, t0, post_t$shape, post_t$rate,
                      post_c$shape, post_c$rate, lower.tail = FALSE)
      nogo <- psurvdiff(crit$theta_mav, t0, post_t$shape, post_t$rate,
                        post_c$shape, post_c$rate)
    }
    return(list(go = go, nogo = nogo))
  }
  n <- length(post_t$shape)
  ca <- conduct_args(design$conduct_future)
  prior_c <- if (is.null(fixed)) {
    c(design$prior_c$shape, design$prior_c$rate)
  } else {
    c(1, 1)
  }
  go <- cpp_pred_tte(
    as.numeric(post_t$shape), as.numeric(post_t$rate),
    if (is.null(fixed)) rep_len(as.numeric(post_c$shape), n) else numeric(0),
    if (is.null(fixed)) rep_len(as.numeric(post_c$rate), n) else numeric(0),
    if (is.null(fixed)) 0 else fixed,
    c(design$prior_t$shape, design$prior_t$rate), prior_c, crit$m_t,
    crit$m_c, ca$acc_time, ca$acc_cum, ca$drop_start, ca$drop_hazard,
    ca$cut_time, ca$cut_events, t0, crit$theta_null, design$nMC
  )
  list(go = go, nogo = 1 - go)
}
