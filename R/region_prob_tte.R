#' Probabilities of the Go and NoGo regions for a time-to-event endpoint
#'
#' For the difference of the milestone survival probabilities, posterior
#' probabilities are computed exactly from the Gamma distribution of the
#' treatment hazard for a fixed control and by \code{psurvdiff()} otherwise.
#' For the hazard ratio, they are computed exactly from the Gamma
#' distribution of the treatment hazard for a fixed control (whose hazard is
#' \eqn{\log(2)} divided by the fixed median survival time) and by
#' \code{phazardratio()} otherwise. Posterior predictive probabilities are
#' estimated by simulating \code{design$nMC} future studies for each
#' posterior distribution (\code{cpp_pred_tte()}); the Go and NoGo
#' probabilities are estimated from the same simulations and add up to 1.
#'
#' With \code{direction = "less"}, the Go region is
#' \eqn{\theta < \theta_{TV}} (\eqn{\tilde\theta < \theta_{NULL}}) and the
#' NoGo region is \eqn{\theta \ge \theta_{MAV}}
#' (\eqn{\tilde\theta \ge \theta_{NULL}}). The distributions of the effects
#' are continuous, so the probabilities are the other tails of the same
#' distribution functions.
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
  less <- crit$direction == "less"
  hr <- identical(design$effect, "hazard_ratio")
  t0 <- design$milestone
  fixed <- design$fixed_c
  if (crit$type == "posterior") {
    # Distribution function of the effect, P(theta <= q), or its upper tail
    cdf <- function(q, lower.tail) {
      if (hr && !is.null(fixed)) {
        # lambda_t / lambda_c0 <= q if and only if lambda_t <= q lambda_c0
        stats::pgamma(q * log(2) / fixed, post_t$shape, post_t$rate,
                      lower.tail = lower.tail)
      } else if (hr) {
        phazardratio(q, post_t$shape, post_t$rate, post_c$shape,
                     post_c$rate, lower.tail = lower.tail)
      } else if (!is.null(fixed)) {
        # S_t - fixed <= q if and only if S_t <= q + fixed
        psurv_fixed(q + fixed, t0, post_t$shape, post_t$rate,
                    lower.tail = lower.tail)
      } else {
        psurvdiff(q, t0, post_t$shape, post_t$rate, post_c$shape,
                  post_c$rate, lower.tail = lower.tail)
      }
    }
    # Go: P(theta > TV), or P(theta < TV) for direction = "less"; NoGo:
    # P(theta <= MAV), or P(theta >= MAV)
    return(list(go = cdf(crit$theta_tv, lower.tail = less),
                nogo = cdf(crit$theta_mav, lower.tail = !less)))
  }
  n <- length(post_t$shape)
  ca <- conduct_args(design$conduct_future)
  prior_c <- if (is.null(fixed)) {
    c(design$prior_c$shape, design$prior_c$rate)
  } else {
    c(1, 1)
  }
  # Fixed control: milestone survival probability, or hazard for the hazard
  # ratio
  value_c <- if (is.null(fixed)) 0 else if (hr) log(2) / fixed else fixed
  above <- cpp_pred_tte(
    as.numeric(post_t$shape), as.numeric(post_t$rate),
    if (is.null(fixed)) rep_len(as.numeric(post_c$shape), n) else numeric(0),
    if (is.null(fixed)) rep_len(as.numeric(post_c$rate), n) else numeric(0),
    value_c, c(design$prior_t$shape, design$prior_t$rate), prior_c,
    crit$m_t, crit$m_c, ca$acc_time, ca$acc_cum, ca$drop_start,
    ca$drop_hazard, ca$cut_time, ca$cut_events, if (hr) 1 else t0,
    crit$theta_null, design$nMC, if (hr) 1L else 0L
  )
  # P(future effect > NULL); its complement is P(future effect <= NULL)
  if (less) {
    list(go = 1 - above, nogo = above)
  } else {
    list(go = above, nogo = 1 - above)
  }
}
