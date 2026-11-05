#' Posterior distribution function of the treatment effect
#'
#' The distribution function of the treatment effect of one endpoint given
#' the data, from the parameters of the posterior distributions of a
#' decision: the difference of Beta variables (binary), of t variables
#' (continuous, by numerical integration) or of milestone survival
#' probabilities (time-to-event), or the ratio of Gamma variables (hazard
#' ratio). For a fixed control, the effect is the parameter of the treatment
#' group minus the fixed value (or its hazard divided by the fixed hazard).
#'
#' @param design A design for one endpoint.
#' @param post_t,post_c The elements \code{posterior_t} and
#'   \code{posterior_c} of a result of \code{qdm_decide()}.
#'
#' @return A function of a numeric vector \code{q} that gives
#'   \eqn{\Pr(\theta \le q \mid D)}.
#' @noRd
decision_cdf <- function(design, post_t, post_c) {
  fixed <- design$fixed_c
  if (design$endpoint == "binary") {
    if (!is.null(fixed)) {
      return(function(q) stats::pbeta(q + fixed, post_t[["a"]], post_t[["b"]]))
    }
    return(function(q) {
      pbetadiff(q, post_t[["a"]], post_t[["b"]], post_c[["a"]], post_c[["b"]])
    })
  }
  if (design$endpoint == "continuous") {
    # The mean has the distribution t_nu(mu, sigma^2 / kappa)
    scale_t <- post_t[["sigma"]] / sqrt(post_t[["kappa"]])
    if (!is.null(fixed)) {
      return(function(q) {
        stats::pt((q + fixed - post_t[["mu"]]) / scale_t, post_t[["nu"]])
      })
    }
    scale_c <- post_c[["sigma"]] / sqrt(post_c[["kappa"]])
    return(function(q) {
      ptdiff(q, post_t[["mu"]], scale_t, post_t[["nu"]], post_c[["mu"]],
             scale_c, post_c[["nu"]], method = "NI")
    })
  }
  shape_t <- post_t[["shape"]]
  rate_t <- post_t[["rate"]]
  if (identical(design$effect, "hazard_ratio")) {
    if (!is.null(fixed)) {
      # lambda_t / lambda_c0 <= q if and only if lambda_t <= q lambda_c0
      return(function(q) {
        stats::pgamma(pmax(q, 0) * log(2) / fixed, shape_t, rate_t)
      })
    }
    return(function(q) {
      phazardratio(q, shape_t, rate_t, post_c[["shape"]], post_c[["rate"]])
    })
  }
  t0 <- design$milestone
  if (!is.null(fixed)) {
    return(function(q) psurv_fixed(q + fixed, t0, shape_t, rate_t))
  }
  function(q) {
    psurvdiff(q, t0, shape_t, rate_t, post_c[["shape"]], post_c[["rate"]])
  }
}
