#' Probabilities of the Go and NoGo regions for a continuous endpoint
#'
#' Computes the posterior (or posterior predictive) probabilities of the Go
#' and NoGo regions of a single continuous endpoint from the parameters of the
#' Normal-Inverse-chi-squared posterior distributions of the two groups. The
#' posterior distribution of the mean of a group is
#' \eqn{t_{\nu}(\mu, \sigma^2 / \kappa)}, and the posterior predictive
#' distribution of the mean of \eqn{m} future subjects is
#' \eqn{t_{\nu}(\mu, \sigma^2 (1 / \kappa + 1 / m))}. With
#' \code{method = "MC"}, both probabilities are estimated from the same draws,
#' so that they add up to 1 for the predictive criterion.
#'
#' @param criterion An object of class \code{"qdm_criterion"} for one
#'   endpoint.
#' @param post_t,post_c Lists with elements \code{mu}, \code{kappa},
#'   \code{nu} and \code{sigma2}.
#' @param method One of \code{"NI"}, \code{"MM"} or \code{"MC"}.
#' @param nMC Number of Monte Carlo draws for \code{method = "MC"}.
#'
#' @return A list with numeric vectors \code{go} and \code{nogo}.
#' @noRd
region_prob_continuous <- function(criterion, post_t, post_c, method,
                                   nMC = 10000L) {
  if (criterion$type == "posterior") {
    scale_t <- sqrt(post_t$sigma2 / post_t$kappa)
    scale_c <- sqrt(post_c$sigma2 / post_c$kappa)
    q_go <- criterion$theta_tv
    q_nogo <- criterion$theta_mav
  } else {
    scale_t <- sqrt(post_t$sigma2 * (1 / post_t$kappa + 1 / criterion$m_t))
    scale_c <- sqrt(post_c$sigma2 * (1 / post_c$kappa + 1 / criterion$m_c))
    q_go <- criterion$theta_null
    q_nogo <- criterion$theta_null
  }
  if (method == "MC") {
    a <- recycle_common(list(q_go = q_go, q_nogo = q_nogo, mu_t = post_t$mu,
                             scale_t = scale_t, df_t = post_t$nu,
                             mu_c = post_c$mu, scale_c = scale_c,
                             df_c = post_c$nu))
    p <- cpp_ptdiff_mc_regions(a$q_go, a$q_nogo, a$mu_t, a$scale_t, a$df_t,
                               a$mu_c, a$scale_c, a$df_c, as.integer(nMC))
    return(list(go = p[, 1L], nogo = p[, 2L]))
  }
  go <- ptdiff(q_go, post_t$mu, scale_t, post_t$nu, post_c$mu, scale_c,
               post_c$nu, method = method, lower.tail = FALSE)
  nogo <- ptdiff(q_nogo, post_t$mu, scale_t, post_t$nu, post_c$mu, scale_c,
                 post_c$nu, method = method)
  list(go = go, nogo = nogo)
}
