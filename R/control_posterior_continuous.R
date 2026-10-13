#' Distribution of the control mean for a continuous endpoint
#'
#' Returns the parameters of the distribution of the control mean of a
#' continuous design in the form returned by \code{posterior_continuous()},
#' so that the mean has the distribution
#' \eqn{t_{\nu}(\mu, \sigma^2 / \kappa)}. For a concurrent control group it is
#' the posterior distribution given the external control data (if any) and
#' the current control data; for an external control it is the posterior
#' distribution given the external data only. For a hypothetical control
#' (Kang et al., 2026, Section 4.2.4) it is
#' \eqn{\mu_c \sim t_{\nu_{nt}}(\mu_{0c}, r \sigma_{nt}^2 / \kappa_{nt})},
#' given by \eqn{\mu = \mu_{0c}}, \eqn{\kappa = \kappa_{nt}},
#' \eqn{\nu = \nu_{nt}} and \eqn{\sigma^2 = r \sigma_{nt}^2}; with these
#' parameters the posterior predictive distribution of the mean of \eqn{m_c}
#' future control subjects is
#' \eqn{t_{\nu_{nt}}(\mu_{0c}, r \sigma_{nt}^2 (1 / \kappa_{nt} + 1 / m_c))}.
#'
#' @param design An object of class \code{"qdm_design_continuous"}.
#' @param post_t Posterior parameters of the treatment group, as returned by
#'   \code{posterior_continuous()}.
#' @param n,mean,sd Current control data (concurrent control only).
#'
#' @return A list with elements \code{mu}, \code{kappa}, \code{nu} and
#'   \code{sigma2}.
#' @noRd
control_posterior_continuous <- function(design, post_t, n = 0L, mean = 0,
                                         sd = 0) {
  switch(
    design$control,
    concurrent = posterior_continuous(design$prior_c, n, mean, sd,
                                      external = design$external_c),
    external = posterior_continuous(design$prior_c, 0L,
                                    external = design$external_c),
    hypothetical = list(mu = design$hypothetical_c$mean,
                        kappa = post_t$kappa, nu = post_t$nu,
                        sigma2 = design$hypothetical_c$r * post_t$sigma2)
  )
}
