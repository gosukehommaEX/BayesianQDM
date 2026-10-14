#' Parameters of the distribution of the control group of two continuous
#' endpoints
#'
#' For a concurrent control group, the posterior parameters given the
#' external and current control data; for external control data only, the
#' parameters given the external data; for a hypothetical control, the
#' pseudo-parameters \eqn{(\mu_{0c}, \kappa_{nt}, \nu_{nt}, r \Lambda_{nt})}
#' from the posterior parameters of the treatment group (Yamaguchi et al.,
#' submitted, Section A.2.6).
#'
#' @param design An object of class \code{"qdm_design_continuous2"}.
#' @param post_t Posterior parameters of the treatment group (see
#'   \code{posterior_continuous2()}).
#' @param n,mean,cov Current control data (concurrent control only).
#'
#' @return A list as returned by \code{posterior_continuous2()}.
#' @noRd
control_posterior_continuous2 <- function(design, post_t, n = 0L,
                                          mean = c(0, 0),
                                          cov = matrix(0, 2L, 2L)) {
  switch(
    design$control,
    concurrent = posterior_continuous2(design$prior_c, n, mean, cov,
                                       external = design$external_c),
    external = posterior_continuous2(design$prior_c, 0L,
                                     external = design$external_c),
    hypothetical = list(
      mu = matrix(design$hypothetical_c$mean, nrow(post_t$mu), 2L,
                  byrow = TRUE),
      kappa = post_t$kappa, nu = post_t$nu,
      Lambda = design$hypothetical_c$r * post_t$Lambda
    )
  )
}
