#' Decision criterion based on posterior probability
#'
#' Specifies a Go/NoGo decision criterion based on the posterior probability
#' of the treatment effect \eqn{\theta} (the difference between the treatment
#' and control groups) given the proof-of-concept data \eqn{D}.
#'
#' For one endpoint, the Go region is \eqn{\theta > \theta_{TV}} and the NoGo
#' region is \eqn{\theta \le \theta_{MAV}}. The Go criterion is
#' \eqn{\Pr(\theta > \theta_{TV} \mid D) \ge \gamma_{go}} and the NoGo
#' criterion is \eqn{\Pr(\theta \le \theta_{MAV} \mid D) \ge \gamma_{nogo}}.
#' In the notation of Kang et al. (2026), where the NoGo criterion is
#' \eqn{\Pr(\theta > \theta_{MAV} \mid D) \le \gamma_2}, this corresponds to
#' \eqn{\gamma_{nogo} = 1 - \gamma_2}.
#'
#' For two endpoints, the target values and minimum acceptable values divide
#' the plane of \eqn{(\theta_1, \theta_2)} into nine regions (see
#' \code{\link{qdm_regions}}). The Go criterion is that the posterior
#' probability of the union of \code{go_regions} is at least
#' \eqn{\gamma_{go}}, and the NoGo criterion is that the posterior probability
#' of the union of \code{nogo_regions} is at least \eqn{\gamma_{nogo}}
#' (Yamaguchi et al., submitted).
#'
#' The decision is Go if only the Go criterion is met, NoGo if only the NoGo
#' criterion is met, Gray if neither is met, and Miss if both are met. For one
#' endpoint, Miss can occur only when \eqn{\gamma_{go} + \gamma_{nogo} \le 1}.
#'
#' @param theta_tv Target value (TV) of the treatment effect: a number for one
#'   endpoint or a numeric vector of length 2 for two endpoints.
#' @param theta_mav Minimum acceptable value (MAV) of the treatment effect, of
#'   the same length as \code{theta_tv}. Each element must be smaller than the
#'   corresponding element of \code{theta_tv}.
#' @param gamma_go Probability threshold for the Go criterion, strictly
#'   between 0 and 1.
#' @param gamma_nogo Probability threshold for the NoGo criterion, strictly
#'   between 0 and 1.
#' @param go_regions Region numbers (from 1 to 9) forming the Go region for two
#'   endpoints. Ignored for one endpoint. The default, region 1, requires both
#'   effects to exceed their target values.
#' @param nogo_regions Region numbers (from 1 to 9) forming the NoGo region for
#'   two endpoints. Ignored for one endpoint. The default, region 9, requires
#'   both effects to be at most their minimum acceptable values.
#'
#' @return An object of class
#'   \code{c("qdm_criterion_posterior", "qdm_criterion")}.
#'
#' @references
#' Kang, R., Yamaguchi, Y. and Han, C. (2026). Bayesian quantitative
#' decision-making framework in cell and gene therapy development for rare
#' diseases. \emph{Journal of Biopharmaceutical Statistics}.
#' \doi{10.1080/10543406.2026.2655410}
#'
#' @family criteria
#' @examples
#' # One endpoint
#' criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                     gamma_go = 0.5, gamma_nogo = 0.8)
#'
#' # Two endpoints
#' criterion_posterior(theta_tv = c(0.2, 0.15), theta_mav = c(0.1, 0.05),
#'                     gamma_go = 0.3, gamma_nogo = 0.4)
#' @export
criterion_posterior <- function(theta_tv, theta_mav, gamma_go, gamma_nogo,
                                go_regions = 1L, nogo_regions = 9L) {
  check_numeric(theta_tv, "theta_tv", len = NULL)
  n_endpoints <- length(theta_tv)
  if (!(n_endpoints %in% c(1L, 2L))) {
    stop("'theta_tv' must have length 1 (one endpoint) or 2 (two endpoints).",
         call. = FALSE)
  }
  check_numeric(theta_mav, "theta_mav", len = n_endpoints)
  if (any(theta_tv <= theta_mav)) {
    stop("Each element of 'theta_tv' must be greater than that of 'theta_mav'.",
         call. = FALSE)
  }
  check_numeric(gamma_go, "gamma_go", lower = 0, upper = 1,
                lower_open = TRUE, upper_open = TRUE)
  check_numeric(gamma_nogo, "gamma_nogo", lower = 0, upper = 1,
                lower_open = TRUE, upper_open = TRUE)
  if (n_endpoints == 2L) {
    regions <- check_regions(go_regions, nogo_regions, n_regions = 9L)
  } else {
    regions <- list(go = NULL, nogo = NULL)
  }
  structure(list(type = "posterior", n_endpoints = n_endpoints,
                 theta_tv = theta_tv, theta_mav = theta_mav,
                 gamma_go = gamma_go, gamma_nogo = gamma_nogo,
                 go_regions = regions$go, nogo_regions = regions$nogo),
            class = c("qdm_criterion_posterior", "qdm_criterion"))
}
