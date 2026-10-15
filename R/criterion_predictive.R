#' Decision criterion based on posterior predictive probability
#'
#' Specifies a Go/NoGo decision criterion based on the posterior predictive
#' probability of the treatment effect \eqn{\tilde\theta} that will be
#' observed in a future study with \code{m_t} and \code{m_c} subjects in the
#' treatment and control groups, given the proof-of-concept data \eqn{D}.
#' The future effect is the difference between the future sample proportions
#' (binary endpoints) or sample means (continuous endpoints).
#'
#' For one endpoint, the Go region is \eqn{\tilde\theta > \theta_{NULL}} and
#' the NoGo region is \eqn{\tilde\theta \le \theta_{NULL}}. The Go criterion
#' is \eqn{\Pr(\tilde\theta > \theta_{NULL} \mid D) \ge \gamma_{go}} and the
#' NoGo criterion is
#' \eqn{\Pr(\tilde\theta \le \theta_{NULL} \mid D) \ge \gamma_{nogo}}. Because
#' the two probabilities add up to 1, a Gray decision is possible only when
#' \eqn{\gamma_{go} + \gamma_{nogo} > 1}, and Miss decisions can occur only when
#' \eqn{\gamma_{go} + \gamma_{nogo} \le 1}.
#'
#' For two endpoints, the null values divide the plane of
#' \eqn{(\tilde\theta_1, \tilde\theta_2)} into four regions (see
#' \code{\link{qdm_regions}}), and the criteria are based on the predictive
#' probabilities of the unions of \code{go_regions} and \code{nogo_regions}
#' (Yamaguchi et al., submitted, equations 7 and 8). With \code{rule = "and"}
#' or \code{"or"}, they are based on the marginal predictive probabilities of
#' each endpoint, \eqn{\Pr(\tilde\theta_k > \theta_{NULLk} \mid D)} for Go
#' and \eqn{\Pr(\tilde\theta_k \le \theta_{NULLk} \mid D)} for NoGo,
#' combined with AND or OR (equations 9 and 10).
#'
#' @param theta_null Null value of the treatment effect in the future study: a
#'   number for one endpoint or a numeric vector of length 2 for two
#'   endpoints.
#' @param m_t,m_c Numbers of subjects in the treatment and control groups of
#'   the future study (whole numbers of at least 1).
#' @param gamma_go Probability threshold for the Go criterion, strictly
#'   between 0 and 1. With \code{rule = "and"} or \code{"or"}, one threshold
#'   per endpoint (a single number is used for both).
#' @param gamma_nogo Probability threshold for the NoGo criterion, strictly
#'   between 0 and 1, given in the same way as \code{gamma_go}.
#' @param go_regions Region numbers (from 1 to 4) forming the Go region for two
#'   endpoints with \code{rule = "joint"}. Ignored for one endpoint. The
#'   default, region 1, requires both future effects to exceed their null
#'   values.
#' @param nogo_regions Region numbers (from 1 to 4) forming the NoGo region for
#'   two endpoints with \code{rule = "joint"}. Ignored for one endpoint. The
#'   default, region 4, requires both future effects to be at most their null
#'   values.
#' @param rule For two endpoints, \code{"joint"} (the default; probabilities
#'   of unions of regions), \code{"and"} or \code{"or"} (marginal
#'   probabilities of each endpoint, combined with AND or OR). Ignored for one
#'   endpoint.
#'
#' @return An object of class
#'   \code{c("qdm_criterion_predictive", "qdm_criterion")}.
#'
#' @references
#' Yamaguchi, Y., Hu, H., Huang, Y., Han, C., Homma, G. and Wang, A.
#' (submitted). Bivariate Bayesian quantitative decision-making for two
#' endpoints in rare disease drug development.
#'
#' @family criteria
#' @examples
#' criterion_predictive(theta_null = 0, m_t = 50, gamma_go = 0.9,
#'                      gamma_nogo = 0.7)
#' @export
criterion_predictive <- function(theta_null, m_t, m_c = m_t, gamma_go,
                                 gamma_nogo, go_regions = 1L,
                                 nogo_regions = 4L,
                                 rule = c("joint", "and", "or")) {
  check_numeric(theta_null, "theta_null", len = NULL)
  n_endpoints <- length(theta_null)
  if (!(n_endpoints %in% c(1L, 2L))) {
    stop("'theta_null' must have length 1 (one endpoint) or 2 (two endpoints).",
         call. = FALSE)
  }
  check_numeric(m_t, "m_t", lower = 1, integer = TRUE)
  check_numeric(m_c, "m_c", lower = 1, integer = TRUE)
  rule <- match.arg(rule)
  given <- !missing(go_regions) || !missing(nogo_regions)
  r <- check_rule(rule, n_endpoints, gamma_go, gamma_nogo, go_regions,
                  nogo_regions, regions_given = given, n_regions = 4L)
  structure(list(type = "predictive", n_endpoints = n_endpoints,
                 theta_null = theta_null, m_t = as.integer(round(m_t)),
                 m_c = as.integer(round(m_c)), gamma_go = r$gamma_go,
                 gamma_nogo = r$gamma_nogo, go_regions = r$go,
                 nogo_regions = r$nogo, rule = r$rule),
            class = c("qdm_criterion_predictive", "qdm_criterion"))
}
