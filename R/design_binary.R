#' Design with one binary endpoint
#'
#' Specifies a proof-of-concept study with one binary endpoint: the sample
#' sizes, the prior distributions, the external data, the form of the control
#' group and the decision criterion. The design is used by
#' \code{\link{qdm_decide}}, \code{\link{qdm_oc}} and
#' \code{\link{qdm_calibrate}}.
#'
#' The response rates \eqn{\pi_t} and \eqn{\pi_c} of the treatment and control
#' groups have Beta prior distributions, and the treatment effect is
#' \eqn{\theta = \pi_t - \pi_c}. The control group takes one of three forms.
#' \describe{
#'   \item{Concurrent control}{\code{n_c} is at least 1. External control data
#'     can be added through \code{external_c}.}
#'   \item{External control only}{\code{n_c = 0} and \code{external_c} is
#'     given. The distribution of \eqn{\pi_c} is the posterior distribution
#'     given the external control data.}
#'   \item{Hypothetical control}{\code{n_c = 0} and \code{hypothetical_c} is
#'     given. The distribution of \eqn{\pi_c} is the posterior distribution
#'     given the hypothetical control data. The hypothetical control
#'     \eqn{\pi_c \sim \textrm{Beta}(\alpha_t + z, \beta_t + n_t - z)} of
#'     Kang et al. (2026, Section 4.1.4) is obtained with \code{prior_c} equal
#'     to \code{prior_t} and \code{hypothetical(data_binary(n_t, z))}.}
#' }
#' External treatment data (\code{external_t}) can be used with any form of
#' the control group. External data enter through a power prior: with
#' external data of \eqn{n_e} subjects and \eqn{y_e} responders and weight
#' \eqn{a_0}, the posterior distribution of the response rate is
#' \eqn{\textrm{Beta}(\alpha + a_0 y_e + y, \beta + a_0 (n_e - y_e) + n - y)}.
#'
#' @param n_t Number of subjects in the treatment group (a whole number of at
#'   least 1).
#' @param n_c Number of subjects in the concurrent control group (a whole
#'   number of at least 0).
#' @param prior_t,prior_c Beta prior distributions of the response rates of
#'   the treatment and control groups, created by \code{\link{prior_beta}}.
#' @param criterion Decision criterion for one endpoint, created by
#'   \code{\link{criterion_posterior}} or \code{\link{criterion_predictive}}.
#' @param external_t,external_c External data of the treatment and control
#'   groups, created by \code{\link{borrow}} with \code{\link{data_binary}},
#'   or \code{NULL} (no external data).
#' @param hypothetical_c Hypothetical control, created by
#'   \code{\link{hypothetical}} with \code{\link{data_binary}}, or
#'   \code{NULL}. It can be used only when \code{n_c = 0}.
#'
#' @return An object of class \code{c("qdm_design_binary", "qdm_design")}, a
#'   list with the arguments and the element \code{control}
#'   (\code{"concurrent"}, \code{"external"} or \code{"hypothetical"}).
#'
#' @references
#' Kang, R., Yamaguchi, Y. and Han, C. (2026). Bayesian quantitative
#' decision-making framework in cell and gene therapy development for rare
#' diseases. \emph{Journal of Biopharmaceutical Statistics}.
#' \doi{10.1080/10543406.2026.2655410}
#'
#' @family designs
#' @examples
#' # Concurrent control with external control data
#' design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8),
#'   external_c = borrow(data_binary(n = 20, y = 8), a0 = 0.5)
#' )
#'
#' # Uncontrolled study with a hypothetical control
#' design_binary(
#'   n_t = 12, n_c = 0,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_predictive(theta_null = 0, m_t = 50,
#'                                    gamma_go = 0.9, gamma_nogo = 0.7),
#'   hypothetical_c = hypothetical(data = data_binary(n = 12, y = 4))
#' )
#' @export
design_binary <- function(n_t, n_c, prior_t, prior_c, criterion,
                          external_t = NULL, external_c = NULL,
                          hypothetical_c = NULL) {
  check_numeric(n_t, "n_t", lower = 1, integer = TRUE)
  check_numeric(n_c, "n_c", lower = 0, integer = TRUE)
  check_prior_beta(prior_t, "prior_t")
  check_prior_beta(prior_c, "prior_c")
  if (!inherits(criterion, "qdm_criterion") || criterion$n_endpoints != 1L) {
    stop(paste("'criterion' must be created by criterion_posterior() or",
               "criterion_predictive() for one endpoint."), call. = FALSE)
  }
  check_external_binary(external_t, "external_t")
  check_external_binary(external_c, "external_c")
  if (!is.null(hypothetical_c) &&
      (!inherits(hypothetical_c, "qdm_hypothetical") ||
       !inherits(hypothetical_c$data, "qdm_data_binary"))) {
    stop(paste("'hypothetical_c' must be NULL or created by hypothetical()",
               "with data_binary()."), call. = FALSE)
  }
  control <- control_form(n_c, external_c, hypothetical_c)
  structure(list(endpoint = "binary", n_t = as.integer(round(n_t)),
                 n_c = as.integer(round(n_c)), prior_t = prior_t,
                 prior_c = prior_c, criterion = criterion,
                 external_t = external_t, external_c = external_c,
                 hypothetical_c = hypothetical_c, control = control),
            class = c("qdm_design_binary", "qdm_design"))
}
