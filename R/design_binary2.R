#' Design with two binary endpoints
#'
#' Specifies a proof-of-concept study with two binary endpoints: the sample
#' sizes, the Dirichlet prior distributions, the external data, the form of
#' the control group, the decision criterion for two endpoints and the
#' computation settings. The design is used by \code{\link{qdm_decide}},
#' \code{\link{qdm_oc}} and \code{\link{qdm_calibrate}}.
#'
#' The responses of a subject on the two endpoints fall into one of the four
#' cells \code{"00"}, \code{"01"}, \code{"10"} and \code{"11"} (first digit:
#' endpoint 1), whose probabilities have Dirichlet prior distributions in each
#' group (Yamaguchi et al., submitted, Appendix A.1). The treatment effects
#' are the differences of the marginal response rates,
#' \eqn{\theta_k = \pi_{tk} - \pi_{ck}} with \eqn{\pi_1 = p_{10} + p_{11}} and
#' \eqn{\pi_2 = p_{01} + p_{11}}. External data with counts \eqn{x_e} enter
#' through a power prior with weight \eqn{a_0}, which adds \eqn{a_0 x_e} to the
#' Dirichlet parameters. The control group takes one of four forms, as in
#' \code{\link{design_binary}}: a concurrent control group (\code{n_c} of at
#' least 1), external control data only (\code{n_c = 0} and
#' \code{external_c}), a hypothetical control (\code{n_c = 0} and
#' \code{hypothetical_c} created from \code{\link{data_binary2}}, whose counts
#' are added to \code{prior_c}; the hypothetical control distribution
#' Dir(beta) of Yamaguchi et al., Appendix A.1.6, is obtained with
#' \code{prior_c = prior_dirichlet(beta - counts)}), or a fixed control
#' (\code{n_c = 0} and \code{fixed_c}, the known response rates of the two
#' endpoints in the control group; \code{prior_c} must then be \code{NULL}).
#' With a fixed control, the posterior probabilities are estimated from draws
#' of the treatment distribution only, and the future effects of the
#' predictive criterion are \eqn{S_{tk} / m_t - \pi_{ck}} (a future effect
#' that differs from the null value by less than \eqn{10^{-10}} is treated as
#' equal to it; \code{m_c} is not used).
#'
#' The posterior probabilities of the regions (posterior criteria) are
#' estimated by Monte Carlo simulation with \code{nMC} draws from the
#' Dirichlet posterior distributions. The posterior predictive probabilities
#' (predictive criteria) are computed exactly: the future effects depend only
#' on the numbers of responders on each endpoint, whose joint distribution in
#' each group is obtained from the Dirichlet-multinomial distribution, and a
#' future effect equal to a null value is not counted as larger.
#'
#' For \code{\link{qdm_oc}} and \code{\link{qdm_calibrate}}, with
#' \code{enumerate = TRUE} all outcomes of the study are enumerated: there are
#' \code{choose(n + 3, 3)} outcomes in a group of \code{n} subjects (286 for
#' \code{n = 10}), and the probabilities of the regions are computed for
#' every pair of treatment and control outcomes, with the same Monte Carlo
#' draws for all pairs. For a posterior criterion the computing time is
#' proportional to the number of pairs times \code{nMC}, and enumeration is
#' limited to 1e7 pairs (for example, 24 subjects in each group). With
#' \code{enumerate = FALSE}, \code{nsim} studies are simulated under each
#' scenario, which is faster for large sample sizes.
#'
#' @param n_t Number of subjects in the treatment group (a whole number of at
#'   least 1).
#' @param n_c Number of subjects in the concurrent control group (a whole
#'   number of at least 0).
#' @param prior_t,prior_c Dirichlet prior distributions of the four cell
#'   probabilities of the treatment and control groups, created by
#'   \code{\link{prior_dirichlet}}. \code{prior_c} must be \code{NULL} (the
#'   default) for a fixed control and is required otherwise.
#' @param criterion Decision criterion for two endpoints, created by
#'   \code{\link{criterion_posterior}} or \code{\link{criterion_predictive}}.
#'   Its argument \code{direction} sets whether larger or smaller effects
#'   are favorable on each endpoint.
#' @param external_t,external_c External data of the treatment and control
#'   groups, created by \code{\link{borrow}} with
#'   \code{\link{data_binary2}}, or \code{NULL} (no external data).
#' @param hypothetical_c Hypothetical control, created by
#'   \code{\link{hypothetical}} with \code{\link{data_binary2}}, or
#'   \code{NULL}. It can be used only when \code{n_c = 0}.
#' @param fixed_c Fixed control response rates of endpoints 1 and 2, a numeric
#'   vector of length 2 with values from 0 to 1, or \code{NULL}. It can be
#'   used only when \code{n_c = 0}.
#' @param nMC Number of Monte Carlo draws for posterior probabilities.
#' @param enumerate Logical; if \code{TRUE} (the default), the operating
#'   characteristics are computed by enumerating all outcomes, otherwise by
#'   simulation.
#'
#' @return An object of class \code{c("qdm_design_binary2", "qdm_design")},
#'   a list with the arguments and the element \code{control}
#'   (\code{"concurrent"}, \code{"external"}, \code{"hypothetical"} or
#'   \code{"fixed"}).
#'
#' @references
#' Yamaguchi, Y., Hu, H., Huang, Y., Han, C., Homma, G. and Wang, A.
#' (submitted). Bivariate Bayesian quantitative decision-making for two
#' endpoints in rare disease drug development.
#'
#' @family designs
#' @examples
#' # Case study 1 of Yamaguchi et al. (submitted), with a small nMC for speed
#' des <- design_binary2(
#'   n_t = 10, n_c = 10,
#'   prior_t = prior_dirichlet(rep(0.25, 4)),
#'   prior_c = prior_dirichlet(rep(0.25, 4)),
#'   criterion = criterion_posterior(theta_tv = c(0.15, 0.15),
#'                                   theta_mav = c(0.10, 0.10),
#'                                   gamma_go = 0.263, gamma_nogo = 0.347),
#'   nMC = 2000
#' )
#' des
#' qdm_decide(des, data_t = data_binary2(c(3, 2, 2, 3)),
#'            data_c = data_binary2(c(6, 2, 1, 1)), seed = 1)
#' @export
design_binary2 <- function(n_t, n_c, prior_t, prior_c = NULL, criterion,
                           external_t = NULL, external_c = NULL,
                           hypothetical_c = NULL, fixed_c = NULL,
                           nMC = 10000L, enumerate = TRUE) {
  check_numeric(n_t, "n_t", lower = 1, integer = TRUE)
  check_numeric(n_c, "n_c", lower = 0, integer = TRUE)
  check_prior_dirichlet(prior_t, "prior_t")
  if (!inherits(criterion, "qdm_criterion") || criterion$n_endpoints != 2L) {
    stop(paste("'criterion' must be created by criterion_posterior() or",
               "criterion_predictive() for two endpoints."), call. = FALSE)
  }
  check_external_binary2(external_t, "external_t")
  check_external_binary2(external_c, "external_c")
  if (!is.null(hypothetical_c) &&
      (!inherits(hypothetical_c, "qdm_hypothetical") ||
       !inherits(hypothetical_c$data, "qdm_data_binary2"))) {
    stop(paste("'hypothetical_c' must be NULL or created by hypothetical()",
               "with data_binary2()."), call. = FALSE)
  }
  fixed_c <- check_fixed_c(fixed_c, 2L, probability = TRUE)
  control <- control_form(n_c, external_c, hypothetical_c, fixed_c)
  check_prior_control(prior_c, control, check_prior_dirichlet)
  check_numeric(nMC, "nMC", lower = 1, upper = .Machine$integer.max,
                integer = TRUE)
  check_flag(enumerate, "enumerate")
  structure(list(endpoint = "binary2", n_t = as.integer(round(n_t)),
                 n_c = as.integer(round(n_c)), prior_t = prior_t,
                 prior_c = prior_c, criterion = criterion,
                 external_t = external_t, external_c = external_c,
                 hypothetical_c = hypothetical_c, fixed_c = fixed_c,
                 control = control, nMC = as.integer(round(nMC)),
                 enumerate = enumerate),
            class = c("qdm_design_binary2", "qdm_design"))
}
