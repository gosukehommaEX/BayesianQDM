#' Design with one time-to-event endpoint
#'
#' Specifies a proof-of-concept study with one time-to-event endpoint
#' evaluated by the milestone survival probability: the sample sizes, the
#' prior distributions, the external data, the form of the control group,
#' the decision criterion, the milestone time and the conduct of the study.
#' The design is used by \code{\link{qdm_decide}}, \code{\link{qdm_oc}},
#' \code{\link{qdm_calibrate}} and \code{\link{qdm_boundary}}.
#'
#' The event times of group \eqn{k} (\eqn{k = t, c}) are exponential with
#' hazard \eqn{\lambda_k}, which has a Gamma prior distribution
#' (\code{\link{prior_gamma}}). The data of a group are summarized by the
#' number of events \eqn{d} and the total observed time \eqn{E}
#' (\code{\link{data_tte}}), and the posterior distribution of the hazard is
#' Gamma with shape \eqn{a + d} and rate \eqn{b + E}. External data enter
#' through a power prior with weight \eqn{a_0}, which adds \eqn{a_0 d_e} to the
#' shape and \eqn{a_0 E_e} to the rate. The treatment effect is the
#' difference of the milestone survival probabilities at time
#' \eqn{t_0 = } \code{milestone},
#' \eqn{\theta = S_t(t_0) - S_c(t_0)} with \eqn{S_k(t_0) = \exp(-\lambda_k t_0)}.
#'
#' The control group takes one of four forms.
#' \describe{
#'   \item{Concurrent control}{\code{n_c} is at least 1. External control
#'     data can be added through \code{external_c}.}
#'   \item{External control only}{\code{n_c = 0} and \code{external_c} is
#'     given. The distribution of \eqn{\lambda_c} is the posterior
#'     distribution given the external control data.}
#'   \item{Hypothetical control}{\code{n_c = 0} and \code{hypothetical_c} is
#'     given by \code{hypothetical(data_tte(events, exposure))}, whose events
#'     and total observed time are added to \code{prior_c}.}
#'   \item{Fixed control}{\code{n_c = 0} and \code{fixed_c} is given: the
#'     control milestone survival probability is a known value, as in a
#'     single-arm study compared with a benchmark (\code{fixed_c = 0} makes
#'     the effect the treatment milestone survival itself). \code{prior_c}
#'     is not used and must be \code{NULL}.}
#' }
#'
#' Posterior probabilities of the Go and NoGo regions are computed exactly
#' for a fixed control (from the Gamma distribution of \eqn{\lambda_t}) and by
#' numerical integration otherwise (see \code{\link{psurvdiff}}).
#'
#' For the posterior predictive criterion (\code{\link{criterion_predictive}}),
#' a future study with \code{m_t} treatment and \code{m_c} control subjects
#' (no control group for a fixed control) is conducted as specified by
#' \code{conduct_future}. The milestone survival of each group in the future
#' study is estimated by \eqn{\exp(-t_0 (a + d) / (b + E))}, the value at the
#' posterior mean of the hazard given the future data alone with the prior
#' \eqn{(a, b)} of the group (\code{prior_t} or \code{prior_c}); with a weak
#' prior this is close to the maximum likelihood estimate and remains
#' defined when the future study has no events. The future effect is the
#' difference of the two estimates (or the treatment estimate minus the
#' fixed value). The predictive probabilities are estimated by simulating
#' \code{nMC} future studies, each with hazards drawn from the posterior
#' distributions given the proof-of-concept data (see
#' \code{\link{conduct_tte}} for the simulation of a study). They are
#' therefore Monte Carlo estimates, and their computing time is proportional
#' to \code{nMC} times the number of future subjects.
#'
#' The operating characteristics (\code{\link{qdm_oc}},
#' \code{\link{qdm_calibrate}}) are computed by simulating \code{nsim}
#' studies conducted as specified by \code{conduct}, with exponential event
#' times whose hazards are given by the median survival times of the
#' scenarios. With a predictive criterion, every simulated study requires
#' \code{nMC} simulated future studies, so small values of \code{nsim} and
#' \code{nMC} are advisable for exploration.
#'
#' @param n_t Number of subjects in the treatment group (a whole number of at
#'   least 1).
#' @param n_c Number of subjects in the concurrent control group (a whole
#'   number of at least 0).
#' @param prior_t,prior_c Gamma prior distributions of the hazards of the
#'   treatment and control groups, created by \code{\link{prior_gamma}}.
#'   \code{prior_c} must be \code{NULL} (the default) for a fixed control and
#'   is required otherwise.
#' @param criterion Decision criterion for one endpoint, created by
#'   \code{\link{criterion_posterior}} or \code{\link{criterion_predictive}},
#'   with thresholds on the scale of the milestone survival difference (or
#'   of the milestone survival itself with \code{fixed_c = 0}).
#' @param milestone Milestone time \eqn{t_0} (positive), in the same unit as
#'   the total observed time of the data.
#' @param conduct Conduct of the study, created by \code{\link{conduct_tte}},
#'   used to simulate studies for the operating characteristics; \code{NULL}
#'   (the default) if they are not needed. With an analysis at an event
#'   (\code{cut_events}), the number of events must not exceed the number of
#'   subjects (\code{n_t} plus \code{n_c}).
#' @param external_t,external_c External data of the treatment and control
#'   groups, created by \code{\link{borrow}} with \code{\link{data_tte}}, or
#'   \code{NULL} (no external data).
#' @param hypothetical_c Hypothetical control, created by
#'   \code{\link{hypothetical}} with \code{\link{data_tte}}, or \code{NULL}.
#'   It can be used only when \code{n_c = 0}.
#' @param fixed_c Fixed control milestone survival probability, a number from
#'   0 to 1, or \code{NULL}. It can be used only when \code{n_c = 0}.
#' @param conduct_future Conduct of the future study for a predictive
#'   criterion, created by \code{\link{conduct_tte}}; by default the same as
#'   \code{conduct}. Required for a predictive criterion. With an analysis
#'   at an event, the number of events must not exceed the number of future
#'   subjects.
#' @param nMC Number of simulated future studies for posterior predictive
#'   probabilities.
#'
#' @return An object of class \code{c("qdm_design_tte", "qdm_design")}, a
#'   list with the arguments and the element \code{control}
#'   (\code{"concurrent"}, \code{"external"}, \code{"hypothetical"} or
#'   \code{"fixed"}).
#'
#' @family designs
#' @examples
#' # Single-arm study compared with a fixed control 3-month survival of 0.72
#' des <- design_tte(
#'   n_t = 30, n_c = 0, prior_t = prior_gamma(0.1, 1),
#'   criterion = criterion_posterior(theta_tv = 0.15, theta_mav = 0,
#'                                   gamma_go = 0.6, gamma_nogo = 0.8),
#'   milestone = 3,
#'   conduct = conduct_tte(accrual_time = c(0, 5, 10, 15),
#'                         accrual_rate = c(0.4, 2, 4.5),
#'                         dropout_hazard = 0.00226, cut_time = 15.44),
#'   fixed_c = 0.72
#' )
#' des
#' qdm_decide(des, data_t = data_tte(events = 6, exposure = 150))
#'
#' # Concurrent control
#' design_tte(
#'   n_t = 30, n_c = 30, prior_t = prior_gamma(0.1, 1),
#'   prior_c = prior_gamma(0.1, 1),
#'   criterion = criterion_posterior(theta_tv = 0.15, theta_mav = 0,
#'                                   gamma_go = 0.6, gamma_nogo = 0.8),
#'   milestone = 3,
#'   conduct = conduct_tte(accrual_time = c(0, 12), cut_events = 25)
#' )
#' @export
design_tte <- function(n_t, n_c, prior_t, prior_c = NULL, criterion,
                       milestone, conduct = NULL, external_t = NULL,
                       external_c = NULL, hypothetical_c = NULL,
                       fixed_c = NULL, conduct_future = conduct,
                       nMC = 10000L) {
  check_numeric(n_t, "n_t", lower = 1, integer = TRUE)
  check_numeric(n_c, "n_c", lower = 0, integer = TRUE)
  check_prior_gamma(prior_t, "prior_t")
  if (!inherits(criterion, "qdm_criterion") || criterion$n_endpoints != 1L) {
    stop(paste("'criterion' must be created by criterion_posterior() or",
               "criterion_predictive() for one endpoint."), call. = FALSE)
  }
  check_numeric(milestone, "milestone", lower = 0, lower_open = TRUE)
  check_external_tte(external_t, "external_t")
  check_external_tte(external_c, "external_c")
  if (!is.null(hypothetical_c) &&
      (!inherits(hypothetical_c, "qdm_hypothetical") ||
       !inherits(hypothetical_c$data, "qdm_data_tte"))) {
    stop(paste("'hypothetical_c' must be NULL or created by hypothetical()",
               "with data_tte()."), call. = FALSE)
  }
  fixed_c <- check_fixed_c(fixed_c, 1L, probability = TRUE)
  control <- control_form(n_c, external_c, hypothetical_c, fixed_c)
  check_prior_control(prior_c, control, check_prior_gamma)
  n_t <- as.integer(round(n_t))
  n_c <- as.integer(round(n_c))
  check_conduct(conduct, "conduct", n_t + n_c, allow_null = TRUE)
  if (criterion$type == "predictive") {
    if (is.null(conduct_future)) {
      stop(paste("A predictive criterion needs 'conduct_future' (or",
                 "'conduct'), the conduct of the future study."),
           call. = FALSE)
    }
    n_future <- criterion$m_t + if (control == "fixed") 0L else criterion$m_c
    check_conduct(conduct_future, "conduct_future", n_future)
  } else {
    conduct_future <- NULL
  }
  check_numeric(nMC, "nMC", lower = 1, upper = .Machine$integer.max,
                integer = TRUE)
  structure(list(endpoint = "tte", n_t = n_t, n_c = n_c, prior_t = prior_t,
                 prior_c = prior_c, criterion = criterion,
                 milestone = milestone, conduct = conduct,
                 external_t = external_t, external_c = external_c,
                 hypothetical_c = hypothetical_c, fixed_c = fixed_c,
                 control = control, conduct_future = conduct_future,
                 nMC = as.integer(round(nMC))),
            class = c("qdm_design_tte", "qdm_design"))
}
