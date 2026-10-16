#' Design with one time-to-event endpoint
#'
#' Specifies a proof-of-concept study with one time-to-event endpoint
#' evaluated by the milestone survival probability or by the hazard ratio:
#' the sample sizes, the prior distributions, the external data, the form of
#' the control group, the decision criterion, the effect measure and the
#' conduct of the study.
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
#' shape and \eqn{a_0 E_e} to the rate. The treatment effect is either
#' \describe{
#'   \item{\code{effect = "milestone"}}{the difference of the milestone
#'     survival probabilities at time \eqn{t_0 = } \code{milestone},
#'     \eqn{\theta = S_t(t_0) - S_c(t_0)} with
#'     \eqn{S_k(t_0) = \exp(-\lambda_k t_0)}, or}
#'   \item{\code{effect = "hazard_ratio"}}{the hazard ratio
#'     \eqn{\theta = \lambda_t / \lambda_c}, as in Huang et al. (2019), whose
#'     decision rules are based on the normal approximation of the estimated
#'     log hazard ratio. Smaller hazard ratios are favorable, so the
#'     criterion is usually specified with \code{direction = "less"}, for
#'     example \code{theta_tv = 0.67} and \code{theta_mav = 0.8}. Under the
#'     exponential model the hazard ratio is also the ratio of the median
#'     survival times of the control and treatment groups.}
#' }
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
#'     control milestone survival probability (\code{effect = "milestone"})
#'     or the control median survival time (\code{effect = "hazard_ratio"},
#'     control hazard \eqn{\log(2)} / \code{fixed_c}) is a known value, as
#'     in a single-arm study compared with a benchmark (with the milestone
#'     survival, \code{fixed_c = 0} makes the effect the treatment milestone
#'     survival itself). \code{prior_c} is not used and must be
#'     \code{NULL}.}
#' }
#'
#' Posterior probabilities of the Go and NoGo regions are computed exactly
#' for a fixed control (from the Gamma distribution of \eqn{\lambda_t}). With
#' another control, they are computed by numerical integration for the
#' milestone survival (see \code{\link{psurvdiff}}) and exactly from the F
#' distribution for the hazard ratio (see \code{\link{phazardratio}}).
#'
#' For the posterior predictive criterion (\code{\link{criterion_predictive}}),
#' a future study with \code{m_t} treatment and \code{m_c} control subjects
#' (no control group for a fixed control) is conducted as specified by
#' \code{conduct_future}. The hazard of each group in the future study is
#' estimated by \eqn{(a + d) / (b + E)}, the posterior mean of the hazard
#' given the future data alone with the prior \eqn{(a, b)} of the group
#' (\code{prior_t} or \code{prior_c}); with a weak prior this is close to the
#' maximum likelihood estimate \eqn{d / E} and, unlike it, is positive when a
#' group has no events, so that the hazard ratio is defined. The milestone
#' survival is estimated by \eqn{\exp(-t_0 (a + d) / (b + E))}. The future
#' effect is the difference
#' of the two milestone survival estimates (or the treatment estimate minus
#' the fixed value), or the ratio of the two hazard estimates (or the
#' treatment estimate divided by the fixed control hazard). The predictive
#' probabilities are estimated by simulating \code{nMC} future studies,
#' each with hazards drawn from the posterior
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
#'   with thresholds on the scale of the effect: the milestone survival
#'   difference (or the milestone survival itself with \code{fixed_c = 0}),
#'   or the hazard ratio (positive thresholds).
#' @param milestone Milestone time \eqn{t_0} (positive), in the same unit as
#'   the total observed time of the data. Required for
#'   \code{effect = "milestone"}; must be \code{NULL} (the default) for the
#'   hazard ratio.
#' @param effect Effect measure: \code{"milestone"} (the default; difference
#'   of the milestone survival probabilities) or \code{"hazard_ratio"}.
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
#'   0 to 1 (\code{effect = "milestone"}), or fixed control median survival
#'   time, a positive number (\code{effect = "hazard_ratio"}); or \code{NULL}.
#'   It can be used only when \code{n_c = 0}.
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
#' @references
#' Huang, B., Talukder, E., Han, L. and Kuan, P.-F. (2019). Quantitative
#' decision-making in randomized Phase II studies with a time-to-event
#' endpoint. \emph{Journal of Biopharmaceutical Statistics}, 29(1), 189--202.
#' \doi{10.1080/10543406.2018.1489400}
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
#'
#' # Hazard ratio in the randomized study of the example of Huang et al.
#' # (2019) (2:1 randomization, analysis at 127 events): Go if HR < 0.67
#' # with probability at least 0.5, NoGo if HR >= 0.8 with probability at
#' # least 0.6
#' design_tte(
#'   n_t = 112, n_c = 56, prior_t = prior_gamma(0.01, 0.01),
#'   prior_c = prior_gamma(0.01, 0.01),
#'   criterion = criterion_posterior(theta_tv = 0.67, theta_mav = 0.8,
#'                                   gamma_go = 0.5, gamma_nogo = 0.6,
#'                                   direction = "less"),
#'   effect = "hazard_ratio",
#'   conduct = conduct_tte(accrual_time = c(0, 13), dropout_hazard = 0.003,
#'                         cut_events = 127)
#' )
#' @export
design_tte <- function(n_t, n_c, prior_t, prior_c = NULL, criterion,
                       milestone = NULL,
                       effect = c("milestone", "hazard_ratio"),
                       conduct = NULL, external_t = NULL, external_c = NULL,
                       hypothetical_c = NULL, fixed_c = NULL,
                       conduct_future = conduct, nMC = 10000L) {
  check_numeric(n_t, "n_t", lower = 1, integer = TRUE)
  check_numeric(n_c, "n_c", lower = 0, integer = TRUE)
  check_prior_gamma(prior_t, "prior_t")
  if (!inherits(criterion, "qdm_criterion") || criterion$n_endpoints != 1L) {
    stop(paste("'criterion' must be created by criterion_posterior() or",
               "criterion_predictive() for one endpoint."), call. = FALSE)
  }
  effect <- match.arg(effect)
  if (effect == "milestone") {
    if (is.null(milestone)) {
      stop("'milestone' must be given for effect = \"milestone\".",
           call. = FALSE)
    }
    check_numeric(milestone, "milestone", lower = 0, lower_open = TRUE)
  } else {
    if (!is.null(milestone)) {
      stop("'milestone' is not used for effect = \"hazard_ratio\".",
           call. = FALSE)
    }
    thresholds <- if (criterion$type == "posterior") {
      c(criterion$theta_tv, criterion$theta_mav)
    } else {
      criterion$theta_null
    }
    if (any(thresholds <= 0)) {
      stop(paste("The thresholds of the criterion must be positive hazard",
                 "ratios for effect = \"hazard_ratio\"."), call. = FALSE)
    }
  }
  check_external_tte(external_t, "external_t")
  check_external_tte(external_c, "external_c")
  if (!is.null(hypothetical_c) &&
      (!inherits(hypothetical_c, "qdm_hypothetical") ||
       !inherits(hypothetical_c$data, "qdm_data_tte"))) {
    stop(paste("'hypothetical_c' must be NULL or created by hypothetical()",
               "with data_tte()."), call. = FALSE)
  }
  if (effect == "hazard_ratio" && !is.null(fixed_c)) {
    # Fixed control median survival time
    check_numeric(fixed_c, "fixed_c", lower = 0, lower_open = TRUE)
    fixed_c <- as.numeric(fixed_c)
  } else {
    fixed_c <- check_fixed_c(fixed_c, 1L, probability = TRUE)
  }
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
                 milestone = milestone, effect = effect, conduct = conduct,
                 external_t = external_t, external_c = external_c,
                 hypothetical_c = hypothetical_c, fixed_c = fixed_c,
                 control = control, conduct_future = conduct_future,
                 nMC = as.integer(round(nMC))),
            class = c("qdm_design_tte", "qdm_design"))
}
