#' Design with one continuous endpoint
#'
#' Specifies a proof-of-concept study with one normally distributed endpoint:
#' the sample sizes, the prior distributions, the external data, the form of
#' the control group, the decision criterion and the method used to compute
#' the probabilities. The design is used by \code{\link{qdm_decide}},
#' \code{\link{qdm_oc}} and \code{\link{qdm_calibrate}}.
#'
#' The outcomes of group \eqn{k} (\eqn{k = t, c}) follow
#' \eqn{N(\mu_k, \sigma_k^2)}, and the treatment effect is
#' \eqn{\theta = \mu_t - \mu_c}. With the Normal-Inverse-chi-squared prior
#' (\code{\link{prior_nix}}), the posterior distribution of
#' \eqn{(\mu_k, \sigma_k^2)} is Normal-Inverse-chi-squared with parameters
#' \eqn{(\mu_{nk}, \kappa_{nk}, \nu_{nk}, \sigma_{nk}^2)}, and the posterior
#' distribution of the mean is
#' \eqn{t_{\nu_{nk}}(\mu_{nk}, \sigma_{nk}^2 / \kappa_{nk})}. With the vague
#' prior (\code{\link{prior_vague}}) and \eqn{n_k} subjects with sample mean
#' \eqn{\bar y_k} and standard deviation \eqn{s_k}, it is
#' \eqn{t_{n_k - 1}(\bar y_k, s_k^2 / n_k)}.
#'
#' External data enter through a power prior with weight \eqn{a_0}: data with
#' \eqn{n_e} subjects, mean \eqn{\bar y_e} and standard deviation \eqn{s_e}
#' are combined with the prior as if they were \eqn{a_0 n_e} subjects with the
#' same mean and the sum of squares \eqn{a_0 (n_e - 1) s_e^2}. The posterior
#' distribution then has a closed form (Huang et al., submitted). With the
#' vague prior and external data, the posterior distribution of the mean has
#' \eqn{a_0 n_e + n_k - 1} degrees of freedom.
#'
#' The control group takes one of four forms.
#' \describe{
#'   \item{Concurrent control}{\code{n_c} is at least 2. External control
#'     data can be added through \code{external_c}.}
#'   \item{External control only}{\code{n_c = 0} and \code{external_c} is
#'     given. The distribution of \eqn{\mu_c} is the posterior distribution
#'     given the external control data. With the vague prior, \eqn{a_0 n_e}
#'     must be greater than 1.}
#'   \item{Hypothetical control}{\code{n_c = 0} and \code{hypothetical_c} is
#'     given by \code{hypothetical(mean = mu0c, r = r)}. Following Kang et al.
#'     (2026, Section 4.2.4), the control mean has the same t distribution as
#'     the posterior distribution of the treatment mean, with location
#'     \eqn{\mu_{0c}} and squared scale multiplied by \eqn{r}:
#'     \eqn{\mu_c \sim t_{\nu_{nt}}(\mu_{0c}, r \sigma_{nt}^2 / \kappa_{nt})}.
#'     \code{prior_c} is not used and must be \code{NULL}.}
#'   \item{Fixed control}{\code{n_c = 0} and \code{fixed_c} is given: the
#'     control mean is a known value \eqn{\mu_c}, as in a single-arm study
#'     compared with a benchmark. The probabilities of the Go and NoGo regions
#'     are computed exactly from the t distribution of the treatment mean
#'     (posterior) or of the mean of \eqn{m_t} future subjects (posterior
#'     predictive), shifted by \eqn{\mu_c}; \code{method}, \code{nMC} and
#'     \code{m_c} are not used. \code{prior_c} is not used and must be
#'     \code{NULL}.}
#' }
#'
#' For the posterior predictive criterion (\code{\link{criterion_predictive}}),
#' the future effect is the difference between the sample means of \eqn{m_t}
#' and \eqn{m_c} future subjects. The posterior predictive distribution of the
#' mean of \eqn{m} future subjects of a group is
#' \eqn{t_{\nu_{n}}(\mu_{n}, \sigma_{n}^2 (1 / \kappa_{n} + 1 / m))}; for a
#' hypothetical control it is
#' \eqn{t_{\nu_{nt}}(\mu_{0c}, r \sigma_{nt}^2 (1 / \kappa_{nt} + 1 / m_c))}.
#'
#' The probabilities of the Go and NoGo regions are probabilities of the
#' difference of two t variables, computed by \code{\link{ptdiff}} with the
#' chosen \code{method}. With \code{method = "MC"}, the probabilities of the
#' two regions are estimated from the same \code{nMC} draws.
#'
#' @param n_t Number of subjects in the treatment group (a whole number of at
#'   least 2).
#' @param n_c Number of subjects in the concurrent control group: 0 (no
#'   concurrent control group) or a whole number of at least 2.
#' @param prior_t,prior_c Prior distributions of the mean and variance of the
#'   treatment and control groups, created by \code{\link{prior_nix}} or
#'   \code{\link{prior_vague}}. \code{prior_c} must be \code{NULL} (the
#'   default) for a hypothetical or fixed control and is required otherwise.
#' @param criterion Decision criterion for one endpoint, created by
#'   \code{\link{criterion_posterior}} or \code{\link{criterion_predictive}}.
#'   Its argument \code{direction} sets whether larger or smaller effects
#'   are favorable.
#' @param external_t,external_c External data of the treatment and control
#'   groups, created by \code{\link{borrow}} with
#'   \code{\link{data_continuous}}, or \code{NULL} (no external data).
#' @param hypothetical_c Hypothetical control, created by
#'   \code{\link{hypothetical}} with \code{mean} (of length 1) and \code{r},
#'   or \code{NULL}. It can be used only when \code{n_c = 0}.
#' @param fixed_c Fixed control mean, a number, or \code{NULL}. It can be used
#'   only when \code{n_c = 0}.
#' @param method Method used to compute the probabilities: \code{"NI"}
#'   (numerical integration, the default), \code{"MM"} (moment-matching
#'   approximation, which requires the distributions of both group means to
#'   have more than 4 degrees of freedom) or \code{"MC"} (Monte Carlo
#'   simulation). See \code{\link{ptdiff}}. Numerical integration is exact
#'   up to a small numerical error but is the slowest of the deterministic
#'   methods: for operating characteristics, which need the probabilities for
#'   many simulated studies, \code{"MM"} is faster by about two orders of
#'   magnitude and is suited to exploring designs, and \code{"NI"} can be
#'   used to confirm the chosen design. The computing time of \code{"MC"} is
#'   proportional to \code{nMC}.
#' @param nMC Number of Monte Carlo draws for \code{method = "MC"}.
#'
#' @return An object of class \code{c("qdm_design_continuous", "qdm_design")},
#'   a list with the arguments and the element \code{control}
#'   (\code{"concurrent"}, \code{"external"}, \code{"hypothetical"} or
#'   \code{"fixed"}).
#'
#' @references
#' Kang, R., Yamaguchi, Y. and Han, C. (2026). Bayesian quantitative
#' decision-making framework in cell and gene therapy development for rare
#' diseases. \emph{Journal of Biopharmaceutical Statistics}.
#' \doi{10.1080/10543406.2026.2655410}
#'
#' Huang, Y., Yamaguchi, Y., Homma, G., Takeda, K. and Maruo, K. (submitted).
#' Conjugate representation of power priors for efficient Bayesian analysis
#' of normal data.
#'
#' @family designs
#' @examples
#' # Concurrent control with external treatment and control data
#' design_continuous(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_vague(), prior_c = prior_vague(),
#'   criterion = criterion_posterior(theta_tv = 30, theta_mav = 15,
#'                                   gamma_go = 0.6, gamma_nogo = 0.7),
#'   external_t = borrow(data_continuous(n = 6, mean = 20, sd = 30), a0 = 0.5),
#'   external_c = borrow(data_continuous(n = 30, mean = -12, sd = 20),
#'                       a0 = 0.4)
#' )
#'
#' # Uncontrolled study with a hypothetical control
#' design_continuous(
#'   n_t = 12, n_c = 0, prior_t = prior_vague(),
#'   criterion = criterion_predictive(theta_null = 0, m_t = 50,
#'                                    gamma_go = 0.9, gamma_nogo = 0.7),
#'   hypothetical_c = hypothetical(mean = -10, r = 1),
#'   method = "MM"
#' )
#'
#' # Single-arm study compared with a fixed control mean of -10
#' design_continuous(
#'   n_t = 12, n_c = 0, prior_t = prior_vague(),
#'   criterion = criterion_posterior(theta_tv = 30, theta_mav = 15,
#'                                   gamma_go = 0.6, gamma_nogo = 0.7),
#'   fixed_c = -10
#' )
#' @export
design_continuous <- function(n_t, n_c, prior_t, prior_c = NULL, criterion,
                              external_t = NULL, external_c = NULL,
                              hypothetical_c = NULL, fixed_c = NULL,
                              method = c("NI", "MM", "MC"), nMC = 10000L) {
  check_numeric(n_t, "n_t", lower = 2, integer = TRUE)
  check_numeric(n_c, "n_c", lower = 0, integer = TRUE)
  if (round(n_c) == 1) {
    stop("'n_c' must be 0 or a whole number of at least 2.", call. = FALSE)
  }
  check_prior_continuous(prior_t, "prior_t")
  if (!inherits(criterion, "qdm_criterion") || criterion$n_endpoints != 1L) {
    stop(paste("'criterion' must be created by criterion_posterior() or",
               "criterion_predictive() for one endpoint."), call. = FALSE)
  }
  check_external_continuous(external_t, "external_t")
  check_external_continuous(external_c, "external_c")
  if (!is.null(hypothetical_c) &&
      (!inherits(hypothetical_c, "qdm_hypothetical") ||
       hypothetical_c$type != "location" ||
       length(hypothetical_c$mean) != 1L)) {
    stop(paste("'hypothetical_c' must be NULL or created by hypothetical()",
               "with 'mean' of length 1 and 'r'."), call. = FALSE)
  }
  fixed_c <- check_fixed_c(fixed_c, 1L, probability = FALSE)
  control <- control_form(n_c, external_c, hypothetical_c, fixed_c)
  check_prior_control(prior_c, control, check_prior_continuous,
                      unused = c("hypothetical", "fixed"))
  method <- match.arg(method)
  check_numeric(nMC, "nMC", lower = 1, upper = .Machine$integer.max,
                integer = TRUE)
  des <- structure(list(endpoint = "continuous",
                        n_t = as.integer(round(n_t)),
                        n_c = as.integer(round(n_c)), prior_t = prior_t,
                        prior_c = prior_c, criterion = criterion,
                        external_t = external_t, external_c = external_c,
                        hypothetical_c = hypothetical_c, fixed_c = fixed_c,
                        control = control, method = method,
                        nMC = as.integer(round(nMC))),
                   class = c("qdm_design_continuous", "qdm_design"))
  if (control == "fixed") {
    # The probabilities are computed exactly, whatever the method
    return(des)
  }
  # Degrees of freedom of the distributions of the group means for the
  # sample sizes of the design (they do not depend on the data values)
  post_t <- posterior_continuous(prior_t, des$n_t, external = external_t)
  post_c <- control_posterior_continuous(des, post_t, des$n_c)
  if (post_c$nu <= 0) {
    stop(paste("With prior_vague() for the control group and external",
               "control data only, 'a0' times the number of external",
               "subjects must be greater than 1, so that the distribution",
               "of the control mean is proper."), call. = FALSE)
  }
  check_df_mm(post_t$nu, post_c$nu, method)
  des
}
