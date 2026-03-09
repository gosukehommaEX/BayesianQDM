#' Find Optimal Go/NoGo Thresholds for a Single Continuous Endpoint
#'
#' Computes the optimal Go threshold \eqn{\gamma_{\mathrm{go}}} and NoGo
#' threshold \eqn{\gamma_{\mathrm{nogo}}} for a single continuous endpoint by
#' searching over a grid of candidate values.  The two thresholds are
#' calibrated independently under separate scenarios:
#' \itemize{
#'   \item \eqn{\gamma_{\mathrm{go}}} is the \strong{smallest} value in
#'         \code{gamma_grid} such that the marginal Go probability
#'         \eqn{\Pr(g_{\mathrm{Go}} \ge \gamma_{\mathrm{go}})} is strictly
#'         less than \code{target_go} under the Go-calibration scenario
#'         (\code{mu_t_go}, \code{mu_c_go}, \code{sigma_t_go},
#'         \code{sigma_c_go}); typically the Null scenario.
#'   \item \eqn{\gamma_{\mathrm{nogo}}} is the \strong{smallest} value in
#'         \code{gamma_grid} such that the marginal NoGo probability
#'         \eqn{\Pr(g_{\mathrm{NoGo}} \ge \gamma_{\mathrm{nogo}})} is
#'         strictly less than \code{target_nogo} under the NoGo-calibration
#'         scenario (\code{mu_t_nogo}, \code{mu_c_nogo}, \code{sigma_t_nogo},
#'         \code{sigma_c_nogo}); typically the Alternative scenario.
#' }
#' Here \eqn{g_{\mathrm{Go}} = P(\theta > \theta_{\mathrm{TV}} \mid \bar{y}_t, s_t, \bar{y}_c, s_c)}
#' and \eqn{g_{\mathrm{NoGo}} = P(\theta \le \theta_{\mathrm{MAV}} \mid \bar{y}_t, s_t, \bar{y}_c, s_c)}
#' for \code{prob = 'posterior'}, consistent with the decision rule in
#' \code{\link{pbayesdecisionprob1cont}}.
#'
#' @param nsim A positive integer giving the number of Monte Carlo simulation
#'        replicates used to approximate the operating characteristics.
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or
#'        \code{'external'}.
#' @param prior A character string specifying the prior distribution.
#'        Must be \code{'vague'} or \code{'N-Inv-Chisq'}.
#' @param CalcMethod A character string specifying the computation method
#'        for \code{\link{pbayespostpred1cont}}.
#'        Must be \code{'NI'}, \code{'MC'}, or \code{'MM'}.
#' @param theta_TV A numeric scalar giving the Target Value (TV) threshold
#'        for the treatment effect.  Required when \code{prob = 'posterior'};
#'        set to \code{NULL} otherwise.
#' @param theta_MAV A numeric scalar giving the Minimum Acceptable Value (MAV)
#'        threshold.  Must satisfy \code{theta_TV > theta_MAV}.  Required when
#'        \code{prob = 'posterior'}; set to \code{NULL} otherwise.
#' @param theta_NULL A numeric scalar giving the null hypothesis threshold
#'        for the predictive probability.  Required when
#'        \code{prob = 'predictive'}; set to \code{NULL} otherwise.
#' @param nMC A positive integer giving the number of inner Monte Carlo draws
#'        passed to \code{\link{pbayespostpred1cont}}.  Required when
#'        \code{CalcMethod = 'MC'}; set to \code{NULL} otherwise.
#' @param mu_t_go A numeric scalar giving the true mean for the treatment group
#'        under the Go-calibration scenario (typically Null).
#' @param mu_c_go A numeric scalar giving the true mean for the control group
#'        under the Go-calibration scenario.  Set to \code{NULL} for
#'        \code{design = 'uncontrolled'}.
#' @param sigma_t_go A positive numeric scalar giving the true standard
#'        deviation for the treatment group under the Go-calibration scenario.
#' @param sigma_c_go A positive numeric scalar giving the true standard
#'        deviation for the control group under the Go-calibration scenario.
#'        Set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param mu_t_nogo A numeric scalar giving the true mean for the treatment
#'        group under the NoGo-calibration scenario (typically Alternative).
#' @param mu_c_nogo A numeric scalar giving the true mean for the control
#'        group under the NoGo-calibration scenario.  Set to \code{NULL} for
#'        \code{design = 'uncontrolled'}.
#' @param sigma_t_nogo A positive numeric scalar giving the true standard
#'        deviation for the treatment group under the NoGo-calibration scenario.
#' @param sigma_c_nogo A positive numeric scalar giving the true standard
#'        deviation for the control group under the NoGo-calibration scenario.
#'        Set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param target_go A numeric scalar in \code{(0, 1)} giving the upper bound
#'        on the marginal Go probability under the Go-calibration scenario.
#'        The optimal \eqn{\gamma_{\mathrm{go}}} is the smallest grid value
#'        satisfying \eqn{\Pr(\mathrm{Go}) < \code{target\_go}}.
#' @param target_nogo A numeric scalar in \code{(0, 1)} giving the upper bound
#'        on the marginal NoGo probability under the NoGo-calibration scenario.
#'        The optimal \eqn{\gamma_{\mathrm{nogo}}} is the smallest grid value
#'        satisfying \eqn{\Pr(\mathrm{NoGo}) < \code{target\_nogo}}.
#' @param n_t A positive integer giving the number of patients in the
#'        treatment group in the PoC trial.
#' @param n_c A positive integer giving the number of patients in the
#'        control group in the PoC trial.  Required for
#'        \code{design = 'controlled'} or \code{'external'}; set to
#'        \code{NULL} for \code{design = 'uncontrolled'}.
#' @param m_t A positive integer giving the future sample size for the
#'        treatment group. Required when \code{prob = 'predictive'};
#'        set to \code{NULL} otherwise.
#' @param m_c A positive integer giving the future sample size for the
#'        control group. Required when \code{prob = 'predictive'};
#'        set to \code{NULL} otherwise.
#' @param kappa0_t A positive numeric scalar giving the prior precision
#'        parameter for the treatment group under the N-Inv-Chi-squared prior.
#'        Required when \code{prior = 'N-Inv-Chisq'}; set to \code{NULL} otherwise.
#' @param kappa0_c A positive numeric scalar giving the prior precision
#'        parameter for the control group.  Required when
#'        \code{prior = 'N-Inv-Chisq'} and \code{design != 'uncontrolled'};
#'        set to \code{NULL} otherwise.
#' @param nu0_t A positive numeric scalar giving the prior degrees of freedom
#'        for the treatment group under the N-Inv-Chi-squared prior.  Required
#'        when \code{prior = 'N-Inv-Chisq'}; set to \code{NULL} otherwise.
#' @param nu0_c A positive numeric scalar giving the prior degrees of freedom
#'        for the control group.  Required when \code{prior = 'N-Inv-Chisq'}
#'        and \code{design != 'uncontrolled'}; set to \code{NULL} otherwise.
#' @param mu0_t A numeric scalar giving the prior mean for the treatment group
#'        under the N-Inv-Chi-squared prior.  Required when
#'        \code{prior = 'N-Inv-Chisq'}; set to \code{NULL} otherwise.
#' @param mu0_c A numeric scalar giving the prior mean for the control group
#'        (or the fixed hypothetical control mean for uncontrolled design).
#'        Required when \code{prior = 'N-Inv-Chisq'}; set to \code{NULL}
#'        otherwise.
#' @param sigma0_t A positive numeric scalar giving the prior scale parameter
#'        for the treatment group under the N-Inv-Chi-squared prior.  Required
#'        when \code{prior = 'N-Inv-Chisq'}; set to \code{NULL} otherwise.
#' @param sigma0_c A positive numeric scalar giving the prior scale parameter
#'        for the control group.  Required when \code{prior = 'N-Inv-Chisq'}
#'        and \code{design != 'uncontrolled'}; set to \code{NULL} otherwise.
#' @param r A positive numeric scalar giving the variance ratio
#'        \eqn{\sigma_c^2 / \sigma_t^2} used for uncontrolled design.
#'        Required when \code{design = 'uncontrolled'}; set to \code{NULL}
#'        otherwise.
#' @param ne_t A positive integer giving the number of patients in the
#'        treatment group of the external data set.  Required when
#'        \code{design = 'external'} and external treatment data are
#'        available; otherwise set to \code{NULL}.
#' @param ne_c A positive integer giving the number of patients in the
#'        control group of the external data set.  Required when
#'        \code{design = 'external'} and external control data are available;
#'        otherwise set to \code{NULL}.
#' @param alpha0e_t A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for external treatment group data.  Required when
#'        \code{design = 'external'} and \code{ne_t} is non-NULL; set to
#'        \code{NULL} otherwise.
#' @param alpha0e_c A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for external control group data.  Required when
#'        \code{design = 'external'} and \code{ne_c} is non-NULL; set to
#'        \code{NULL} otherwise.
#' @param bar_ye_t A numeric scalar giving the sample mean of the external
#'        treatment group data.  Required when \code{ne_t} is non-NULL; set to
#'        \code{NULL} otherwise.
#' @param bar_ye_c A numeric scalar giving the sample mean of the external
#'        control group data.  Required when \code{ne_c} is non-NULL; set to
#'        \code{NULL} otherwise.
#' @param se_t A positive numeric scalar giving the sample standard deviation
#'        of the external treatment group data.  Required when \code{ne_t} is
#'        non-NULL; set to \code{NULL} otherwise.
#' @param se_c A positive numeric scalar giving the sample standard deviation
#'        of the external control group data.  Required when \code{ne_c} is
#'        non-NULL; set to \code{NULL} otherwise.
#' @param gamma_grid A numeric vector of candidate threshold values in
#'        \code{(0, 1)} to search over.  Defaults to
#'        \code{seq(0.01, 0.99, by = 0.01)}.
#' @param seed A numeric scalar for reproducible random number generation.
#'        The Go-calibration simulation uses \code{seed} and the
#'        NoGo-calibration simulation uses \code{seed + 1} to ensure
#'        independence between the two scenarios.
#'
#' @return A list of class \code{getgamma1cont} with the following elements:
#' \describe{
#'   \item{gamma_go}{Optimal Go threshold: the smallest value in
#'         \code{gamma_grid} for which \eqn{\Pr(\mathrm{Go}) < \code{target\_go}}
#'         under the Go-calibration scenario.
#'         \code{NA} if no such value exists.}
#'   \item{gamma_nogo}{Optimal NoGo threshold: the smallest value in
#'         \code{gamma_grid} for which \eqn{\Pr(\mathrm{NoGo}) < \code{target\_nogo}}
#'         under the NoGo-calibration scenario.
#'         \code{NA} if no such value exists.}
#'   \item{PrGo_opt}{Marginal \eqn{\Pr(g_{\mathrm{Go}} \ge \gamma_{\mathrm{go}})}
#'         at the optimal \eqn{\gamma_{\mathrm{go}}} under the Go-calibration
#'         scenario.  \code{NA} if \code{gamma_go} is \code{NA}.}
#'   \item{PrNoGo_opt}{Marginal \eqn{\Pr(g_{\mathrm{NoGo}} \ge \gamma_{\mathrm{nogo}})}
#'         at the optimal \eqn{\gamma_{\mathrm{nogo}}} under the NoGo-calibration
#'         scenario.  \code{NA} if \code{gamma_nogo} is \code{NA}.}
#'   \item{grid_results}{A data frame with columns \code{gamma_grid},
#'         \code{PrGo_grid} (marginal Go probability under the Go-calibration
#'         scenario), and \code{PrNoGo_grid} (marginal NoGo probability under
#'         the NoGo-calibration scenario).}
#' }
#'
#' @details
#' The function uses a two-stage simulate-then-sweep strategy:
#' \enumerate{
#'   \item \strong{Simulation}: \code{nsim} datasets are generated
#'         independently for each calibration scenario.  For the
#'         Go-calibration scenario, standardised residuals are shifted by
#'         \code{mu_t_go} and \code{mu_c_go}; for the NoGo-calibration
#'         scenario, by \code{mu_t_nogo} and \code{mu_c_nogo}.
#'         \code{\link{pbayespostpred1cont}} is called once per scenario to
#'         obtain \eqn{g_{\mathrm{Go}}} (\code{lower.tail = FALSE} at
#'         \code{theta_TV}) and \eqn{g_{\mathrm{NoGo}}} (\code{lower.tail
#'         = TRUE} at \code{theta_MAV}).
#'   \item \strong{Gamma sweep}: Marginal probabilities are estimated as the
#'         proportion of simulated datasets satisfying the respective
#'         indicator: \eqn{g_{\mathrm{Go}} \ge \gamma} for
#'         \code{PrGo_grid}, and \eqn{g_{\mathrm{NoGo}} \ge \gamma} for
#'         \code{PrNoGo_grid}.  Both are monotone non-increasing in
#'         \eqn{\gamma}.
#' }
#' The optimal \eqn{\gamma_{\mathrm{go}}} and \eqn{\gamma_{\mathrm{nogo}}}
#' are each the \emph{smallest} grid value crossing below the respective
#' target probability.
#'
#' @examples
#' # Example 1: Controlled design, vague prior, posterior probability
#' # gamma_go  : smallest gamma s.t. Pr(Go)   < 0.05 under Null (mu_t = mu_c = 1.0)
#' # gamma_nogo: smallest gamma s.t. Pr(NoGo) < 0.20 under Alt  (mu_t = 2.5, mu_c = 1.0)
#' getgamma1cont(
#'   nsim = 1000L, prob = 'posterior', design = 'controlled',
#'   prior = 'vague', CalcMethod = 'MM',
#'   theta_TV = 1.5, theta_MAV = 0.0, theta_NULL = NULL, nMC = NULL,
#'   mu_t_go = 1.0, mu_c_go = 1.0, sigma_t_go = 2.0, sigma_c_go = 2.0,
#'   mu_t_nogo = 2.5, mu_c_nogo = 1.0, sigma_t_nogo = 2.0, sigma_c_nogo = 2.0,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 15L, n_c = 15L, m_t = NULL, m_c = NULL,
#'   kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
#'   mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
#'   r = NULL, ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   gamma_grid = seq(0.01, 0.99, by = 0.01), seed = 1L
#' )
#'
#' # Example 2: Uncontrolled design, N-Inv-Chisq prior, posterior probability
#' getgamma1cont(
#'   nsim = 1000L, prob = 'posterior', design = 'uncontrolled',
#'   prior = 'N-Inv-Chisq', CalcMethod = 'NI',
#'   theta_TV = 1.0, theta_MAV = 0.0, theta_NULL = NULL, nMC = NULL,
#'   mu_t_go = 1.5, mu_c_go = NULL, sigma_t_go = 1.5, sigma_c_go = NULL,
#'   mu_t_nogo = 3.0, mu_c_nogo = NULL, sigma_t_nogo = 1.5, sigma_c_nogo = NULL,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 20L, n_c = NULL, m_t = NULL, m_c = NULL,
#'   kappa0_t = 2, kappa0_c = NULL, nu0_t = 5, nu0_c = NULL,
#'   mu0_t = 3.0, mu0_c = 1.5, sigma0_t = 1.5, sigma0_c = NULL,
#'   r = 1.0, ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   gamma_grid = seq(0.01, 0.99, by = 0.01), seed = 2L
#' )
#'
#' # Example 3: External design, vague prior, posterior probability
#' getgamma1cont(
#'   nsim = 1000L, prob = 'posterior', design = 'external',
#'   prior = 'vague', CalcMethod = 'MM',
#'   theta_TV = 1.0, theta_MAV = 0.0, theta_NULL = NULL, nMC = NULL,
#'   mu_t_go = 1.0, mu_c_go = 1.0, sigma_t_go = 1.5, sigma_c_go = 1.5,
#'   mu_t_nogo = 2.5, mu_c_nogo = 1.0, sigma_t_nogo = 1.5, sigma_c_nogo = 1.5,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 15L, n_c = 15L, m_t = NULL, m_c = NULL,
#'   kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
#'   mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
#'   r = NULL, ne_t = NULL, ne_c = 20L, alpha0e_t = NULL, alpha0e_c = 0.5,
#'   bar_ye_t = NULL, bar_ye_c = 0.0, se_t = NULL, se_c = 1.5,
#'   gamma_grid = seq(0.01, 0.99, by = 0.01), seed = 4L
#' )
#'
#' # Example 4: Controlled design, vague prior, predictive probability
#' getgamma1cont(
#'   nsim = 1000L, prob = 'predictive', design = 'controlled',
#'   prior = 'vague', CalcMethod = 'MM',
#'   theta_TV = NULL, theta_MAV = NULL, theta_NULL = 1.0, nMC = NULL,
#'   mu_t_go = 1.0, mu_c_go = 1.0, sigma_t_go = 2.0, sigma_c_go = 2.0,
#'   mu_t_nogo = 2.5, mu_c_nogo = 1.0, sigma_t_nogo = 2.0, sigma_c_nogo = 2.0,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 15L, n_c = 15L, m_t = 50L, m_c = 50L,
#'   kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
#'   mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
#'   r = NULL, ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   gamma_grid = seq(0.01, 0.99, by = 0.01), seed = 3L
#' )
#'
#' # Example 5: Uncontrolled design, vague prior, predictive probability
#' getgamma1cont(
#'   nsim = 1000L, prob = 'predictive', design = 'uncontrolled',
#'   prior = 'vague', CalcMethod = 'MM',
#'   theta_TV = NULL, theta_MAV = NULL, theta_NULL = 1.0, nMC = NULL,
#'   mu_t_go = 1.0, mu_c_go = NULL, sigma_t_go = 2.0, sigma_c_go = NULL,
#'   mu_t_nogo = 2.5, mu_c_nogo = NULL, sigma_t_nogo = 2.0, sigma_c_nogo = NULL,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 15L, n_c = NULL, m_t = 50L, m_c = 50L,
#'   kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
#'   mu0_t = NULL, mu0_c = 0.0, sigma0_t = NULL, sigma0_c = NULL,
#'   r = 1.0, ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   gamma_grid = seq(0.01, 0.99, by = 0.01), seed = 5L
#' )
#'
#' # Example 6: External design, vague prior, predictive probability
#' getgamma1cont(
#'   nsim = 1000L, prob = 'predictive', design = 'external',
#'   prior = 'vague', CalcMethod = 'MM',
#'   theta_TV = NULL, theta_MAV = NULL, theta_NULL = 1.0, nMC = NULL,
#'   mu_t_go = 1.0, mu_c_go = 1.0, sigma_t_go = 1.5, sigma_c_go = 1.5,
#'   mu_t_nogo = 2.5, mu_c_nogo = 1.0, sigma_t_nogo = 1.5, sigma_c_nogo = 1.5,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 15L, n_c = 15L, m_t = 50L, m_c = 50L,
#'   kappa0_t = NULL, kappa0_c = NULL, nu0_t = NULL, nu0_c = NULL,
#'   mu0_t = NULL, mu0_c = NULL, sigma0_t = NULL, sigma0_c = NULL,
#'   r = NULL, ne_t = NULL, ne_c = 20L, alpha0e_t = NULL, alpha0e_c = 0.5,
#'   bar_ye_t = NULL, bar_ye_c = 0.0, se_t = NULL, se_c = 1.5,
#'   gamma_grid = seq(0.01, 0.99, by = 0.01), seed = 6L
#' )
#'
#' @importFrom stats rnorm
#' @export
getgamma1cont <- function(nsim,
                          prob = 'posterior', design = 'controlled',
                          prior = 'vague', CalcMethod = 'NI',
                          theta_TV = NULL, theta_MAV = NULL, theta_NULL = NULL,
                          nMC = NULL,
                          mu_t_go, mu_c_go = NULL,
                          sigma_t_go, sigma_c_go = NULL,
                          mu_t_nogo, mu_c_nogo = NULL,
                          sigma_t_nogo, sigma_c_nogo = NULL,
                          target_go, target_nogo,
                          n_t, n_c = NULL,
                          m_t = NULL, m_c = NULL,
                          kappa0_t = NULL, kappa0_c = NULL,
                          nu0_t = NULL, nu0_c = NULL,
                          mu0_t = NULL, mu0_c = NULL,
                          sigma0_t = NULL, sigma0_c = NULL,
                          r = NULL,
                          ne_t = NULL, ne_c = NULL,
                          alpha0e_t = NULL, alpha0e_c = NULL,
                          bar_ye_t = NULL, bar_ye_c = NULL,
                          se_t = NULL, se_c = NULL,
                          gamma_grid = seq(0.01, 0.99, by = 0.01),
                          seed) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) ||
      nsim != floor(nsim) || nsim < 1L) {
    stop("'nsim' must be a single positive integer")
  }

  if (!is.character(prob) || length(prob) != 1L ||
      !prob %in% c('posterior', 'predictive')) {
    stop("'prob' must be either 'posterior' or 'predictive'")
  }

  if (!is.character(design) || length(design) != 1L ||
      !design %in% c('controlled', 'uncontrolled', 'external')) {
    stop("'design' must be 'controlled', 'uncontrolled', or 'external'")
  }

  if (!is.character(prior) || length(prior) != 1L ||
      !prior %in% c('vague', 'N-Inv-Chisq')) {
    stop("'prior' must be either 'vague' or 'N-Inv-Chisq'")
  }

  if (!is.character(CalcMethod) || length(CalcMethod) != 1L ||
      !CalcMethod %in% c('NI', 'MC', 'MM')) {
    stop("'CalcMethod' must be 'NI', 'MC', or 'MM'")
  }

  if (CalcMethod == 'MC') {
    if (is.null(nMC) || !is.numeric(nMC) || length(nMC) != 1L ||
        is.na(nMC) || nMC != floor(nMC) || nMC < 1L) {
      stop("'nMC' must be a single positive integer when CalcMethod = 'MC'")
    }
  }

  if (prob == 'posterior') {
    if (is.null(theta_TV) || is.null(theta_MAV)) {
      stop("'theta_TV' and 'theta_MAV' must be non-NULL when prob = 'posterior'")
    }
    if (!is.numeric(theta_TV) || length(theta_TV) != 1L || is.na(theta_TV)) {
      stop("'theta_TV' must be a single numeric value")
    }
    if (!is.numeric(theta_MAV) || length(theta_MAV) != 1L || is.na(theta_MAV)) {
      stop("'theta_MAV' must be a single numeric value")
    }
    if (theta_TV <= theta_MAV) {
      stop("'theta_TV' must be strictly greater than 'theta_MAV'")
    }
  } else {
    if (is.null(theta_NULL)) {
      stop("'theta_NULL' must be non-NULL when prob = 'predictive'")
    }
    if (!is.numeric(theta_NULL) || length(theta_NULL) != 1L ||
        is.na(theta_NULL)) {
      stop("'theta_NULL' must be a single numeric value")
    }
    if (is.null(m_t) || is.null(m_c)) {
      stop("'m_t' and 'm_c' must be non-NULL when prob = 'predictive'")
    }
    for (nm in c("m_t", "m_c")) {
      val <- get(nm)
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 1L) {
        stop(paste0("'", nm, "' must be a single positive integer"))
      }
    }
  }

  # Validate Go-calibration scenario parameters
  if (!is.numeric(mu_t_go) || length(mu_t_go) != 1L || is.na(mu_t_go)) {
    stop("'mu_t_go' must be a single numeric value")
  }
  if (!is.numeric(sigma_t_go) || length(sigma_t_go) != 1L ||
      is.na(sigma_t_go) || sigma_t_go <= 0) {
    stop("'sigma_t_go' must be a single positive numeric value")
  }

  # Validate NoGo-calibration scenario parameters
  if (!is.numeric(mu_t_nogo) || length(mu_t_nogo) != 1L || is.na(mu_t_nogo)) {
    stop("'mu_t_nogo' must be a single numeric value")
  }
  if (!is.numeric(sigma_t_nogo) || length(sigma_t_nogo) != 1L ||
      is.na(sigma_t_nogo) || sigma_t_nogo <= 0) {
    stop("'sigma_t_nogo' must be a single positive numeric value")
  }

  if (design != 'uncontrolled') {
    if (is.null(mu_c_go) || !is.numeric(mu_c_go) || length(mu_c_go) != 1L ||
        is.na(mu_c_go)) {
      stop("'mu_c_go' must be a single numeric value for controlled or external design")
    }
    if (is.null(sigma_c_go) || !is.numeric(sigma_c_go) ||
        length(sigma_c_go) != 1L || is.na(sigma_c_go) || sigma_c_go <= 0) {
      stop("'sigma_c_go' must be a single positive numeric value for controlled or external design")
    }
    if (is.null(mu_c_nogo) || !is.numeric(mu_c_nogo) || length(mu_c_nogo) != 1L ||
        is.na(mu_c_nogo)) {
      stop("'mu_c_nogo' must be a single numeric value for controlled or external design")
    }
    if (is.null(sigma_c_nogo) || !is.numeric(sigma_c_nogo) ||
        length(sigma_c_nogo) != 1L || is.na(sigma_c_nogo) || sigma_c_nogo <= 0) {
      stop("'sigma_c_nogo' must be a single positive numeric value for controlled or external design")
    }
  }

  for (nm in c("target_go", "target_nogo")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val <= 0 || val >= 1) {
      stop(paste0("'", nm, "' must be a single numeric value in (0, 1)"))
    }
  }

  if (!is.numeric(n_t) || length(n_t) != 1L || is.na(n_t) ||
      n_t != floor(n_t) || n_t < 1L) {
    stop("'n_t' must be a single positive integer")
  }

  if (design != 'uncontrolled') {
    if (is.null(n_c) || !is.numeric(n_c) || length(n_c) != 1L ||
        is.na(n_c) || n_c != floor(n_c) || n_c < 1L) {
      stop("'n_c' must be a single positive integer for controlled or external design")
    }
  }

  if (prior == 'N-Inv-Chisq') {
    for (nm in c("kappa0_t", "nu0_t", "mu0_t", "sigma0_t")) {
      val <- get(nm)
      if (is.null(val)) {
        stop(paste0("'", nm, "' must be non-NULL when prior = 'N-Inv-Chisq'"))
      }
    }
    if (!is.numeric(kappa0_t) || length(kappa0_t) != 1L || is.na(kappa0_t) ||
        kappa0_t <= 0) {
      stop("'kappa0_t' must be a single positive numeric value")
    }
    if (!is.numeric(nu0_t) || length(nu0_t) != 1L || is.na(nu0_t) ||
        nu0_t <= 0) {
      stop("'nu0_t' must be a single positive numeric value")
    }
    if (!is.numeric(mu0_t) || length(mu0_t) != 1L || is.na(mu0_t)) {
      stop("'mu0_t' must be a single numeric value")
    }
    if (!is.numeric(sigma0_t) || length(sigma0_t) != 1L || is.na(sigma0_t) ||
        sigma0_t <= 0) {
      stop("'sigma0_t' must be a single positive numeric value")
    }
    if (design != 'uncontrolled') {
      for (nm in c("kappa0_c", "nu0_c", "mu0_c", "sigma0_c")) {
        val <- get(nm)
        if (is.null(val)) {
          stop(paste0("'", nm, "' must be non-NULL when prior = 'N-Inv-Chisq' and design != 'uncontrolled'"))
        }
      }
      if (!is.numeric(kappa0_c) || length(kappa0_c) != 1L || is.na(kappa0_c) ||
          kappa0_c <= 0) {
        stop("'kappa0_c' must be a single positive numeric value")
      }
      if (!is.numeric(nu0_c) || length(nu0_c) != 1L || is.na(nu0_c) ||
          nu0_c <= 0) {
        stop("'nu0_c' must be a single positive numeric value")
      }
      if (!is.numeric(mu0_c) || length(mu0_c) != 1L || is.na(mu0_c)) {
        stop("'mu0_c' must be a single numeric value")
      }
      if (!is.numeric(sigma0_c) || length(sigma0_c) != 1L || is.na(sigma0_c) ||
          sigma0_c <= 0) {
        stop("'sigma0_c' must be a single positive numeric value")
      }
    } else {
      if (is.null(mu0_c) || !is.numeric(mu0_c) || length(mu0_c) != 1L ||
          is.na(mu0_c)) {
        stop("'mu0_c' must be a single numeric value when design = 'uncontrolled'")
      }
    }
  }

  if (design == 'uncontrolled') {
    if (is.null(r) || !is.numeric(r) || length(r) != 1L || is.na(r) ||
        r <= 0) {
      stop("'r' must be a single positive numeric value when design = 'uncontrolled'")
    }
  }

  if (design == 'external') {
    has_ext_t <- !is.null(ne_t) && !is.null(alpha0e_t) && !is.null(bar_ye_t) &&
      !is.null(se_t)
    has_ext_c <- !is.null(ne_c) && !is.null(alpha0e_c) && !is.null(bar_ye_c) &&
      !is.null(se_c)
    if (!has_ext_t && !has_ext_c) {
      stop("For design = 'external', at least one of (ne_t, alpha0e_t, bar_ye_t, se_t) or (ne_c, alpha0e_c, bar_ye_c, se_c) must be fully specified")
    }
    if (!is.null(ne_t)) {
      if (!is.numeric(ne_t) || length(ne_t) != 1L || is.na(ne_t) ||
          ne_t != floor(ne_t) || ne_t < 1L) {
        stop("'ne_t' must be a single positive integer")
      }
      if (!is.numeric(alpha0e_t) || length(alpha0e_t) != 1L ||
          is.na(alpha0e_t) || alpha0e_t <= 0 || alpha0e_t > 1) {
        stop("'alpha0e_t' must be a single numeric value in (0, 1]")
      }
      if (!is.numeric(bar_ye_t) || length(bar_ye_t) != 1L || is.na(bar_ye_t)) {
        stop("'bar_ye_t' must be a single numeric value")
      }
      if (!is.numeric(se_t) || length(se_t) != 1L || is.na(se_t) || se_t <= 0) {
        stop("'se_t' must be a single positive numeric value")
      }
    }
    if (!is.null(ne_c)) {
      if (!is.numeric(ne_c) || length(ne_c) != 1L || is.na(ne_c) ||
          ne_c != floor(ne_c) || ne_c < 1L) {
        stop("'ne_c' must be a single positive integer")
      }
      if (!is.numeric(alpha0e_c) || length(alpha0e_c) != 1L ||
          is.na(alpha0e_c) || alpha0e_c <= 0 || alpha0e_c > 1) {
        stop("'alpha0e_c' must be a single numeric value in (0, 1]")
      }
      if (!is.numeric(bar_ye_c) || length(bar_ye_c) != 1L || is.na(bar_ye_c)) {
        stop("'bar_ye_c' must be a single numeric value")
      }
      if (!is.numeric(se_c) || length(se_c) != 1L || is.na(se_c) || se_c <= 0) {
        stop("'se_c' must be a single positive numeric value")
      }
    }
  }

  if (!is.numeric(gamma_grid) || length(gamma_grid) < 1L ||
      any(is.na(gamma_grid)) || any(gamma_grid <= 0) || any(gamma_grid >= 1)) {
    stop("'gamma_grid' must be a numeric vector with all values in (0, 1)")
  }
  gamma_grid <- sort(unique(gamma_grid))

  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("'seed' must be a single numeric value")
  }

  # ---------------------------------------------------------------------------
  # Internal helper: simulate nsim datasets and return g_Go and g_NoGo vectors
  # ---------------------------------------------------------------------------
  .simulate_g <- function(mu_t_s, sigma_t_s, mu_c_s, sigma_c_s, seed_val) {
    set.seed(seed_val)

    Z_t         <- matrix(rnorm(nsim * n_t, mean = 0, sd = sigma_t_s), nrow = nsim)
    Z_t_rowmean <- rowSums(Z_t) / n_t
    bar_y_t_sim <- Z_t_rowmean + mu_t_s
    s_t_sim     <- sqrt(rowSums((Z_t - Z_t_rowmean) ^ 2) / (n_t - 1))

    if (design %in% c('controlled', 'external')) {
      Z_c         <- matrix(rnorm(nsim * n_c, mean = 0, sd = sigma_c_s), nrow = nsim)
      Z_c_rowmean <- rowSums(Z_c) / n_c
      bar_y_c_sim <- Z_c_rowmean + mu_c_s
      s_c_sim     <- sqrt(rowSums((Z_c - Z_c_rowmean) ^ 2) / (n_c - 1))
    } else {
      bar_y_c_sim <- NULL
      s_c_sim     <- NULL
    }

    theta0_go   <- if (prob == 'posterior') theta_TV  else theta_NULL
    theta0_nogo <- if (prob == 'posterior') theta_MAV else theta_NULL

    g_Go <- pbayespostpred1cont(
      prob = prob, design = design, prior = prior, CalcMethod = CalcMethod,
      theta0 = theta0_go, nMC = nMC,
      n_t = n_t, n_c = n_c, m_t = m_t, m_c = m_c,
      kappa0_t = kappa0_t, kappa0_c = kappa0_c,
      nu0_t = nu0_t, nu0_c = nu0_c,
      mu0_t = mu0_t, mu0_c = mu0_c,
      sigma0_t = sigma0_t, sigma0_c = sigma0_c,
      bar_y_t = bar_y_t_sim, bar_y_c = bar_y_c_sim,
      s_t = s_t_sim, s_c = s_c_sim,
      r = r,
      ne_t = ne_t, ne_c = ne_c,
      alpha0e_t = alpha0e_t, alpha0e_c = alpha0e_c,
      bar_ye_t = bar_ye_t, bar_ye_c = bar_ye_c,
      se_t = se_t, se_c = se_c,
      lower.tail = FALSE
    )

    g_NoGo <- pbayespostpred1cont(
      prob = prob, design = design, prior = prior, CalcMethod = CalcMethod,
      theta0 = theta0_nogo, nMC = nMC,
      n_t = n_t, n_c = n_c, m_t = m_t, m_c = m_c,
      kappa0_t = kappa0_t, kappa0_c = kappa0_c,
      nu0_t = nu0_t, nu0_c = nu0_c,
      mu0_t = mu0_t, mu0_c = mu0_c,
      sigma0_t = sigma0_t, sigma0_c = sigma0_c,
      bar_y_t = bar_y_t_sim, bar_y_c = bar_y_c_sim,
      s_t = s_t_sim, s_c = s_c_sim,
      r = r,
      ne_t = ne_t, ne_c = ne_c,
      alpha0e_t = alpha0e_t, alpha0e_c = alpha0e_c,
      bar_ye_t = bar_ye_t, bar_ye_c = bar_ye_c,
      se_t = se_t, se_c = se_c,
      lower.tail = TRUE
    )

    list(g_Go = g_Go, g_NoGo = g_NoGo)
  }

  # ---------------------------------------------------------------------------
  # Stage 1: Simulate nsim datasets for each calibration scenario independently
  # ---------------------------------------------------------------------------
  # Go-calibration scenario (typically Null): use seed
  sim_go   <- .simulate_g(mu_t_go,   sigma_t_go,   mu_c_go,   sigma_c_go,
                          seed_val = seed)
  # NoGo-calibration scenario (typically Alternative): use seed + 1
  sim_nogo <- .simulate_g(mu_t_nogo, sigma_t_nogo, mu_c_nogo, sigma_c_nogo,
                          seed_val = seed + 1L)

  gPost_Go   <- sim_go$g_Go
  gPost_NoGo <- sim_nogo$g_NoGo

  # ---------------------------------------------------------------------------
  # Stage 2: Sweep gamma_grid using vectorised outer comparison
  # ---------------------------------------------------------------------------
  # Marginal Go   criterion: g_Go   >= gamma  (lower.tail = FALSE -> P(theta > TV))
  # Marginal NoGo criterion: g_NoGo >= gamma  (lower.tail = TRUE  -> P(theta <= MAV))
  # Both Pr(Go) and Pr(NoGo) are monotone non-increasing in gamma.
  go_ind   <- outer(gPost_Go,   gamma_grid, ">=")  # nsim x n_gamma
  nogo_ind <- outer(gPost_NoGo, gamma_grid, ">=")  # nsim x n_gamma

  PrGo_grid   <- colMeans(go_ind)
  PrNoGo_grid <- colMeans(nogo_ind)

  # ---------------------------------------------------------------------------
  # Stage 3: Select optimal gamma values
  # ---------------------------------------------------------------------------
  # gamma_go  : smallest gamma with Pr(Go)   < target_go
  # gamma_nogo: smallest gamma with Pr(NoGo) < target_nogo
  go_mask   <- PrGo_grid   < target_go
  nogo_mask <- PrNoGo_grid < target_nogo

  idx_go <- which(go_mask)
  if (length(idx_go) == 0L) {
    gamma_go <- NA_real_
    PrGo_opt <- NA_real_
  } else {
    opt1     <- min(idx_go)
    gamma_go <- gamma_grid[opt1]
    PrGo_opt <- PrGo_grid[opt1]
  }

  idx_nogo <- which(nogo_mask)
  if (length(idx_nogo) == 0L) {
    gamma_nogo <- NA_real_
    PrNoGo_opt <- NA_real_
  } else {
    opt2       <- min(idx_nogo)
    gamma_nogo <- gamma_grid[opt2]
    PrNoGo_opt <- PrNoGo_grid[opt2]
  }

  # ---------------------------------------------------------------------------
  # Build and return result
  # ---------------------------------------------------------------------------
  result <- list(
    gamma_go     = gamma_go,
    gamma_nogo   = gamma_nogo,
    PrGo_opt     = PrGo_opt,
    PrNoGo_opt   = PrNoGo_opt,
    grid_results = data.frame(
      gamma_grid  = gamma_grid,
      PrGo_grid   = PrGo_grid,
      PrNoGo_grid = PrNoGo_grid
    )
  )

  class(result) <- 'getgamma1cont'
  return(result)
}
