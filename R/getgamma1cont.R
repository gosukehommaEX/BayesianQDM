#' Find Optimal Go/NoGo Thresholds for a Single Continuous Endpoint
#'
#' Computes the optimal Go threshold \eqn{\gamma_1} and NoGo threshold
#' \eqn{\gamma_2} for a single continuous endpoint by searching over a fine
#' grid of candidate values.  The search follows a two-stage approach:
#' posterior or predictive probabilities are precomputed once for every
#' Monte Carlo simulated dataset, and the operating characteristics under
#' each candidate \eqn{\gamma} are obtained by a fast indicator sweep without
#' additional probability evaluations.
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
#' @param theta.TV A numeric scalar giving the Target Value (TV) threshold
#'        for the treatment effect.  Required when \code{prob = 'posterior'};
#'        set to \code{NULL} otherwise.
#' @param theta.MAV A numeric scalar giving the Minimum Acceptable Value (MAV)
#'        threshold.  Must satisfy \code{theta.TV > theta.MAV}.  Required when
#'        \code{prob = 'posterior'}; set to \code{NULL} otherwise.
#' @param theta.NULL A numeric scalar giving the null hypothesis threshold
#'        for the predictive probability.  Required when
#'        \code{prob = 'predictive'}; set to \code{NULL} otherwise.
#' @param nMC A positive integer giving the number of inner Monte Carlo draws
#'        passed to \code{\link{pbayespostpred1cont}}.  Required when
#'        \code{CalcMethod = 'MC'}; set to \code{NULL} otherwise.
#' @param mu1 A numeric scalar giving the true mean for group 1 (treatment)
#'        for the scenario used to search for both \eqn{\gamma_1} and
#'        \eqn{\gamma_2}.
#' @param mu2 A numeric scalar giving the true mean for group 2 (control).
#'        Set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param sigma1 A positive numeric scalar giving the true standard deviation
#'        for group 1.
#' @param sigma2 A positive numeric scalar giving the true standard deviation
#'        for group 2.  Set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param target_go A numeric scalar in \code{(0, 1)} giving the target value
#'        for \eqn{\Pr(\mathrm{Go})} used to determine the optimal
#'        \eqn{\gamma_1}.  The comparison operator applied is specified by
#'        \code{crit_go}.
#' @param target_nogo A numeric scalar in \code{(0, 1)} giving the target
#'        value for \eqn{\Pr(\mathrm{NoGo})} used to determine the optimal
#'        \eqn{\gamma_2}.  The comparison operator applied is specified by
#'        \code{crit_nogo}.
#' @param crit_go A character string specifying the comparison operator
#'        applied to \eqn{\Pr(\mathrm{Go})} when searching for
#'        \eqn{\gamma_1}.  Must be one of \code{"<"}, \code{"<="}, \code{">"},
#'        or \code{">="}. Default is \code{"<"}.
#' @param crit_nogo A character string specifying the comparison operator
#'        applied to \eqn{\Pr(\mathrm{NoGo})} when searching for
#'        \eqn{\gamma_2}.  Must be one of \code{"<"}, \code{"<="}, \code{">"},
#'        or \code{">="}. Default is \code{"<"}.
#' @param sel_go A character string specifying whether to select the
#'        \code{"smallest"} or \code{"largest"} value in \code{gamma_grid}
#'        among those satisfying the \code{crit_go} criterion.
#'        Default is \code{"smallest"}.
#' @param sel_nogo A character string specifying whether to select the
#'        \code{"smallest"} or \code{"largest"} value in \code{gamma_grid}
#'        among those satisfying the \code{crit_nogo} criterion.
#'        Default is \code{"largest"}.
#' @param n1 A positive integer giving the number of patients in group 1
#'        (treatment) in the PoC trial.
#' @param n2 A positive integer giving the number of patients in group 2
#'        (control) in the PoC trial.  Required for \code{design = 'controlled'}
#'        or \code{'external'}; set to \code{NULL} for
#'        \code{design = 'uncontrolled'}.
#' @param m1 A positive integer giving the future sample size for group 1.
#'        Required when \code{prob = 'predictive'}; set to \code{NULL}
#'        otherwise.
#' @param m2 A positive integer giving the future sample size for group 2.
#'        Required when \code{prob = 'predictive'}; set to \code{NULL}
#'        otherwise.
#' @param kappa01 A positive numeric scalar giving the prior precision
#'        parameter for group 1 under the N-Inv-Chi-squared prior.  Required
#'        when \code{prior = 'N-Inv-Chisq'}; set to \code{NULL} otherwise.
#' @param kappa02 A positive numeric scalar giving the prior precision
#'        parameter for group 2.  Required when \code{prior = 'N-Inv-Chisq'}
#'        and \code{design != 'uncontrolled'}; set to \code{NULL} otherwise.
#' @param nu01 A positive numeric scalar giving the prior degrees of freedom
#'        for group 1 under the N-Inv-Chi-squared prior.  Required when
#'        \code{prior = 'N-Inv-Chisq'}; set to \code{NULL} otherwise.
#' @param nu02 A positive numeric scalar giving the prior degrees of freedom
#'        for group 2.  Required when \code{prior = 'N-Inv-Chisq'} and
#'        \code{design != 'uncontrolled'}; set to \code{NULL} otherwise.
#' @param mu01 A numeric scalar giving the prior mean for group 1 under the
#'        N-Inv-Chi-squared prior.  Required when \code{prior = 'N-Inv-Chisq'};
#'        set to \code{NULL} otherwise.
#' @param mu02 A numeric scalar giving the prior mean for group 2 (or the
#'        fixed hypothetical control mean for uncontrolled design).  Required
#'        when \code{prior = 'N-Inv-Chisq'}; set to \code{NULL} otherwise.
#' @param sigma01 A positive numeric scalar giving the prior scale parameter
#'        for group 1 under the N-Inv-Chi-squared prior.  Required when
#'        \code{prior = 'N-Inv-Chisq'}; set to \code{NULL} otherwise.
#' @param sigma02 A positive numeric scalar giving the prior scale parameter
#'        for group 2.  Required when \code{prior = 'N-Inv-Chisq'} and
#'        \code{design != 'uncontrolled'}; set to \code{NULL} otherwise.
#' @param r A positive numeric scalar giving the variance ratio
#'        \eqn{\sigma_2^2 / \sigma_1^2} used for uncontrolled design.
#'        Required when \code{design = 'uncontrolled'}; set to \code{NULL}
#'        otherwise.
#' @param ne1 A positive integer giving the number of patients in group 1
#'        of the external data set.  Required when \code{design = 'external'}
#'        and external treatment data are available; otherwise set to
#'        \code{NULL}.
#' @param ne2 A positive integer giving the number of patients in group 2
#'        of the external data set.  Required when \code{design = 'external'}
#'        and external control data are available; otherwise set to
#'        \code{NULL}.
#' @param alpha01 A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for external group 1 data.  Required when
#'        \code{design = 'external'} and \code{ne1} is non-NULL; set to
#'        \code{NULL} otherwise.
#' @param alpha02 A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for external group 2 data.  Required when
#'        \code{design = 'external'} and \code{ne2} is non-NULL; set to
#'        \code{NULL} otherwise.
#' @param bar.ye1 A numeric scalar giving the sample mean of the external
#'        group 1 data.  Required when \code{ne1} is non-NULL; set to
#'        \code{NULL} otherwise.
#' @param bar.ye2 A numeric scalar giving the sample mean of the external
#'        group 2 data.  Required when \code{ne2} is non-NULL; set to
#'        \code{NULL} otherwise.
#' @param se1 A positive numeric scalar giving the sample standard deviation
#'        of the external group 1 data.  Required when \code{ne1} is non-NULL;
#'        set to \code{NULL} otherwise.
#' @param se2 A positive numeric scalar giving the sample standard deviation
#'        of the external group 2 data.  Required when \code{ne2} is non-NULL;
#'        set to \code{NULL} otherwise.
#' @param gamma_grid A numeric vector of candidate threshold values in
#'        \code{(0, 1)} to search over.  Defaults to
#'        \code{seq(0.01, 0.99, by = 0.01)}.
#' @param seed A numeric scalar for reproducible random number generation.
#'
#' @return A list of class \code{getgamma1cont} with the following elements:
#' \describe{
#'   \item{gamma1}{Optimal Go threshold selected from \code{gamma_grid}
#'         according to \code{crit_go} and \code{sel_go}.
#'         \code{NA} if no value satisfies the criterion.}
#'   \item{gamma2}{Optimal NoGo threshold selected from \code{gamma_grid}
#'         according to \code{crit_nogo} and \code{sel_nogo}.
#'         \code{NA} if no value satisfies the criterion.}
#'   \item{PrGo_at_gamma1}{Pr(Go) evaluated at the optimal \eqn{\gamma_1}.
#'         \code{NA} if \code{gamma1} is \code{NA}.}
#'   \item{PrNoGo_at_gamma2}{Pr(NoGo) evaluated at the optimal \eqn{\gamma_2}.
#'         \code{NA} if \code{gamma2} is \code{NA}.}
#'   \item{gamma_grid}{The candidate grid used for the search.}
#'   \item{PrGo_grid}{Numeric vector of Pr(Go) values over \code{gamma_grid}.}
#'   \item{PrNoGo_grid}{Numeric vector of Pr(NoGo) values over
#'         \code{gamma_grid}.}
#' }
#'
#' @details
#' The function uses a two-stage simulate-then-sweep strategy:
#' \enumerate{
#'   \item \strong{Simulation}: \code{nsim} datasets are generated from the
#'         normal model with the specified true parameters.  Standardised
#'         residuals are generated once and shifted by \code{mu1} (and
#'         \code{mu2} for non-uncontrolled designs).  For each simulated
#'         dataset, \code{\link{pbayespostpred1cont}} is called twice to
#'         obtain the Go probability \eqn{g_{\mathrm{Go},i}} (at
#'         \code{theta.TV} or \code{theta.NULL}, \code{lower.tail = FALSE})
#'         and the NoGo probability \eqn{g_{\mathrm{NoGo},i}} (at
#'         \code{theta.MAV} or \code{theta.NULL}, \code{lower.tail = TRUE}).
#'         Both calls are fully vectorised over the \code{nsim} replicates.
#'         This step is independent of \eqn{\gamma}.
#'   \item \strong{Gamma sweep}: For each candidate \eqn{\gamma} in
#'         \code{gamma_grid}, Pr(Go) and Pr(NoGo) are estimated as the
#'         proportion of simulated datasets satisfying the respective
#'         decision criterion.  No further probability evaluations are
#'         required at this stage.
#' }
#' The optimal \eqn{\gamma_1} is the \code{sel_go} (\code{"smallest"} or
#' \code{"largest"}) value in \code{gamma_grid} for which
#' \eqn{\Pr(\mathrm{Go})} satisfies the \code{crit_go} comparison against
#' \code{target_go}.  Analogously, the optimal \eqn{\gamma_2} is the
#' \code{sel_nogo} value satisfying \code{crit_nogo} against
#' \code{target_nogo}.
#'
#' @examples
#' # Example 1: Controlled design, vague prior, posterior probability
#' # gamma1: smallest gamma such that Pr(Go) < 0.05
#' # gamma2: largest  gamma such that Pr(NoGo) < 0.20
#' getgamma1cont(
#'   nsim = 1000L, prob = 'posterior', design = 'controlled',
#'   prior = 'vague', CalcMethod = 'MM',
#'   theta.TV = 1.5, theta.MAV = 0.0, theta.NULL = NULL, nMC = NULL,
#'   mu1 = 1.0, mu2 = 1.0, sigma1 = 2.0, sigma2 = 2.0,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n1 = 15L, n2 = 15L, m1 = NULL, m2 = NULL,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   gamma_grid = seq(0.01, 0.99, by = 0.01), seed = 1L
#' )
#'
#' # Example 2: Uncontrolled design, N-Inv-Chisq prior, posterior probability
#' getgamma1cont(
#'   nsim = 1000L, prob = 'posterior', design = 'uncontrolled',
#'   prior = 'N-Inv-Chisq', CalcMethod = 'NI',
#'   theta.TV = 1.0, theta.MAV = 0.0, theta.NULL = NULL, nMC = NULL,
#'   mu1 = 1.5, mu2 = NULL, sigma1 = 1.5, sigma2 = NULL,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n1 = 20L, n2 = 20L, m1 = NULL, m2 = NULL,
#'   kappa01 = 2, kappa02 = 20, nu01 = 5, nu02 = 20,
#'   mu01 = 3.0, mu02 = 1.5, sigma01 = 1.5, sigma02 = 1.2,
#'   r = 1.0, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   gamma_grid = seq(0.01, 0.99, by = 0.01), seed = 2L
#' )
#'
#' # Example 3: Controlled design, vague prior, predictive probability
#' getgamma1cont(
#'   nsim = 1000L, prob = 'predictive', design = 'controlled',
#'   prior = 'vague', CalcMethod = 'MM',
#'   theta.TV = NULL, theta.MAV = NULL, theta.NULL = 1.0, nMC = NULL,
#'   mu1 = 1.0, mu2 = 1.0, sigma1 = 2.0, sigma2 = 2.0,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n1 = 15L, n2 = 15L, m1 = 50L, m2 = 50L,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   gamma_grid = seq(0.01, 0.99, by = 0.01), seed = 3L
#' )
#'
#' # Example 4: External design, vague prior, posterior probability
#' \dontrun{
#' getgamma1cont(
#'   nsim = 1000L, prob = 'posterior', design = 'external',
#'   prior = 'vague', CalcMethod = 'MM',
#'   theta.TV = 1.0, theta.MAV = 0.0, theta.NULL = NULL, nMC = NULL,
#'   mu1 = 1.0, mu2 = 1.0, sigma1 = 1.5, sigma2 = 1.5,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n1 = 12L, n2 = 12L, m1 = NULL, m2 = NULL,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   r = NULL, ne1 = NULL, ne2 = 20L, alpha01 = NULL, alpha02 = 0.5,
#'   bar.ye1 = NULL, bar.ye2 = 0.0, se1 = NULL, se2 = 1.5,
#'   gamma_grid = seq(0.01, 0.99, by = 0.01), seed = 4L
#' )
#' }
#'
#' @importFrom stats rnorm
#' @export
getgamma1cont <- function(nsim,
                          prob = 'posterior', design = 'controlled',
                          prior = 'vague', CalcMethod = 'NI',
                          theta.TV = NULL, theta.MAV = NULL, theta.NULL = NULL,
                          nMC = NULL,
                          mu1, mu2 = NULL, sigma1, sigma2 = NULL,
                          target_go, target_nogo,
                          crit_go  = '<',        crit_nogo  = '<',
                          sel_go   = 'smallest', sel_nogo   = 'largest',
                          n1, n2 = NULL,
                          m1 = NULL, m2 = NULL,
                          kappa01 = NULL, kappa02 = NULL,
                          nu01 = NULL, nu02 = NULL,
                          mu01 = NULL, mu02 = NULL,
                          sigma01 = NULL, sigma02 = NULL,
                          r = NULL,
                          ne1 = NULL, ne2 = NULL,
                          alpha01 = NULL, alpha02 = NULL,
                          bar.ye1 = NULL, bar.ye2 = NULL,
                          se1 = NULL, se2 = NULL,
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
    if (is.null(theta.TV) || is.null(theta.MAV)) {
      stop("'theta.TV' and 'theta.MAV' must be non-NULL when prob = 'posterior'")
    }
    if (!is.numeric(theta.TV) || length(theta.TV) != 1L || is.na(theta.TV)) {
      stop("'theta.TV' must be a single numeric value")
    }
    if (!is.numeric(theta.MAV) || length(theta.MAV) != 1L || is.na(theta.MAV)) {
      stop("'theta.MAV' must be a single numeric value")
    }
    if (theta.TV <= theta.MAV) {
      stop("'theta.TV' must be strictly greater than 'theta.MAV'")
    }
  } else {
    if (is.null(theta.NULL)) {
      stop("'theta.NULL' must be non-NULL when prob = 'predictive'")
    }
    if (!is.numeric(theta.NULL) || length(theta.NULL) != 1L ||
        is.na(theta.NULL)) {
      stop("'theta.NULL' must be a single numeric value")
    }
    if (is.null(m1) || is.null(m2)) {
      stop("'m1' and 'm2' must be non-NULL when prob = 'predictive'")
    }
    for (nm in c("m1", "m2")) {
      val <- get(nm)
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 1L) {
        stop(paste0("'", nm, "' must be a single positive integer"))
      }
    }
  }

  if (!is.numeric(mu1) || length(mu1) != 1L || is.na(mu1)) {
    stop("'mu1' must be a single numeric value")
  }

  if (!is.numeric(sigma1) || length(sigma1) != 1L || is.na(sigma1) ||
      sigma1 <= 0) {
    stop("'sigma1' must be a single positive numeric value")
  }

  if (design != 'uncontrolled') {
    if (is.null(mu2) || !is.numeric(mu2) || length(mu2) != 1L ||
        is.na(mu2)) {
      stop("'mu2' must be a single numeric value for controlled or external design")
    }
    if (is.null(sigma2) || !is.numeric(sigma2) || length(sigma2) != 1L ||
        is.na(sigma2) || sigma2 <= 0) {
      stop("'sigma2' must be a single positive numeric value for controlled or external design")
    }
  }

  for (nm in c("target_go", "target_nogo")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val <= 0 || val >= 1) {
      stop(paste0("'", nm, "' must be a single numeric value in (0, 1)"))
    }
  }

  valid_crit <- c('<', '<=', '>', '>=')
  if (!is.character(crit_go) || length(crit_go) != 1L ||
      !crit_go %in% valid_crit) {
    stop("'crit_go' must be one of '<', '<=', '>', '>='")
  }
  if (!is.character(crit_nogo) || length(crit_nogo) != 1L ||
      !crit_nogo %in% valid_crit) {
    stop("'crit_nogo' must be one of '<', '<=', '>', '>='")
  }

  valid_sel <- c('smallest', 'largest')
  if (!is.character(sel_go) || length(sel_go) != 1L ||
      !sel_go %in% valid_sel) {
    stop("'sel_go' must be either 'smallest' or 'largest'")
  }
  if (!is.character(sel_nogo) || length(sel_nogo) != 1L ||
      !sel_nogo %in% valid_sel) {
    stop("'sel_nogo' must be either 'smallest' or 'largest'")
  }

  if (!is.numeric(n1) || length(n1) != 1L || is.na(n1) ||
      n1 != floor(n1) || n1 < 1L) {
    stop("'n1' must be a single positive integer")
  }

  if (design != 'uncontrolled') {
    if (is.null(n2) || !is.numeric(n2) || length(n2) != 1L ||
        is.na(n2) || n2 != floor(n2) || n2 < 1L) {
      stop("'n2' must be a single positive integer for controlled or external design")
    }
  }

  if (prior == 'N-Inv-Chisq') {
    for (nm in c("kappa01", "nu01", "mu01", "sigma01")) {
      val <- get(nm)
      if (is.null(val)) {
        stop(paste0("'", nm, "' must be non-NULL when prior = 'N-Inv-Chisq'"))
      }
    }
    if (!is.numeric(kappa01) || length(kappa01) != 1L || is.na(kappa01) ||
        kappa01 <= 0) {
      stop("'kappa01' must be a single positive numeric value")
    }
    if (!is.numeric(nu01) || length(nu01) != 1L || is.na(nu01) ||
        nu01 <= 0) {
      stop("'nu01' must be a single positive numeric value")
    }
    if (!is.numeric(mu01) || length(mu01) != 1L || is.na(mu01)) {
      stop("'mu01' must be a single numeric value")
    }
    if (!is.numeric(sigma01) || length(sigma01) != 1L || is.na(sigma01) ||
        sigma01 <= 0) {
      stop("'sigma01' must be a single positive numeric value")
    }
    if (design != 'uncontrolled') {
      for (nm in c("kappa02", "nu02", "mu02", "sigma02")) {
        val <- get(nm)
        if (is.null(val)) {
          stop(paste0("'", nm, "' must be non-NULL when prior = 'N-Inv-Chisq' and design != 'uncontrolled'"))
        }
      }
      if (!is.numeric(kappa02) || length(kappa02) != 1L || is.na(kappa02) ||
          kappa02 <= 0) {
        stop("'kappa02' must be a single positive numeric value")
      }
      if (!is.numeric(nu02) || length(nu02) != 1L || is.na(nu02) ||
          nu02 <= 0) {
        stop("'nu02' must be a single positive numeric value")
      }
      if (!is.numeric(mu02) || length(mu02) != 1L || is.na(mu02)) {
        stop("'mu02' must be a single numeric value")
      }
      if (!is.numeric(sigma02) || length(sigma02) != 1L || is.na(sigma02) ||
          sigma02 <= 0) {
        stop("'sigma02' must be a single positive numeric value")
      }
    } else {
      if (is.null(mu02) || !is.numeric(mu02) || length(mu02) != 1L ||
          is.na(mu02)) {
        stop("'mu02' must be a single numeric value when design = 'uncontrolled'")
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
    has_ext1 <- !is.null(ne1) && !is.null(alpha01) && !is.null(bar.ye1) &&
      !is.null(se1)
    has_ext2 <- !is.null(ne2) && !is.null(alpha02) && !is.null(bar.ye2) &&
      !is.null(se2)
    if (!has_ext1 && !has_ext2) {
      stop("For design = 'external', at least one of (ne1, alpha01, bar.ye1, se1) or (ne2, alpha02, bar.ye2, se2) must be fully specified")
    }
    if (!is.null(ne1)) {
      if (!is.numeric(ne1) || length(ne1) != 1L || is.na(ne1) ||
          ne1 != floor(ne1) || ne1 < 1L) {
        stop("'ne1' must be a single positive integer")
      }
      if (!is.numeric(alpha01) || length(alpha01) != 1L || is.na(alpha01) ||
          alpha01 <= 0 || alpha01 > 1) {
        stop("'alpha01' must be a single numeric value in (0, 1]")
      }
      if (!is.numeric(bar.ye1) || length(bar.ye1) != 1L || is.na(bar.ye1)) {
        stop("'bar.ye1' must be a single numeric value")
      }
      if (!is.numeric(se1) || length(se1) != 1L || is.na(se1) || se1 <= 0) {
        stop("'se1' must be a single positive numeric value")
      }
    }
    if (!is.null(ne2)) {
      if (!is.numeric(ne2) || length(ne2) != 1L || is.na(ne2) ||
          ne2 != floor(ne2) || ne2 < 1L) {
        stop("'ne2' must be a single positive integer")
      }
      if (!is.numeric(alpha02) || length(alpha02) != 1L || is.na(alpha02) ||
          alpha02 <= 0 || alpha02 > 1) {
        stop("'alpha02' must be a single numeric value in (0, 1]")
      }
      if (!is.numeric(bar.ye2) || length(bar.ye2) != 1L || is.na(bar.ye2)) {
        stop("'bar.ye2' must be a single numeric value")
      }
      if (!is.numeric(se2) || length(se2) != 1L || is.na(se2) || se2 <= 0) {
        stop("'se2' must be a single positive numeric value")
      }
    }
  }

  if (!is.numeric(gamma_grid) || length(gamma_grid) < 1L ||
      any(is.na(gamma_grid)) || any(gamma_grid <= 0) || any(gamma_grid >= 1)) {
    stop("'gamma_grid' must be a numeric vector with all values in (0, 1)")
  }
  gamma_grid <- sort(unique(gamma_grid))
  n_gamma    <- length(gamma_grid)

  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("'seed' must be a single numeric value")
  }

  # ---------------------------------------------------------------------------
  # Stage 1: Simulate nsim datasets and precompute g(bar.y1, s1, bar.y2, s2)
  # ---------------------------------------------------------------------------
  set.seed(seed)

  Z1         <- matrix(rnorm(nsim * n1, mean = 0, sd = sigma1), nrow = nsim)
  Z1_rowmean <- rowSums(Z1) / n1
  bar.y1     <- Z1_rowmean + mu1
  s1_sim     <- sqrt(rowSums((Z1 - Z1_rowmean) ^ 2) / (n1 - 1))

  if (design %in% c('controlled', 'external')) {
    Z2         <- matrix(rnorm(nsim * n2, mean = 0, sd = sigma2), nrow = nsim)
    Z2_rowmean <- rowSums(Z2) / n2
    bar.y2_sim <- Z2_rowmean + mu2
    s2_sim     <- sqrt(rowSums((Z2 - Z2_rowmean) ^ 2) / (n2 - 1))
  } else {
    bar.y2_sim <- NULL
    s2_sim     <- NULL
  }

  theta0_go   <- if (prob == 'posterior') theta.TV  else theta.NULL
  theta0_nogo <- if (prob == 'posterior') theta.MAV else theta.NULL

  gPost_Go <- pbayespostpred1cont(
    prob = prob, design = design, prior = prior, CalcMethod = CalcMethod,
    theta0 = theta0_go, nMC = nMC,
    n1 = n1, n2 = n2, m1 = m1, m2 = m2,
    kappa01 = kappa01, kappa02 = kappa02,
    nu01 = nu01, nu02 = nu02,
    mu01 = mu01, mu02 = mu02,
    sigma01 = sigma01, sigma02 = sigma02,
    bar.y1 = bar.y1, bar.y2 = bar.y2_sim,
    s1 = s1_sim, s2 = s2_sim,
    r = r,
    ne1 = ne1, ne2 = ne2,
    alpha01 = alpha01, alpha02 = alpha02,
    bar.ye1 = bar.ye1, bar.ye2 = bar.ye2,
    se1 = se1, se2 = se2,
    lower.tail = FALSE
  )

  gPost_NoGo <- pbayespostpred1cont(
    prob = prob, design = design, prior = prior, CalcMethod = CalcMethod,
    theta0 = theta0_nogo, nMC = nMC,
    n1 = n1, n2 = n2, m1 = m1, m2 = m2,
    kappa01 = kappa01, kappa02 = kappa02,
    nu01 = nu01, nu02 = nu02,
    mu01 = mu01, mu02 = mu02,
    sigma01 = sigma01, sigma02 = sigma02,
    bar.y1 = bar.y1, bar.y2 = bar.y2_sim,
    s1 = s1_sim, s2 = s2_sim,
    r = r,
    ne1 = ne1, ne2 = ne2,
    alpha01 = alpha01, alpha02 = alpha02,
    bar.ye1 = bar.ye1, bar.ye2 = bar.ye2,
    se1 = se1, se2 = se2,
    lower.tail = TRUE
  )

  # ---------------------------------------------------------------------------
  # Stage 2: Sweep gamma_grid using vectorised outer comparison
  # ---------------------------------------------------------------------------
  # go_ind[i, k]   = 1 if gPost_Go[i]   >= gamma_grid[k]
  # nogo_ind[i, k] = 1 if gPost_NoGo[i] <= gamma_grid[k]
  go_ind   <- outer(gPost_Go,   gamma_grid, ">=")  # nsim x n_gamma
  nogo_ind <- outer(gPost_NoGo, gamma_grid, "<=")  # nsim x n_gamma

  PrGo_grid   <- colMeans(go_ind)
  PrNoGo_grid <- colMeans(nogo_ind)

  # ---------------------------------------------------------------------------
  # Stage 3: Apply comparison operators and select optimal gamma values
  # ---------------------------------------------------------------------------
  # Helper: apply comparison operator element-wise
  .compare <- function(x, op, val) {
    switch(op,
           '<'  = x <  val,
           '<=' = x <= val,
           '>'  = x >  val,
           '>=' = x >= val)
  }

  # Helper: select smallest or largest qualifying index
  .select_idx <- function(mask, sel) {
    idx <- which(mask)
    if (length(idx) == 0L) return(NA_integer_)
    if (sel == 'smallest') min(idx) else max(idx)
  }

  # gamma1
  opt1 <- .select_idx(.compare(PrGo_grid,   crit_go,   target_go),   sel_go)
  if (is.na(opt1)) {
    gamma1         <- NA_real_
    PrGo_at_gamma1 <- NA_real_
  } else {
    gamma1         <- gamma_grid[opt1]
    PrGo_at_gamma1 <- PrGo_grid[opt1]
  }

  # gamma2
  opt2 <- .select_idx(.compare(PrNoGo_grid, crit_nogo, target_nogo), sel_nogo)
  if (is.na(opt2)) {
    gamma2            <- NA_real_
    PrNoGo_at_gamma2  <- NA_real_
  } else {
    gamma2            <- gamma_grid[opt2]
    PrNoGo_at_gamma2  <- PrNoGo_grid[opt2]
  }

  # ---------------------------------------------------------------------------
  # Build and return result
  # ---------------------------------------------------------------------------
  result <- list(
    gamma1           = gamma1,
    gamma2           = gamma2,
    PrGo_at_gamma1   = PrGo_at_gamma1,
    PrNoGo_at_gamma2 = PrNoGo_at_gamma2,
    gamma_grid       = gamma_grid,
    PrGo_grid        = PrGo_grid,
    PrNoGo_grid      = PrNoGo_grid
  )

  class(result) <- 'getgamma1cont'
  return(result)
}
