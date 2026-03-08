#' Go/NoGo/Gray Decision Probabilities for Two Continuous Endpoints
#'
#' Computes the operating characteristics (Go, NoGo, Gray, and optionally Miss
#' probabilities) for a two-continuous-endpoint Bayesian Go/NoGo decision
#' framework by Monte Carlo simulation over treatment scenarios.
#'
#' @param nsim A positive integer. Number of simulated datasets per scenario.
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or
#'        \code{'external'}.
#' @param prior A character string specifying the prior distribution.
#'        Must be \code{'vague'} or \code{'N-Inv-Wishart'}.
#' @param GoRegions An integer vector specifying which of the nine posterior
#'        regions (R1--R9) or four predictive regions (R1--R4) constitute a
#'        Go decision.  For \code{prob = 'posterior'}, valid values are
#'        integers in 1--9; for \code{prob = 'predictive'}, in 1--4.
#'        A common choice is \code{GoRegions = 1} (both endpoints exceed TV
#'        or NULL for posterior/predictive, respectively).
#' @param NoGoRegions An integer vector specifying which regions constitute a
#'        NoGo decision.  A common choice is \code{NoGoRegions = 9} (both
#'        endpoints below MAV) for posterior, or \code{NoGoRegions = 4} for
#'        predictive.  Must be disjoint from \code{GoRegions}.
#' @param gamma_go A numeric scalar in \code{(0, 1)}. Go threshold:
#'        a Go decision is made if \eqn{P(\mathrm{GoRegions}) \ge \gamma_1}.
#' @param gamma_nogo A numeric scalar in \code{(0, 1)}. NoGo threshold:
#'        a NoGo decision is made if \eqn{P(\mathrm{NoGoRegions}) \ge \gamma_2}.
#'        No ordering constraint on \code{gamma_go} and \code{gamma_nogo} is
#'        imposed; their combination determines the frequency of Miss outcomes.
#' @param theta_TV1 A numeric scalar giving the target value (TV) threshold
#'        for Endpoint 1. Required when \code{prob = 'posterior'}; must
#'        satisfy \code{theta_TV1 > theta_MAV1}. Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta_MAV1 A numeric scalar giving the minimum acceptable value
#'        (MAV) threshold for Endpoint 1. Required when
#'        \code{prob = 'posterior'}; must satisfy
#'        \code{theta_TV1 > theta_MAV1}. Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta_TV2 A numeric scalar giving the target value (TV) threshold
#'        for Endpoint 2. Required when \code{prob = 'posterior'}; must
#'        satisfy \code{theta_TV2 > theta_MAV2}. Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta_MAV2 A numeric scalar giving the minimum acceptable value
#'        (MAV) threshold for Endpoint 2. Required when
#'        \code{prob = 'posterior'}; must satisfy
#'        \code{theta_TV2 > theta_MAV2}. Set to \code{NULL} when
#'        \code{prob = 'predictive'}.
#' @param theta_NULL1 A numeric scalar giving the null hypothesis threshold
#'        for Endpoint 1. Required when \code{prob = 'predictive'}; set to
#'        \code{NULL} when \code{prob = 'posterior'}.
#' @param theta_NULL2 A numeric scalar giving the null hypothesis threshold
#'        for Endpoint 2. Required when \code{prob = 'predictive'}; set to
#'        \code{NULL} when \code{prob = 'posterior'}.
#' @param n_t A positive integer giving the number of patients in the
#'        treatment group in the proof-of-concept (PoC) trial.
#' @param n_c A positive integer giving the number of patients in the
#'        control group in the PoC trial. For \code{design = 'uncontrolled'},
#'        this is the hypothetical control sample size (required for
#'        consistency with other designs).
#' @param m_t A positive integer giving the number of patients in the
#'        treatment group for the future trial. Required when
#'        \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param m_c A positive integer giving the number of patients in the
#'        control group for the future trial. Required when
#'        \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param mu_t Numeric matrix with 2 columns. Each row gives the true treatment
#'        mean vector for one scenario.  A length-2 vector is coerced to a
#'        1-row matrix.
#' @param Sigma_t A 2x2 positive-definite matrix. True treatment covariance.
#' @param mu_c Numeric matrix with 2 columns or a length-2 vector. True control
#'        (or hypothetical control) mean vector(s).
#' @param Sigma_c A 2x2 positive-definite matrix. True control covariance.
#' @param kappa0_t A positive numeric scalar giving the NIW prior
#'        concentration parameter for the treatment group. Required when
#'        \code{prior = 'N-Inv-Wishart'}; otherwise set to \code{NULL}.
#' @param nu0_t A numeric scalar giving the NIW prior degrees of freedom for
#'        the treatment group. Must be greater than 3. Required when
#'        \code{prior = 'N-Inv-Wishart'}; otherwise set to \code{NULL}.
#' @param mu0_t A numeric vector of length 2 giving the NIW prior mean for
#'        the treatment group. Required when
#'        \code{prior = 'N-Inv-Wishart'}; otherwise set to \code{NULL}.
#' @param Lambda0_t A 2x2 positive-definite numeric matrix giving the NIW
#'        prior scale matrix for the treatment group. Required when
#'        \code{prior = 'N-Inv-Wishart'}; otherwise set to \code{NULL}.
#' @param kappa0_c A positive numeric scalar giving the NIW prior
#'        concentration parameter for the control group. Required when
#'        \code{prior = 'N-Inv-Wishart'} and
#'        \code{design != 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param nu0_c A numeric scalar giving the NIW prior degrees of freedom for
#'        the control group. Must be greater than 3. Required when
#'        \code{prior = 'N-Inv-Wishart'} and
#'        \code{design != 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param mu0_c A numeric vector of length 2 giving the NIW prior mean for
#'        the control group, or the hypothetical control location when
#'        \code{design = 'uncontrolled'}. Required when
#'        \code{prior = 'N-Inv-Wishart'}; otherwise set to \code{NULL}.
#' @param Lambda0_c A 2x2 positive-definite numeric matrix giving the NIW
#'        prior scale matrix for the control group. Required when
#'        \code{prior = 'N-Inv-Wishart'} and
#'        \code{design != 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param r A positive numeric scalar giving the variance scaling factor for
#'        the hypothetical control distribution. Required when
#'        \code{design = 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param ne_t A positive integer giving the external treatment group sample
#'        size. Required when \code{design = 'external'} and external
#'        treatment data are used; otherwise set to \code{NULL}.
#' @param ne_c A positive integer giving the external control group sample
#'        size. Required when \code{design = 'external'} and external
#'        control data are used; otherwise set to \code{NULL}.
#' @param alpha0e_t A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for the external treatment data. Required when external
#'        treatment data are used; otherwise set to \code{NULL}.
#' @param alpha0e_c A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for the external control data. Required when external
#'        control data are used; otherwise set to \code{NULL}.
#' @param bar_ye_t A numeric vector of length 2 giving the external treatment
#'        group sample mean. Required when external treatment data are used;
#'        otherwise set to \code{NULL}.
#' @param bar_ye_c A numeric vector of length 2 giving the external control
#'        group sample mean. Required when external control data are used;
#'        otherwise set to \code{NULL}.
#' @param se_t A 2x2 numeric matrix giving the external treatment group
#'        sum-of-squares matrix. Required when external treatment data are
#'        used; otherwise set to \code{NULL}.
#' @param se_c A 2x2 numeric matrix giving the external control group
#'        sum-of-squares matrix. Required when external control data are
#'        used; otherwise set to \code{NULL}.
#' @param nMC A positive integer giving the number of Monte Carlo draws used
#'        to estimate region probabilities. Default is \code{10000}. Required
#'        when \code{CalcMethod = 'MC'}. May be set to \code{NULL} when
#'        \code{CalcMethod = 'MM'} and \eqn{\nu_k > 4} (the MM method uses
#'        \code{mvtnorm::pmvt} analytically); if \code{CalcMethod = 'MM'} but
#'        \eqn{\nu_k \le 4} causes a fallback to MC, \code{nMC} must be a
#'        positive integer.
#' @param CalcMethod A character string specifying the computation method.
#'        Must be \code{'MC'} (Monte Carlo, default) or \code{'MM'}
#'        (Moment-Matching via \code{mvtnorm::pmvt}). When
#'        \code{CalcMethod = 'MM'} and \eqn{\nu_k \le 4}, a warning is issued
#'        and the function falls back to \code{CalcMethod = 'MC'}.
#' @param error_if_Miss Logical. If \code{TRUE} (default), the function stops
#'        with an error when positive Miss probability is obtained. If
#'        \code{FALSE}, Miss probability is handled according to
#'        \code{Gray_inc_Miss}.
#' @param Gray_inc_Miss Logical. If \code{TRUE}, Miss probability is included
#'        in Gray probability. If \code{FALSE} (default), Miss probability is
#'        reported as a separate column. Active only when
#'        \code{error_if_Miss = FALSE}.
#' @param seed A single numeric value. Seed for reproducible random number
#'        generation.
#'
#' @return A data frame of class \code{c("pbayesdecisionprob2cont", "data.frame")}
#'         with columns for the scenario parameters and Go, NoGo, Gray
#'         (and optionally Miss) probabilities.  All input parameters are
#'         attached as attributes.
#'
#' @details
#' The function follows the same structure as
#' \code{\link{pbayesdecisionprob1cont}}:
#' \enumerate{
#'   \item For each scenario \eqn{s}, \code{nsim} datasets are simulated by
#'         generating treatment (and control) observations from
#'         \eqn{N_2(\mu_k^{(s)}, \Sigma_k)}.  To minimise overhead, raw
#'         standardised residuals are generated \emph{once} (scenario-
#'         invariant) and shifted by the scenario mean.
#'   \item All \code{nsim} simulated sufficient statistics
#'         \eqn{(\bar{y}_{1,i}, S_{1,i})} (and \eqn{(\bar{y}_{2,i}, S_{2,i})}
#'         for controlled/external designs) are passed to
#'         \code{\link{pbayespostpred2cont}} in a \emph{single vectorised
#'         call}, returning an \eqn{nsim \times n_{\rm regions}} matrix of
#'         region probabilities.
#'   \item Go/NoGo/Miss probabilities are obtained as the column means of
#'         indicator matrices derived from the region probability matrix.
#' }
#'
#' @examples
#' # Example 1: Controlled design, posterior probability, vague prior
#' Sigma <- matrix(c(4.0, 0.8, 0.8, 1.0), 2, 2)
#' pbayesdecisionprob2cont(
#'   nsim = 100L, prob = 'posterior', design = 'controlled',
#'   prior = 'vague',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   gamma_go = 0.8, gamma_nogo = 0.8,
#'   theta_TV1 = 1.5, theta_MAV1 = 0.5,
#'   theta_TV2 = 1.0, theta_MAV2 = 0.3,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   n_t = 20L, n_c = 20L, m_t = NULL, m_c = NULL,
#'   mu_t = rbind(c(1.0, 0.5), c(2.5, 1.5), c(4.0, 2.5)),
#'   Sigma_t = Sigma,
#'   mu_c = rbind(c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)),
#'   Sigma_c = Sigma,
#'   kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
#'   kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
#'   r = NULL,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   nMC = 500L, CalcMethod = 'MC',
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 1L
#' )
#'
#' # Example 2: Uncontrolled design, posterior probability, NIW prior
#' Sigma <- matrix(c(4.0, 0.8, 0.8, 1.0), 2, 2)
#' L0    <- matrix(c(8.0, 0.0, 0.0, 2.0), 2, 2)
#' pbayesdecisionprob2cont(
#'   nsim = 100L, prob = 'posterior', design = 'uncontrolled',
#'   prior = 'N-Inv-Wishart',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   gamma_go = 0.8, gamma_nogo = 0.8,
#'   theta_TV1 = 1.5, theta_MAV1 = 0.5,
#'   theta_TV2 = 1.0, theta_MAV2 = 0.3,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   n_t = 20L, n_c = NULL, m_t = NULL, m_c = NULL,
#'   mu_t = rbind(c(1.0, 0.5), c(2.5, 1.5), c(4.0, 2.5)),
#'   Sigma_t = Sigma,
#'   mu_c = rbind(c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)),
#'   Sigma_c = Sigma,
#'   kappa0_t = 2.0, nu0_t = 5.0, mu0_t = c(2.0, 1.0), Lambda0_t = L0,
#'   kappa0_c = NULL, nu0_c = NULL, mu0_c = c(0.0, 0.0), Lambda0_c = NULL,
#'   r = 1.0,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   nMC = 500L, CalcMethod = 'MC',
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 3L
#' )
#'
#' # Example 3: External design (control only), posterior probability, NIW prior
#' Sigma  <- matrix(c(4.0, 0.8, 0.8, 1.0), 2, 2)
#' L0     <- matrix(c(8.0, 0.0, 0.0, 2.0), 2, 2)
#' se_mat <- matrix(c(7.0, 1.2, 1.2, 1.8), 2, 2)
#' pbayesdecisionprob2cont(
#'   nsim = 100L, prob = 'posterior', design = 'external',
#'   prior = 'N-Inv-Wishart',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   gamma_go = 0.8, gamma_nogo = 0.8,
#'   theta_TV1 = 1.5, theta_MAV1 = 0.5,
#'   theta_TV2 = 1.0, theta_MAV2 = 0.3,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   n_t = 20L, n_c = 20L, m_t = NULL, m_c = NULL,
#'   mu_t = rbind(c(1.0, 0.5), c(2.5, 1.5), c(4.0, 2.5)),
#'   Sigma_t = Sigma,
#'   mu_c = rbind(c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)),
#'   Sigma_c = Sigma,
#'   kappa0_t = 2.0, nu0_t = 5.0, mu0_t = c(2.0, 1.0), Lambda0_t = L0,
#'   kappa0_c = 2.0, nu0_c = 5.0, mu0_c = c(0.0, 0.0), Lambda0_c = L0,
#'   r = NULL,
#'   ne_t = NULL, ne_c = 15L, alpha0e_t = NULL, alpha0e_c = 0.5,
#'   bar_ye_t = NULL, bar_ye_c = c(0.2, 0.1), se_t = NULL, se_c = se_mat,
#'   nMC = 500L, CalcMethod = 'MC',
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 5L
#' )
#'
#' # Example 4: Controlled design, predictive probability, NIW prior
#' Sigma <- matrix(c(4.0, 0.8, 0.8, 1.0), 2, 2)
#' L0    <- matrix(c(8.0, 0.0, 0.0, 2.0), 2, 2)
#' pbayesdecisionprob2cont(
#'   nsim = 100L, prob = 'predictive', design = 'controlled',
#'   prior = 'N-Inv-Wishart',
#'   GoRegions = 1L, NoGoRegions = 4L,
#'   gamma_go = 0.8, gamma_nogo = 0.8,
#'   theta_TV1 = NULL, theta_MAV1 = NULL,
#'   theta_TV2 = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 0.5, theta_NULL2 = 0.3,
#'   n_t = 20L, n_c = 20L, m_t = 60L, m_c = 60L,
#'   mu_t = rbind(c(1.0, 0.5), c(2.5, 1.5), c(4.0, 2.5)),
#'   Sigma_t = Sigma,
#'   mu_c = rbind(c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)),
#'   Sigma_c = Sigma,
#'   kappa0_t = 2.0, nu0_t = 5.0, mu0_t = c(2.0, 1.0), Lambda0_t = L0,
#'   kappa0_c = 2.0, nu0_c = 5.0, mu0_c = c(0.0, 0.0), Lambda0_c = L0,
#'   r = NULL,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   nMC = 500L, CalcMethod = 'MC',
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 4L
#' )
#'
#' # Example 5: Uncontrolled design, predictive probability, vague prior
#' Sigma <- matrix(c(4.0, 0.8, 0.8, 1.0), 2, 2)
#' pbayesdecisionprob2cont(
#'   nsim = 100L, prob = 'predictive', design = 'uncontrolled',
#'   prior = 'vague',
#'   GoRegions = 1L, NoGoRegions = 4L,
#'   gamma_go = 0.8, gamma_nogo = 0.8,
#'   theta_TV1 = NULL, theta_MAV1 = NULL,
#'   theta_TV2 = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 0.5, theta_NULL2 = 0.3,
#'   n_t = 20L, n_c = NULL, m_t = 60L, m_c = 60L,
#'   mu_t = rbind(c(1.0, 0.5), c(2.5, 1.5), c(4.0, 2.5)),
#'   Sigma_t = Sigma,
#'   mu_c = NULL,
#'   Sigma_c = NULL,
#'   kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
#'   kappa0_c = NULL, nu0_c = NULL, mu0_c = c(0.0, 0.0), Lambda0_c = NULL,
#'   r = 1.0,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   nMC = 500L, CalcMethod = 'MC',
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 8L
#' )
#'
#' # Example 6: External design (control only), predictive probability, NIW prior
#' Sigma  <- matrix(c(4.0, 0.8, 0.8, 1.0), 2, 2)
#' L0     <- matrix(c(8.0, 0.0, 0.0, 2.0), 2, 2)
#' se_mat <- matrix(c(7.0, 1.2, 1.2, 1.8), 2, 2)
#' pbayesdecisionprob2cont(
#'   nsim = 100L, prob = 'predictive', design = 'external',
#'   prior = 'N-Inv-Wishart',
#'   GoRegions = 1L, NoGoRegions = 4L,
#'   gamma_go = 0.8, gamma_nogo = 0.8,
#'   theta_TV1 = NULL, theta_MAV1 = NULL,
#'   theta_TV2 = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 0.5, theta_NULL2 = 0.3,
#'   n_t = 20L, n_c = 20L, m_t = 60L, m_c = 60L,
#'   mu_t = rbind(c(1.0, 0.5), c(2.5, 1.5), c(4.0, 2.5)),
#'   Sigma_t = Sigma,
#'   mu_c = rbind(c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)),
#'   Sigma_c = Sigma,
#'   kappa0_t = 2.0, nu0_t = 5.0, mu0_t = c(2.0, 1.0), Lambda0_t = L0,
#'   kappa0_c = 2.0, nu0_c = 5.0, mu0_c = c(0.0, 0.0), Lambda0_c = L0,
#'   r = NULL,
#'   ne_t = NULL, ne_c = 15L, alpha0e_t = NULL, alpha0e_c = 0.5,
#'   bar_ye_t = NULL, bar_ye_c = c(0.2, 0.1), se_t = NULL, se_c = se_mat,
#'   nMC = 500L, CalcMethod = 'MC',
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 9L
#' )
#'
#' @importFrom stats rnorm
#' @export
pbayesdecisionprob2cont <- function(nsim,
                                    prob,
                                    design,
                                    prior,
                                    GoRegions, NoGoRegions,
                                    gamma_go, gamma_nogo,
                                    theta_TV1     = NULL, theta_MAV1  = NULL,
                                    theta_TV2     = NULL, theta_MAV2  = NULL,
                                    theta_NULL1   = NULL, theta_NULL2 = NULL,
                                    n_t, n_c      = NULL,
                                    m_t           = NULL, m_c         = NULL,
                                    mu_t, Sigma_t,
                                    mu_c = NULL, Sigma_c = NULL,
                                    kappa0_t      = NULL, nu0_t       = NULL,
                                    mu0_t         = NULL, Lambda0_t   = NULL,
                                    kappa0_c      = NULL, nu0_c       = NULL,
                                    mu0_c         = NULL, Lambda0_c   = NULL,
                                    r             = NULL,
                                    ne_t          = NULL, ne_c        = NULL,
                                    alpha0e_t     = NULL, alpha0e_c   = NULL,
                                    bar_ye_t      = NULL, bar_ye_c    = NULL,
                                    se_t          = NULL, se_c        = NULL,
                                    nMC           = NULL,
                                    CalcMethod    = 'MC',
                                    error_if_Miss = TRUE,
                                    Gray_inc_Miss = FALSE,
                                    seed) {

  # ---------------------------------------------------------------------------
  # Section 1: Input validation
  # ---------------------------------------------------------------------------

  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) ||
      nsim != floor(nsim) || nsim < 1L)
    stop("'nsim' must be a single positive integer")
  nsim <- as.integer(nsim)

  if (!is.character(prob) || length(prob) != 1L ||
      !prob %in% c('posterior', 'predictive'))
    stop("'prob' must be either 'posterior' or 'predictive'")

  if (!is.character(design) || length(design) != 1L ||
      !design %in% c('controlled', 'uncontrolled', 'external'))
    stop("'design' must be 'controlled', 'uncontrolled', or 'external'")

  if (!is.character(prior) || length(prior) != 1L ||
      !prior %in% c('vague', 'N-Inv-Wishart'))
    stop("'prior' must be 'vague' or 'N-Inv-Wishart'")

  max_region <- if (prob == 'posterior') 9L else 4L
  if (!is.numeric(GoRegions) || length(GoRegions) < 1L ||
      any(is.na(GoRegions)) || any(GoRegions != floor(GoRegions)) ||
      any(GoRegions < 1L) || any(GoRegions > max_region))
    stop(sprintf("'GoRegions' must be integer(s) in 1:%d", max_region))
  if (!is.numeric(NoGoRegions) || length(NoGoRegions) < 1L ||
      any(is.na(NoGoRegions)) || any(NoGoRegions != floor(NoGoRegions)) ||
      any(NoGoRegions < 1L) || any(NoGoRegions > max_region))
    stop(sprintf("'NoGoRegions' must be integer(s) in 1:%d", max_region))
  if (any(GoRegions %in% NoGoRegions))
    stop("'GoRegions' and 'NoGoRegions' must be disjoint")

  for (nm in c('gamma_go', 'gamma_nogo')) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val <= 0 || val >= 1)
      stop(paste0("'", nm, "' must be a single numeric value in (0, 1)"))
  }

  if (prob == 'posterior') {
    for (nm in c('theta_TV1', 'theta_MAV1', 'theta_TV2', 'theta_MAV2')) {
      val <- get(nm)
      if (is.null(val))
        stop(paste0("'", nm, "' must be non-NULL when prob = 'posterior'"))
      if (!is.numeric(val) || length(val) != 1L || is.na(val))
        stop(paste0("'", nm, "' must be a single numeric value"))
    }
    if (theta_TV1 <= theta_MAV1)
      stop("'theta_TV1' must be strictly greater than 'theta_MAV1'")
    if (theta_TV2 <= theta_MAV2)
      stop("'theta_TV2' must be strictly greater than 'theta_MAV2'")
  } else {
    for (nm in c('theta_NULL1', 'theta_NULL2')) {
      val <- get(nm)
      if (is.null(val))
        stop(paste0("'", nm, "' must be non-NULL when prob = 'predictive'"))
      if (!is.numeric(val) || length(val) != 1L || is.na(val))
        stop(paste0("'", nm, "' must be a single numeric value"))
    }
    if (is.null(m_t) || is.null(m_c))
      stop("'m_t' and 'm_c' must be non-NULL when prob = 'predictive'")
    for (nm in c('m_t', 'm_c')) {
      val <- get(nm)
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 1L)
        stop(paste0("'", nm, "' must be a single positive integer"))
    }
  }

  if (!is.numeric(n_t) || length(n_t) != 1L || is.na(n_t) ||
      n_t != floor(n_t) || n_t < 1L)
    stop("'n_t' must be a single positive integer")
  n_t <- as.integer(n_t)

  if (design %in% c('controlled', 'external')) {
    if (is.null(n_c))
      stop("'n_c' must be non-NULL when design is 'controlled' or 'external'")
    if (!is.numeric(n_c) || length(n_c) != 1L || is.na(n_c) ||
        n_c != floor(n_c) || n_c < 1L)
      stop("'n_c' must be a single positive integer")
    n_c <- as.integer(n_c)
  }

  # Coerce mu_t to matrix if needed
  if (is.vector(mu_t) && length(mu_t) == 2L) mu_t <- matrix(mu_t, nrow = 1L)
  if (!is.matrix(mu_t) || ncol(mu_t) != 2L)
    stop("'mu_t' must be a matrix with 2 columns (or a length-2 vector)")

  n_scen_t <- nrow(mu_t)

  if (design == 'uncontrolled') {
    # mu_c is not used in simulation; n_scen is determined by mu_t alone
    n_scen <- n_scen_t
  } else {
    if (is.vector(mu_c) && length(mu_c) == 2L) mu_c <- matrix(mu_c, nrow = 1L)
    if (!is.matrix(mu_c) || ncol(mu_c) != 2L)
      stop("'mu_c' must be a matrix with 2 columns (or a length-2 vector)")
    n_scen_c <- nrow(mu_c)
    if (n_scen_c != n_scen_t)
      stop("'mu_t' and 'mu_c' must have the same number of rows")
    n_scen <- n_scen_t
  }

  if (!is.matrix(Sigma_t) || nrow(Sigma_t) != 2L || ncol(Sigma_t) != 2L)
    stop("'Sigma_t' must be a 2x2 numeric matrix")
  if (design != 'uncontrolled') {
    if (!is.matrix(Sigma_c) || nrow(Sigma_c) != 2L || ncol(Sigma_c) != 2L)
      stop("'Sigma_c' must be a 2x2 numeric matrix")
  }

  if (prior == 'N-Inv-Wishart') {
    for (nm in c('kappa0_t', 'nu0_t')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L ||
          is.na(val) || val <= 0)
        stop(paste0("'", nm, "' must be a single positive numeric for NIW prior"))
    }
    if (is.null(mu0_t) || !is.numeric(mu0_t) || length(mu0_t) != 2L)
      stop("'mu0_t' must be a length-2 numeric vector for NIW prior")
    if (is.null(Lambda0_t) || !is.matrix(Lambda0_t) || nrow(Lambda0_t) != 2L)
      stop("'Lambda0_t' must be a 2x2 numeric matrix for NIW prior")
    if (design %in% c('controlled', 'external')) {
      for (nm in c('kappa0_c', 'nu0_c')) {
        val <- get(nm)
        if (is.null(val) || !is.numeric(val) || length(val) != 1L ||
            is.na(val) || val <= 0)
          stop(paste0("'", nm, "' must be a single positive numeric for NIW prior"))
      }
      if (is.null(mu0_c) || !is.numeric(mu0_c) || length(mu0_c) != 2L)
        stop("'mu0_c' must be a length-2 numeric vector for NIW prior")
      if (is.null(Lambda0_c) || !is.matrix(Lambda0_c) || nrow(Lambda0_c) != 2L)
        stop("'Lambda0_c' must be a 2x2 numeric matrix for NIW prior")
    }
  }

  if (design == 'uncontrolled') {
    if (is.null(r) || !is.numeric(r) || length(r) != 1L || is.na(r) || r <= 0)
      stop("'r' must be a single positive numeric when design = 'uncontrolled'")
    if (is.null(mu0_c) || !is.numeric(mu0_c) || length(mu0_c) != 2L)
      stop("'mu0_c' must be a length-2 numeric vector when design = 'uncontrolled'")
  }

  if (design == 'external') {
    has_ext_t <- !is.null(ne_t) && !is.null(alpha0e_t) &&
      !is.null(bar_ye_t) && !is.null(se_t)
    has_ext_c <- !is.null(ne_c) && !is.null(alpha0e_c) &&
      !is.null(bar_ye_c) && !is.null(se_c)
    if (!has_ext_t && !has_ext_c)
      stop(paste0("For design = 'external', at least one complete set of ",
                  "external data must be provided"))
  }

  if (!is.character(CalcMethod) || length(CalcMethod) != 1L ||
      !CalcMethod %in% c('MC', 'MM'))
    stop("'CalcMethod' must be either 'MC' or 'MM'")

  # nMC validation: required for CalcMethod = 'MC', optional for CalcMethod = 'MM'
  if (CalcMethod == 'MC') {
    if (is.null(nMC))
      stop("'nMC' must be non-NULL when CalcMethod = 'MC'")
    if (!is.numeric(nMC) || length(nMC) != 1L || is.na(nMC) ||
        nMC != floor(nMC) || nMC < 1L)
      stop("'nMC' must be a single positive integer")
    nMC <- as.integer(nMC)
  } else {
    # CalcMethod == 'MM': nMC may be NULL or a positive integer
    if (!is.null(nMC)) {
      if (!is.numeric(nMC) || length(nMC) != 1L || is.na(nMC) ||
          nMC != floor(nMC) || nMC < 1L)
        stop("'nMC' must be a single positive integer or NULL")
      nMC <- as.integer(nMC)
    }
  }

  if (!is.logical(error_if_Miss) || length(error_if_Miss) != 1L ||
      is.na(error_if_Miss))
    stop("'error_if_Miss' must be a single logical value")

  if (!is.logical(Gray_inc_Miss) || length(Gray_inc_Miss) != 1L ||
      is.na(Gray_inc_Miss))
    stop("'Gray_inc_Miss' must be a single logical value")

  # ---------------------------------------------------------------------------
  # Section 2: Pre-generate simulation data (scenario-invariant)
  #
  # Raw standardised residuals are generated once.
  # For group k, the sample mean for scenario s and replicate i is:
  #   ybar_k[s, i, ] = (Z_k_colsums[i, ] / n_k) + mu_k[s, ]
  # The scatter matrix S_k[i] depends only on Z_k (not on mu_k[s,]).
  # ---------------------------------------------------------------------------

  set.seed(seed)

  R_Sigma_t <- chol(Sigma_t)
  Z_t_raw   <- matrix(rnorm(nsim * n_t * 2L),
                      nrow = nsim * n_t, ncol = 2L) %*% R_Sigma_t
  block_t   <- rep(seq_len(nsim), each = n_t)

  Z_t_colsums      <- apply(Z_t_raw, 2L, function(col) tapply(col, block_t, sum))
  Z_t_colmeans_rep <- Z_t_colsums[block_t, ] / n_t
  Z_t_centered     <- Z_t_raw - Z_t_colmeans_rep

  # S_t stored as 3 unique elements per replicate (symmetric 2x2)
  S_t_11 <- tapply(Z_t_centered[, 1L] ^ 2,                block_t, sum)
  S_t_12 <- tapply(Z_t_centered[, 1L] * Z_t_centered[, 2L], block_t, sum)
  S_t_22 <- tapply(Z_t_centered[, 2L] ^ 2,                block_t, sum)

  if (design %in% c('controlled', 'external')) {
    R_Sigma_c <- chol(Sigma_c)
    Z_c_raw   <- matrix(rnorm(nsim * n_c * 2L),
                        nrow = nsim * n_c, ncol = 2L) %*% R_Sigma_c
    block_c   <- rep(seq_len(nsim), each = n_c)

    Z_c_colsums      <- apply(Z_c_raw, 2L, function(col) tapply(col, block_c, sum))
    Z_c_colmeans_rep <- Z_c_colsums[block_c, ] / n_c
    Z_c_centered     <- Z_c_raw - Z_c_colmeans_rep

    S_c_11 <- tapply(Z_c_centered[, 1L] ^ 2,                block_c, sum)
    S_c_12 <- tapply(Z_c_centered[, 1L] * Z_c_centered[, 2L], block_c, sum)
    S_c_22 <- tapply(Z_c_centered[, 2L] ^ 2,                block_c, sum)
  }

  # ---------------------------------------------------------------------------
  # Section 3: Scenario loop
  #
  # For each scenario s:
  #   1. Construct nsim sufficient statistics by shifting the pre-generated
  #      residuals by the scenario mean.
  #   2. Call pbayespostpred2cont once in vectorised mode (ybar_t as matrix,
  #      S_t as list) to obtain an nsim x n_regions matrix of Pr_R values.
  #   3. Compute PrGo and PrNoGo for each replicate and classify.
  # ---------------------------------------------------------------------------

  result_mat <- matrix(0, nrow = n_scen, ncol = 3L)

  for (s in seq_len(n_scen)) {

    # Shift residual column sums by scenario mean (vectorised over nsim)
    ybar_t_sim <- sweep(Z_t_colsums / n_t, 2L, mu_t[s, ], '+')
    # ybar_t_sim: nsim x 2 matrix

    # Build S_t list for this scenario (scenario-invariant values)
    S_t_list <- vector('list', nsim)
    for (i in seq_len(nsim)) {
      S_t_list[[i]] <- matrix(c(S_t_11[i], S_t_12[i], S_t_12[i], S_t_22[i]),
                              nrow = 2L, ncol = 2L)
    }

    if (design %in% c('controlled', 'external')) {
      ybar_c_sim <- sweep(Z_c_colsums / n_c, 2L, mu_c[s, ], '+')
      S_c_list <- vector('list', nsim)
      for (i in seq_len(nsim)) {
        S_c_list[[i]] <- matrix(c(S_c_11[i], S_c_12[i], S_c_12[i], S_c_22[i]),
                                nrow = 2L, ncol = 2L)
      }
    } else {
      ybar_c_sim <- NULL
      S_c_list   <- NULL
    }

    # Vectorised call: returns nsim x n_regions matrix
    Pr_R_mat <- pbayespostpred2cont(
      prob        = prob,
      design      = design,
      prior       = prior,
      theta_TV1   = theta_TV1,   theta_MAV1  = theta_MAV1,
      theta_TV2   = theta_TV2,   theta_MAV2  = theta_MAV2,
      theta_NULL1 = theta_NULL1, theta_NULL2 = theta_NULL2,
      n_t = n_t, n_c = n_c,
      ybar_t = ybar_t_sim, S_t = S_t_list,
      ybar_c = ybar_c_sim, S_c = S_c_list,
      m_t = m_t, m_c = m_c,
      kappa0_t  = kappa0_t,  nu0_t     = nu0_t,
      mu0_t     = mu0_t,     Lambda0_t = Lambda0_t,
      kappa0_c  = kappa0_c,  nu0_c     = nu0_c,
      mu0_c     = mu0_c,     Lambda0_c = Lambda0_c,
      r        = r,
      ne_t = ne_t, ne_c = ne_c,
      alpha0e_t = alpha0e_t, alpha0e_c = alpha0e_c,
      bar_ye_t  = bar_ye_t,  bar_ye_c  = bar_ye_c,
      se_t      = se_t,      se_c      = se_c,
      nMC    = nMC,
      CalcMethod = CalcMethod
    )
    # Pr_R_mat: nsim x n_regions (column names R1...R9 or R1...R4)

    # Go/NoGo probability per replicate
    PrGo_vec   <- rowSums(Pr_R_mat[, GoRegions,   drop = FALSE])
    PrNoGo_vec <- rowSums(Pr_R_mat[, NoGoRegions, drop = FALSE])

    # --- Decision indicators (Go, NoGo, Miss are mutually exclusive; Gray is the complement) ---
    ind_Go   <- (PrGo_vec >= gamma_go) & (PrNoGo_vec <  gamma_nogo)
    ind_NoGo <- (PrGo_vec <  gamma_go) & (PrNoGo_vec >= gamma_nogo)
    ind_Miss <- (PrGo_vec >= gamma_go) & (PrNoGo_vec >= gamma_nogo)

    result_mat[s, 1L] <- mean(ind_Go)
    result_mat[s, 2L] <- mean(ind_NoGo)
    result_mat[s, 3L] <- mean(ind_Miss)
  }

  # ---------------------------------------------------------------------------
  # Section 4: Assemble output
  # ---------------------------------------------------------------------------

  result_mat[result_mat < .Machine$double.eps ^ 0.25] <- 0

  if (error_if_Miss && any(result_mat[, 3L] > 0))
    stop("Positive Miss probability detected.  Please re-consider the chosen thresholds.")

  GrayProb <- if (Gray_inc_Miss) {
    1 - result_mat[, 1L] - result_mat[, 2L]
  } else {
    1 - rowSums(result_mat)
  }

  if (design == 'uncontrolled') {
    results <- data.frame(
      mu_t1 = mu_t[, 1L],
      mu_t2 = mu_t[, 2L],
      Go      = result_mat[, 1L],
      Gray    = GrayProb,
      NoGo    = result_mat[, 2L]
    )
  } else {
    results <- data.frame(
      mu_t1 = mu_t[, 1L],
      mu_t2 = mu_t[, 2L],
      mu_c1 = mu_c[, 1L],
      mu_c2 = mu_c[, 2L],
      Go      = result_mat[, 1L],
      Gray    = GrayProb,
      NoGo    = result_mat[, 2L]
    )
  }

  if (!error_if_Miss && !Gray_inc_Miss)
    results$Miss <- result_mat[, 3L]

  # Suppress floating-point noise in output probability columns
  prob_cols <- c('Go', 'Gray', 'NoGo', 'Miss')
  prob_cols <- prob_cols[prob_cols %in% names(results)]
  results[prob_cols] <- lapply(results[prob_cols], function(col) {
    ifelse(col < .Machine$double.eps ^ 0.25, 0, col)
  })

  # Attach metadata as attributes
  attr(results, 'prob')           <- prob
  attr(results, 'design')         <- design
  attr(results, 'prior')          <- prior
  attr(results, 'nsim')           <- nsim
  attr(results, 'nMC')            <- nMC
  attr(results, 'CalcMethod')     <- CalcMethod
  attr(results, 'GoRegions')      <- GoRegions
  attr(results, 'NoGoRegions')    <- NoGoRegions
  attr(results, 'gamma_go')       <- gamma_go
  attr(results, 'gamma_nogo')     <- gamma_nogo
  attr(results, 'theta_TV1')      <- theta_TV1
  attr(results, 'theta_MAV1')     <- theta_MAV1
  attr(results, 'theta_TV2')      <- theta_TV2
  attr(results, 'theta_MAV2')     <- theta_MAV2
  attr(results, 'theta_NULL1')    <- theta_NULL1
  attr(results, 'theta_NULL2')    <- theta_NULL2
  attr(results, 'n_t')            <- n_t
  attr(results, 'n_c')            <- n_c
  attr(results, 'm_t')            <- m_t
  attr(results, 'm_c')            <- m_c
  attr(results, 'Sigma_t')        <- Sigma_t
  attr(results, 'Sigma_c')        <- Sigma_c
  attr(results, 'kappa0_t')       <- kappa0_t
  attr(results, 'nu0_t')          <- nu0_t
  attr(results, 'mu0_t')          <- mu0_t
  attr(results, 'Lambda0_t')      <- Lambda0_t
  attr(results, 'kappa0_c')       <- kappa0_c
  attr(results, 'nu0_c')          <- nu0_c
  attr(results, 'mu0_c')          <- mu0_c
  attr(results, 'Lambda0_c')      <- Lambda0_c
  attr(results, 'r')              <- r
  attr(results, 'ne_t')           <- ne_t
  attr(results, 'ne_c')           <- ne_c
  attr(results, 'alpha0e_t')       <- alpha0e_t
  attr(results, 'alpha0e_c')       <- alpha0e_c
  attr(results, 'bar_ye_t')       <- bar_ye_t
  attr(results, 'bar_ye_c')       <- bar_ye_c
  attr(results, 'se_t')           <- se_t
  attr(results, 'se_c')           <- se_c
  attr(results, 'error_if_Miss')  <- error_if_Miss
  attr(results, 'Gray_inc_Miss')  <- Gray_inc_Miss
  attr(results, 'seed')           <- seed

  class(results) <- c('pbayesdecisionprob2cont', 'data.frame')

  return(results)
}

#' Print Method for pbayesdecisionprob2cont Objects
#'
#' Displays a formatted summary of Go/NoGo/Gray decision probabilities for
#' two-continuous-endpoint results returned by
#' \code{\link{pbayesdecisionprob2cont}}.
#'
#' @param x An object of class \code{pbayesdecisionprob2cont}.
#' @param digits A positive integer specifying the number of decimal places
#'        for probability values. Default is 4.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.pbayesdecisionprob2cont <- function(x, digits = 4, ...) {

  # Helper: format a scalar/vector/NULL as string
  fmt <- function(v) {
    if (is.null(v))    return("NULL")
    if (length(v) > 1) return(paste0("(", paste(v, collapse = ", "), ")"))
    as.character(v)
  }

  # Helper: format a 2x2 matrix as "[r1c1, r1c2; r2c1, r2c2]"
  fmt_mat <- function(m) {
    if (is.null(m)) return("NULL")
    sprintf("[%s, %s; %s, %s]", m[1,1], m[1,2], m[2,1], m[2,2])
  }

  # Extract metadata from attributes
  prob          <- attr(x, "prob")
  design        <- attr(x, "design")
  prior         <- attr(x, "prior")
  nsim          <- attr(x, "nsim")
  nMC           <- attr(x, "nMC")
  CalcMethod    <- attr(x, "CalcMethod")
  GoRegions     <- attr(x, "GoRegions")
  NoGoRegions   <- attr(x, "NoGoRegions")
  gamma_go      <- attr(x, "gamma_go")
  gamma_nogo    <- attr(x, "gamma_nogo")
  theta_TV1     <- attr(x, "theta_TV1")
  theta_MAV1    <- attr(x, "theta_MAV1")
  theta_TV2     <- attr(x, "theta_TV2")
  theta_MAV2    <- attr(x, "theta_MAV2")
  theta_NULL1   <- attr(x, "theta_NULL1")
  theta_NULL2   <- attr(x, "theta_NULL2")
  n_t           <- attr(x, "n_t")
  n_c           <- attr(x, "n_c")
  m_t           <- attr(x, "m_t")
  m_c           <- attr(x, "m_c")
  Sigma_t       <- attr(x, "Sigma_t")
  Sigma_c       <- attr(x, "Sigma_c")
  kappa0_t      <- attr(x, "kappa0_t")
  nu0_t         <- attr(x, "nu0_t")
  mu0_t         <- attr(x, "mu0_t")
  Lambda0_t     <- attr(x, "Lambda0_t")
  kappa0_c      <- attr(x, "kappa0_c")
  nu0_c         <- attr(x, "nu0_c")
  mu0_c         <- attr(x, "mu0_c")
  Lambda0_c     <- attr(x, "Lambda0_c")
  r             <- attr(x, "r")
  ne_t          <- attr(x, "ne_t")
  ne_c          <- attr(x, "ne_c")
  alpha0e_t     <- attr(x, "alpha0e_t")
  alpha0e_c     <- attr(x, "alpha0e_c")
  bar_ye_t      <- attr(x, "bar_ye_t")
  bar_ye_c      <- attr(x, "bar_ye_c")
  se_t          <- attr(x, "se_t")
  se_c          <- attr(x, "se_c")
  error_if_Miss <- attr(x, "error_if_Miss")
  Gray_inc_Miss <- attr(x, "Gray_inc_Miss")
  seed          <- attr(x, "seed")

  # Build info lines with fixed label width (lw) for consistent alignment
  lw  <- 17L   # label field width
  pad <- "  "  # left margin

  lines <- character(0)
  lines <- c(lines, sprintf("%s%-*s: %s",        pad, lw, "Probability type", prob))
  lines <- c(lines, sprintf("%s%-*s: %s",        pad, lw, "Design",           design))
  lines <- c(lines, sprintf("%s%-*s: %s",        pad, lw, "Prior",            prior))
  lines <- c(lines, sprintf("%s%-*s: %s",        pad, lw, "Method",           fmt(CalcMethod)))
  lines <- c(lines, sprintf("%s%-*s: nsim = %s", pad, lw, "Simulations",      fmt(nsim)))
  lines <- c(lines, sprintf("%s%-*s: nMC = %s",  pad, lw, "MC draws",         fmt(nMC)))
  lines <- c(lines, sprintf("%s%-*s: %s",        pad, lw, "Seed",             fmt(seed)))

  # Threshold(s): split posterior across two lines to avoid overflow
  if (prob == "posterior") {
    lines <- c(lines, sprintf("%s%-*s: TV1 = %s, MAV1 = %s",
                              pad, lw, "Threshold(s)",
                              fmt(theta_TV1), fmt(theta_MAV1)))
    lines <- c(lines, sprintf("%s%-*s  TV2 = %s, MAV2 = %s",
                              pad, lw, "",
                              fmt(theta_TV2), fmt(theta_MAV2)))
  } else {
    lines <- c(lines, sprintf("%s%-*s: NULL1 = %s, NULL2 = %s",
                              pad, lw, "Threshold(s)",
                              fmt(theta_NULL1), fmt(theta_NULL2)))
  }

  lines <- c(lines, sprintf("%s%-*s: gamma_go = %s",
                            pad, lw, "Go  threshold",  fmt(gamma_go)))
  lines <- c(lines, sprintf("%s%-*s: gamma_nogo = %s",
                            pad, lw, "NoGo threshold", fmt(gamma_nogo)))
  lines <- c(lines, sprintf("%s%-*s: {%s}",
                            pad, lw, "Go  regions",
                            paste(GoRegions,   collapse = ", ")))
  lines <- c(lines, sprintf("%s%-*s: {%s}",
                            pad, lw, "NoGo regions",
                            paste(NoGoRegions, collapse = ", ")))
  lines <- c(lines, sprintf("%s%-*s: n_t = %s, n_c = %s",
                            pad, lw, "Sample size", fmt(n_t), fmt(n_c)))
  lines <- c(lines, sprintf("%s%-*s: Sigma_t = %s",
                            pad, lw, "True cov (treat.)", fmt_mat(Sigma_t)))
  if (design != "uncontrolled") {
    lines <- c(lines, sprintf("%s%-*s: Sigma_c = %s",
                              pad, lw, "True cov (cont.) ", fmt_mat(Sigma_c)))
  }

  # N-Inv-Wishart prior: treatment and control each split across two lines
  if (prior == "N-Inv-Wishart") {
    lines <- c(lines, sprintf(
      "%s%-*s: kappa0_t = %s, nu0_t = %s, mu0_t = %s",
      pad, lw, "Prior (treatment)", fmt(kappa0_t), fmt(nu0_t), fmt(mu0_t)))
    lines <- c(lines, sprintf("%s%-*s  Lambda0_t = %s",
                              pad, lw, "", fmt_mat(Lambda0_t)))
    if (design %in% c("controlled", "external")) {
      lines <- c(lines, sprintf(
        "%s%-*s: kappa0_c = %s, nu0_c = %s, mu0_c = %s",
        pad, lw, "Prior (control)  ", fmt(kappa0_c), fmt(nu0_c), fmt(mu0_c)))
      lines <- c(lines, sprintf("%s%-*s  Lambda0_c = %s",
                                pad, lw, "", fmt_mat(Lambda0_c)))
    }
  }

  if (design == "uncontrolled") {
    lines <- c(lines, sprintf("%s%-*s: mu0_c = %s, r = %s",
                              pad, lw, "Hyp. control", fmt(mu0_c), fmt(r)))
  }
  if (prob == "predictive") {
    lines <- c(lines, sprintf("%s%-*s: m_t = %s, m_c = %s",
                              pad, lw, "Future trial", fmt(m_t), fmt(m_c)))
  }
  if (design == "external") {
    # External data: treatment and control on separate lines
    lines <- c(lines, sprintf(
      "%s%-*s: ne_t = %s, alpha0e_t = %s",
      pad, lw, "External (treat.)", fmt(ne_t), fmt(alpha0e_t)))
    if (!is.null(bar_ye_t)) {
      lines <- c(lines, sprintf("%s%-*s  bar_ye_t = %s, se_t = %s",
                                pad, lw, "", fmt(bar_ye_t), fmt(se_t)))
    }
    lines <- c(lines, sprintf(
      "%s%-*s: ne_c = %s, alpha0e_c = %s",
      pad, lw, "External (cont.) ", fmt(ne_c), fmt(alpha0e_c)))
    if (!is.null(bar_ye_c)) {
      lines <- c(lines, sprintf("%s%-*s  bar_ye_c = %s, se_c = %s",
                                pad, lw, "", fmt(bar_ye_c), fmt(se_c)))
    }
  }
  lines <- c(lines, sprintf("%s%-*s: error_if_Miss = %s, Gray_inc_Miss = %s",
                            pad, lw, "Miss handling",
                            fmt(error_if_Miss), fmt(Gray_inc_Miss)))

  # Determine separator width dynamically from the longest line
  title     <- "Go/NoGo/Gray Decision Probabilities (Two Continuous Endpoints)"
  sep_width <- max(nchar(title), max(nchar(lines)))
  sep       <- strrep("-", sep_width)

  # Print header block
  cat(title, "\n")
  cat(sep, "\n")
  for (ln in lines) cat(ln, "\n")
  cat(sep, "\n")

  # Format probability columns only (not scenario columns)
  scenario_cols <- c("mu_t1", "mu_t2", "mu_c1", "mu_c2")
  prob_cols     <- names(x)[!names(x) %in% scenario_cols]

  x_print <- x
  x_print[prob_cols] <- lapply(x[prob_cols], function(col) {
    formatC(col, digits = digits, format = "f")
  })

  # Print table without row names (explicit call to avoid recursion)
  print.data.frame(x_print, row.names = FALSE, quote = FALSE)
  cat(sep, "\n")

  invisible(x)
}
