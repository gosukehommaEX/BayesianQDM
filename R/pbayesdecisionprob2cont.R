#' Go/NoGo/Gray Decision Probabilities for a Clinical Trial with Two
#' Continuous Endpoints
#'
#' Evaluates operating characteristics (Go, NoGo, Gray probabilities) for
#' clinical trials with two continuous endpoints under the Bayesian framework.
#' The function supports controlled, uncontrolled, and external-control designs,
#' and uses both posterior probability and posterior predictive probability
#' criteria.
#'
#' For each of the \code{nsim} simulated PoC datasets, sufficient statistics
#' (\eqn{\bar{y}_k}, \eqn{S_k}) are computed and passed to
#' \code{\link{pbayespostpred2cont}}, which returns the nine (posterior) or four
#' (predictive) region probabilities. Go and NoGo probabilities for that
#' replicate are obtained by summing the specified \code{GoRegions} and
#' \code{NoGoRegions}. Operating characteristics are then the proportions
#' of replicates classified as Go, NoGo, or Gray across all \code{nsim}
#' simulations.
#'
#' @param nsim A positive integer giving the number of Monte Carlo simulations
#'        for evaluating operating characteristics.
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or
#'        \code{'external'}.
#' @param prior A character string specifying the prior distribution.
#'        Must be \code{'vague'} or \code{'N-Inv-Wishart'}.
#' @param GoRegions An integer vector specifying which of the nine posterior
#'        regions (R1--R9) or four predictive regions (R1--R4) constitute a
#'        Go decision. For \code{prob = 'posterior'}, valid values are
#'        integers in 1--9; for \code{prob = 'predictive'}, in 1--4.
#'        A common choice is \code{GoRegions = 1} (both endpoints exceed TV
#'        or NULL).
#' @param NoGoRegions An integer vector specifying which regions constitute
#'        a NoGo decision. Must be disjoint from \code{GoRegions}.
#' @param gamma1 A numeric scalar in \code{(0, 1)} giving the minimum
#'        posterior/predictive probability required for a Go decision.
#'        A Go decision is made when the sum of region probabilities for
#'        \code{GoRegions} is at least \code{gamma1}.
#' @param gamma2 A numeric scalar in \code{(0, 1)} giving the minimum
#'        posterior/predictive probability required for a NoGo decision.
#'        Unlike single-endpoint designs, \code{gamma2} may be greater than,
#'        equal to, or less than \code{gamma1}.
#' @param theta.TV1 A numeric scalar giving the TV threshold for Endpoint 1.
#'        Required when \code{prob = 'posterior'}; otherwise \code{NULL}.
#' @param theta.MAV1 A numeric scalar giving the MAV threshold for Endpoint 1.
#'        Required when \code{prob = 'posterior'}; must satisfy
#'        \code{theta.TV1 > theta.MAV1}. Otherwise \code{NULL}.
#' @param theta.TV2 A numeric scalar giving the TV threshold for Endpoint 2.
#'        Required when \code{prob = 'posterior'}; otherwise \code{NULL}.
#' @param theta.MAV2 A numeric scalar giving the MAV threshold for Endpoint 2.
#'        Required when \code{prob = 'posterior'}; must satisfy
#'        \code{theta.TV2 > theta.MAV2}. Otherwise \code{NULL}.
#' @param theta.NULL1 A numeric scalar giving the null hypothesis threshold
#'        for Endpoint 1. Required when \code{prob = 'predictive'};
#'        otherwise \code{NULL}.
#' @param theta.NULL2 A numeric scalar giving the null hypothesis threshold
#'        for Endpoint 2. Required when \code{prob = 'predictive'};
#'        otherwise \code{NULL}.
#' @param n1 A positive integer giving the PoC sample size for group 1
#'        (treatment). Must be >= 3 when \code{prior = 'vague'}.
#' @param n2 A positive integer giving the PoC sample size for group 2
#'        (control). Required for \code{design = 'controlled'} or
#'        \code{'external'}; set to \code{NULL} for \code{'uncontrolled'}.
#'        Must be >= 3 when \code{prior = 'vague'}.
#' @param m1 A positive integer giving the future sample size for group 1.
#'        Required when \code{prob = 'predictive'}; otherwise \code{NULL}.
#' @param m2 A positive integer giving the future sample size for group 2.
#'        Required when \code{prob = 'predictive'}; otherwise \code{NULL}.
#' @param mu1 A numeric vector of length 2 or a matrix with 2 columns giving
#'        the true mean(s) for group 1 (treatment). Multiple rows define
#'        multiple scenarios. Each row is \eqn{(\mu_{t1}, \mu_{t2})}.
#' @param Sigma1 A 2x2 positive-definite numeric matrix giving the true
#'        covariance matrix for group 1.
#' @param mu2 A numeric vector of length 2 or a matrix with 2 columns giving
#'        the true mean(s) for group 2 (control). For
#'        \code{design = 'uncontrolled'}, this parameter is not used in
#'        data generation but is retained in the output for reference.
#'        Must have the same number of rows as \code{mu1}.
#' @param Sigma2 A 2x2 positive-definite numeric matrix giving the true
#'        covariance matrix for group 2. For \code{design = 'uncontrolled'},
#'        not used in data generation.
#' @param kappa01 A positive numeric scalar. Prior precision parameter for
#'        group 1 when \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param nu01 A positive numeric scalar (> 2). Prior degrees of freedom for
#'        group 1 when \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param mu01 A numeric vector of length 2. Prior mean for group 1 when
#'        \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param Lambda01 A 2x2 positive-definite numeric matrix. Prior scale matrix
#'        for group 1 when \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param kappa02 A positive numeric scalar. Prior precision parameter for
#'        group 2 when \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param nu02 A positive numeric scalar (> 2). Prior degrees of freedom for
#'        group 2 when \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param mu02 A numeric vector of length 2. Prior mean for group 2 when
#'        \code{prior = 'N-Inv-Wishart'}, or the hypothetical control location
#'        when \code{design = 'uncontrolled'}; otherwise \code{NULL}.
#' @param Lambda02 A 2x2 positive-definite numeric matrix. Prior scale matrix
#'        for group 2 when \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param r A positive numeric scalar. Variance scaling factor for the
#'        hypothetical control. Required when \code{design = 'uncontrolled'};
#'        otherwise \code{NULL}.
#' @param ne1 A positive integer giving the external sample size for group 1.
#'        Required when \code{design = 'external'} and external treatment
#'        data are used; otherwise \code{NULL}.
#' @param ne2 A positive integer giving the external sample size for group 2.
#'        Required when \code{design = 'external'} and external control data
#'        are used; otherwise \code{NULL}.
#' @param alpha01e A numeric scalar in \code{(0, 1]}. Power prior weight for
#'        group 1 external data. Required when \code{ne1} is non-NULL;
#'        otherwise \code{NULL}.
#' @param alpha02e A numeric scalar in \code{(0, 1]}. Power prior weight for
#'        group 2 external data. Required when \code{ne2} is non-NULL;
#'        otherwise \code{NULL}.
#' @param ybar_e1 A numeric vector of length 2. Sample mean of external data
#'        for group 1. Required when \code{ne1} is non-NULL; otherwise
#'        \code{NULL}.
#' @param ybar_e2 A numeric vector of length 2. Sample mean of external data
#'        for group 2. Required when \code{ne2} is non-NULL; otherwise
#'        \code{NULL}.
#' @param Se1 A 2x2 positive-definite numeric matrix. Sum-of-squares matrix
#'        of external data for group 1. Required when \code{ne1} is non-NULL;
#'        otherwise \code{NULL}.
#' @param Se2 A 2x2 positive-definite numeric matrix. Sum-of-squares matrix
#'        of external data for group 2. Required when \code{ne2} is non-NULL;
#'        otherwise \code{NULL}.
#' @param nMC A positive integer giving the number of Monte Carlo samples
#'        passed to \code{\link{pbayespostpred2cont}} for each simulated dataset.
#'        Used when \code{method = 'MC'} or when \code{method = 'MM'} falls back
#'        to MC.  Default is \code{1000L}.
#' @param method A character string specifying the computation method passed to
#'        \code{\link{pbayespostpred2cont}}.  Must be \code{'MC'} (Monte Carlo;
#'        default) or \code{'MM'} (Moment-Matching approximation).  See
#'        \code{\link{pbayespostpred2cont}} for details on the MM approximation
#'        and its fallback behaviour.
#' @param error_if_Miss A logical scalar; if \code{TRUE} (default), the
#'        function stops with an error if Miss probability is positive.
#' @param Gray_inc_Miss A logical scalar; if \code{TRUE}, Miss probability
#'        is added to Gray. If \code{FALSE} (default), Miss is reported
#'        separately. Active only when \code{error_if_Miss = FALSE}.
#' @param seed A numeric scalar for reproducible random number generation.
#'
#' @return A data frame with one row per scenario and columns:
#' \describe{
#'   \item{mu1_ep1}{True treatment mean for Endpoint 1.}
#'   \item{mu1_ep2}{True treatment mean for Endpoint 2.}
#'   \item{mu2_ep1}{True control mean for Endpoint 1 (omitted for
#'         \code{design = 'uncontrolled'}).}
#'   \item{mu2_ep2}{True control mean for Endpoint 2 (omitted for
#'         \code{design = 'uncontrolled'}).}
#'   \item{Go}{Probability of a Go decision.}
#'   \item{Gray}{Probability of a Gray decision.}
#'   \item{NoGo}{Probability of a NoGo decision.}
#'   \item{Miss}{Probability of a Miss decision (only present when
#'         \code{error_if_Miss = FALSE} and \code{Gray_inc_Miss = FALSE}).}
#' }
#' The object has class \code{c('pbayesdecisionprob2cont', 'data.frame')} and
#' carries metadata as attributes for use in the print method.
#'
#' @details
#' For each of the \code{nsim} simulated PoC datasets:
#' \enumerate{
#'   \item Bivariate Normal data are generated for each arm using the
#'         specified true means and covariance matrices.
#'   \item Sufficient statistics (\eqn{\bar{y}_k}, \eqn{S_k}) are computed.
#'   \item \code{\link{pbayespostpred2cont}} is called to obtain region
#'         probabilities \eqn{(R_1, \ldots, R_9)} or \eqn{(R_1, \ldots, R_4)}.
#'   \item \eqn{PrGo = \sum_{l \in \text{GoRegions}} R_l} and
#'         \eqn{PrNoGo = \sum_{l \in \text{NoGoRegions}} R_l} are evaluated
#'         against \code{gamma1} and \code{gamma2}.
#' }
#' Operating characteristics are estimated as proportions across \code{nsim}
#' replicates.
#'
#' @examples
#' # Example 1: Controlled design, posterior probability, vague prior
#' Sigma <- matrix(c(4.0, 0.8, 0.8, 1.0), 2, 2)
#' pbayesdecisionprob2cont(
#'   nsim = 100L, prob = 'posterior', design = 'controlled',
#'   prior = 'vague',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   gamma1 = 0.8, gamma2 = 0.8,
#'   theta.TV1 = 1.5, theta.MAV1 = 0.5,
#'   theta.TV2 = 1.0, theta.MAV2 = 0.3,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   n1 = 20L, n2 = 20L, m1 = NULL, m2 = NULL,
#'   mu1 = rbind(c(1.0, 0.5), c(2.5, 1.5), c(4.0, 2.5)),
#'   Sigma1 = Sigma,
#'   mu2 = rbind(c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)),
#'   Sigma2 = Sigma,
#'   kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
#'   kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
#'   r = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
#'   ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
#'   nMC = 500L, error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 1L
#' )
#'
#' # Example 2: Controlled design, posterior probability, NIW prior
#' Sigma <- matrix(c(4.0, 0.8, 0.8, 1.0), 2, 2)
#' L0    <- matrix(c(8.0, 0.0, 0.0, 2.0), 2, 2)
#' pbayesdecisionprob2cont(
#'   nsim = 100L, prob = 'posterior', design = 'controlled',
#'   prior = 'N-Inv-Wishart',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   gamma1 = 0.8, gamma2 = 0.8,
#'   theta.TV1 = 1.5, theta.MAV1 = 0.5,
#'   theta.TV2 = 1.0, theta.MAV2 = 0.3,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   n1 = 20L, n2 = 20L, m1 = NULL, m2 = NULL,
#'   mu1 = rbind(c(1.0, 0.5), c(2.5, 1.5), c(4.0, 2.5)),
#'   Sigma1 = Sigma,
#'   mu2 = rbind(c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)),
#'   Sigma2 = Sigma,
#'   kappa01 = 2.0, nu01 = 5.0, mu01 = c(2.0, 1.0), Lambda01 = L0,
#'   kappa02 = 2.0, nu02 = 5.0, mu02 = c(0.0, 0.0), Lambda02 = L0,
#'   r = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
#'   ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
#'   nMC = 500L, error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 2L
#' )
#'
#' # Example 3: Uncontrolled design, posterior probability, NIW prior
#' Sigma <- matrix(c(4.0, 0.8, 0.8, 1.0), 2, 2)
#' L0    <- matrix(c(8.0, 0.0, 0.0, 2.0), 2, 2)
#' pbayesdecisionprob2cont(
#'   nsim = 100L, prob = 'posterior', design = 'uncontrolled',
#'   prior = 'N-Inv-Wishart',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   gamma1 = 0.8, gamma2 = 0.8,
#'   theta.TV1 = 1.5, theta.MAV1 = 0.5,
#'   theta.TV2 = 1.0, theta.MAV2 = 0.3,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   n1 = 20L, n2 = NULL, m1 = NULL, m2 = NULL,
#'   mu1 = rbind(c(1.0, 0.5), c(2.5, 1.5), c(4.0, 2.5)),
#'   Sigma1 = Sigma,
#'   mu2 = rbind(c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)),
#'   Sigma2 = Sigma,
#'   kappa01 = 2.0, nu01 = 5.0, mu01 = c(2.0, 1.0), Lambda01 = L0,
#'   kappa02 = NULL, nu02 = NULL, mu02 = c(0.0, 0.0), Lambda02 = NULL,
#'   r = 1.0,
#'   ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
#'   ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
#'   nMC = 500L, error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 3L
#' )
#'
#' # Example 4: Controlled design, predictive probability, NIW prior
#' Sigma <- matrix(c(4.0, 0.8, 0.8, 1.0), 2, 2)
#' L0    <- matrix(c(8.0, 0.0, 0.0, 2.0), 2, 2)
#' pbayesdecisionprob2cont(
#'   nsim = 100L, prob = 'predictive', design = 'controlled',
#'   prior = 'N-Inv-Wishart',
#'   GoRegions = 1L, NoGoRegions = 4L,
#'   gamma1 = 0.8, gamma2 = 0.8,
#'   theta.TV1 = NULL, theta.MAV1 = NULL,
#'   theta.TV2 = NULL, theta.MAV2 = NULL,
#'   theta.NULL1 = 0.5, theta.NULL2 = 0.3,
#'   n1 = 20L, n2 = 20L, m1 = 60L, m2 = 60L,
#'   mu1 = rbind(c(1.0, 0.5), c(2.5, 1.5), c(4.0, 2.5)),
#'   Sigma1 = Sigma,
#'   mu2 = rbind(c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)),
#'   Sigma2 = Sigma,
#'   kappa01 = 2.0, nu01 = 5.0, mu01 = c(2.0, 1.0), Lambda01 = L0,
#'   kappa02 = 2.0, nu02 = 5.0, mu02 = c(0.0, 0.0), Lambda02 = L0,
#'   r = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
#'   ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
#'   nMC = 500L, error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 4L
#' )
#'
#' # Example 5: External design (control only), posterior probability, NIW prior
#' \dontrun{
#' Sigma <- matrix(c(4.0, 0.8, 0.8, 1.0), 2, 2)
#' L0    <- matrix(c(8.0, 0.0, 0.0, 2.0), 2, 2)
#' Se2   <- matrix(c(7.0, 1.2, 1.2, 1.8), 2, 2)
#' pbayesdecisionprob2cont(
#'   nsim = 100L, prob = 'posterior', design = 'external',
#'   prior = 'N-Inv-Wishart',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   gamma1 = 0.8, gamma2 = 0.8,
#'   theta.TV1 = 1.5, theta.MAV1 = 0.5,
#'   theta.TV2 = 1.0, theta.MAV2 = 0.3,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   n1 = 20L, n2 = 20L, m1 = NULL, m2 = NULL,
#'   mu1 = rbind(c(1.0, 0.5), c(2.5, 1.5), c(4.0, 2.5)),
#'   Sigma1 = Sigma,
#'   mu2 = rbind(c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)),
#'   Sigma2 = Sigma,
#'   kappa01 = 2.0, nu01 = 5.0, mu01 = c(2.0, 1.0), Lambda01 = L0,
#'   kappa02 = 2.0, nu02 = 5.0, mu02 = c(0.0, 0.0), Lambda02 = L0,
#'   r = NULL,
#'   ne1 = NULL, ne2 = 15L, alpha01e = NULL, alpha02e = 0.5,
#'   ybar_e1 = NULL, ybar_e2 = c(0.2, 0.1), Se1 = NULL, Se2 = Se2,
#'   nMC = 500L, error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 5L
#' )
#' }
#'
#' # Example 6: Controlled design, posterior probability, NIW prior, MM method
#' \dontrun{
#' Sigma <- matrix(c(4.0, 0.8, 0.8, 1.0), 2, 2)
#' L0    <- matrix(c(8.0, 0.0, 0.0, 2.0), 2, 2)
#' pbayesdecisionprob2cont(
#'   nsim = 100L, prob = 'posterior', design = 'controlled',
#'   prior = 'N-Inv-Wishart',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   gamma1 = 0.8, gamma2 = 0.8,
#'   theta.TV1 = 1.5, theta.MAV1 = 0.5,
#'   theta.TV2 = 1.0, theta.MAV2 = 0.3,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   n1 = 20L, n2 = 20L, m1 = NULL, m2 = NULL,
#'   mu1 = rbind(c(1.0, 0.5), c(2.5, 1.5), c(4.0, 2.5)),
#'   Sigma1 = Sigma,
#'   mu2 = rbind(c(0.0, 0.0), c(0.0, 0.0), c(0.0, 0.0)),
#'   Sigma2 = Sigma,
#'   kappa01 = 2.0, nu01 = 5.0, mu01 = c(2.0, 1.0), Lambda01 = L0,
#'   kappa02 = 2.0, nu02 = 5.0, mu02 = c(0.0, 0.0), Lambda02 = L0,
#'   r = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
#'   ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
#'   nMC = 500L, method = 'MM',
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE, seed = 6L
#' )
#' }
#'
#' @importFrom stats rnorm
#' @export
pbayesdecisionprob2cont <- function(nsim,
                                    prob, design, prior,
                                    GoRegions, NoGoRegions,
                                    gamma1, gamma2,
                                    theta.TV1   = NULL, theta.MAV1  = NULL,
                                    theta.TV2   = NULL, theta.MAV2  = NULL,
                                    theta.NULL1 = NULL, theta.NULL2 = NULL,
                                    n1, n2      = NULL,
                                    m1          = NULL, m2          = NULL,
                                    mu1, Sigma1,
                                    mu2, Sigma2,
                                    kappa01     = NULL, nu01        = NULL,
                                    mu01        = NULL, Lambda01    = NULL,
                                    kappa02     = NULL, nu02        = NULL,
                                    mu02        = NULL, Lambda02    = NULL,
                                    r           = NULL,
                                    ne1         = NULL, ne2         = NULL,
                                    alpha01e    = NULL, alpha02e    = NULL,
                                    ybar_e1     = NULL, ybar_e2     = NULL,
                                    Se1         = NULL, Se2         = NULL,
                                    nMC         = 1000L,
                                    method      = 'MC',
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

  for (nm in c('gamma1', 'gamma2')) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val <= 0 || val >= 1)
      stop(paste0("'", nm, "' must be a single numeric value in (0, 1)"))
  }

  # Threshold validation (delegated to pbayespostpred2cont; check basics here)
  if (prob == 'posterior') {
    for (nm in c('theta.TV1', 'theta.MAV1', 'theta.TV2', 'theta.MAV2')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val))
        stop(paste0("'", nm, "' must be a single numeric when prob = 'posterior'"))
    }
    if (theta.TV1 <= theta.MAV1)
      stop("'theta.TV1' must be strictly greater than 'theta.MAV1'")
    if (theta.TV2 <= theta.MAV2)
      stop("'theta.TV2' must be strictly greater than 'theta.MAV2'")
  } else {
    for (nm in c('theta.NULL1', 'theta.NULL2')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val))
        stop(paste0("'", nm, "' must be a single numeric when prob = 'predictive'"))
    }
    for (nm in c('m1', 'm2')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 1L)
        stop(paste0("'", nm, "' must be a single positive integer when prob = 'predictive'"))
    }
    m1 <- as.integer(m1); m2 <- as.integer(m2)
  }

  if (!is.numeric(n1) || length(n1) != 1L || is.na(n1) ||
      n1 != floor(n1) || n1 < 1L)
    stop("'n1' must be a single positive integer")
  n1 <- as.integer(n1)

  if (design %in% c('controlled', 'external')) {
    if (is.null(n2) || !is.numeric(n2) || length(n2) != 1L || is.na(n2) ||
        n2 != floor(n2) || n2 < 1L)
      stop("'n2' must be a single positive integer for controlled/external design")
    n2 <- as.integer(n2)
  }

  # mu1: coerce to matrix (n_scen x 2)
  mu1 <- if (is.matrix(mu1)) mu1 else matrix(mu1, ncol = 2L)
  if (ncol(mu1) != 2L)
    stop("'mu1' must have 2 columns (one per endpoint)")
  n_scen <- nrow(mu1)

  # mu2: coerce to matrix and check dimensions
  mu2 <- if (is.matrix(mu2)) mu2 else matrix(mu2, ncol = 2L)
  if (ncol(mu2) != 2L)
    stop("'mu2' must have 2 columns")
  if (nrow(mu2) != n_scen)
    stop("'mu1' and 'mu2' must have the same number of rows (scenarios)")

  # Sigma1 and Sigma2 validation
  for (nm in c('Sigma1', 'Sigma2')) {
    val <- get(nm)
    if (!is.matrix(val) || !is.numeric(val) || nrow(val) != 2L ||
        ncol(val) != 2L || any(is.na(val)))
      stop(paste0("'", nm, "' must be a 2x2 numeric matrix"))
  }

  if (!is.logical(error_if_Miss) || length(error_if_Miss) != 1L ||
      is.na(error_if_Miss))
    stop("'error_if_Miss' must be a single logical value")
  if (!is.logical(Gray_inc_Miss) || length(Gray_inc_Miss) != 1L ||
      is.na(Gray_inc_Miss))
    stop("'Gray_inc_Miss' must be a single logical value")

  if (!is.numeric(nMC) || length(nMC) != 1L || is.na(nMC) ||
      nMC != floor(nMC) || nMC < 1L)
    stop("'nMC' must be a single positive integer")
  nMC <- as.integer(nMC)

  if (!is.character(method) || length(method) != 1L ||
      !method %in% c('MC', 'MM'))
    stop("'method' must be either 'MC' or 'MM'")

  # ---------------------------------------------------------------------------
  # Section 2: Pre-generate shared PoC data components
  #
  # For efficiency across scenarios, we separate the random residuals (shared)
  # from the scenario-specific mean shifts (added per scenario), following
  # the same approach as pGNGsinglecontinuous:
  #
  #   y_{k,i,j} = Z_{k,i,j} + mu_k[s]
  #
  # where Z_{k,i,j} ~ N_2(0, Sigma_k) is scenario-independent.
  # Sufficient statistics:
  #   ybar_k[s] = colMeans(Z_k) + mu_k[s]   (mean shifts by mu_k[s])
  #   S_k       = t(Z_k_centered) %*% Z_k_centered  (invariant to mean)
  # ---------------------------------------------------------------------------

  set.seed(seed)

  # Generate nsim x (n1 * 2) standard normal residuals for group 1,
  # then transform by Cholesky factor of Sigma1.
  # Result: Z1_resid[sim, ] = vectorised n1 x 2 residual matrix for replicate sim
  R1 <- chol(Sigma1)
  # raw standard normals: (nsim * n1) x 2
  Z1_raw <- matrix(rnorm(nsim * n1 * 2L), nrow = nsim * n1, ncol = 2L) %*% R1

  # Pre-compute within-replicate column sums and sum-of-squares for group 1.
  # Z1_raw is laid out as nsim blocks of n1 rows each.
  # block_id[row] indicates which simulation replicate the row belongs to.
  block1 <- rep(seq_len(nsim), each = n1)

  # colSums by replicate: (nsim x 2) matrix of raw (mean=0) sums
  Z1_colsums <- apply(Z1_raw, 2L, function(col) tapply(col, block1, sum))

  # Centered residuals for S1: subtract within-replicate column means
  Z1_colmeans_rep <- Z1_colsums[block1, ] / n1
  Z1_centered <- Z1_raw - Z1_colmeans_rep

  # S1 for each replicate: (nsim x 4) where col 1-4 = c(s11, s12, s21, s22)
  # We need S_k as 2x2; store as the 4 elements
  S1_11 <- tapply(Z1_centered[, 1L] ^ 2,        block1, sum)
  S1_12 <- tapply(Z1_centered[, 1L] * Z1_centered[, 2L], block1, sum)
  S1_22 <- tapply(Z1_centered[, 2L] ^ 2,        block1, sum)

  # Group 2 (only for controlled/external)
  if (design %in% c('controlled', 'external')) {
    R2 <- chol(Sigma2)
    Z2_raw <- matrix(rnorm(nsim * n2 * 2L), nrow = nsim * n2, ncol = 2L) %*% R2
    block2 <- rep(seq_len(nsim), each = n2)

    Z2_colsums  <- apply(Z2_raw, 2L, function(col) tapply(col, block2, sum))
    Z2_colmeans_rep <- Z2_colsums[block2, ] / n2
    Z2_centered <- Z2_raw - Z2_colmeans_rep

    S2_11 <- tapply(Z2_centered[, 1L] ^ 2,              block2, sum)
    S2_12 <- tapply(Z2_centered[, 1L] * Z2_centered[, 2L], block2, sum)
    S2_22 <- tapply(Z2_centered[, 2L] ^ 2,              block2, sum)
  }

  # ---------------------------------------------------------------------------
  # Section 3: Simulation loop over scenarios
  # ---------------------------------------------------------------------------

  # Result matrix: columns = Go, NoGo, Miss (n_scen x 3)
  result_mat <- matrix(0, nrow = n_scen, ncol = 3L)

  for (s in seq_len(n_scen)) {

    # Sample means for group 1 in this scenario: add true mean to colSums / n1
    ybar1_sim <- sweep(Z1_colsums / n1, 2L, mu1[s, ], '+')  # nsim x 2

    # S1 is scenario-invariant (depends only on residuals, not means)
    # -- no recomputation needed

    if (design %in% c('controlled', 'external')) {
      ybar2_sim <- sweep(Z2_colsums / n2, 2L, mu2[s, ], '+')  # nsim x 2
    }

    # Vectors of PrGo and PrNoGo across nsim replicates
    PrGo_vec   <- numeric(nsim)
    PrNoGo_vec <- numeric(nsim)

    for (i in seq_len(nsim)) {

      # Assemble sufficient statistics for this replicate
      yb1 <- ybar1_sim[i, ]
      S1i <- matrix(c(S1_11[i], S1_12[i], S1_12[i], S1_22[i]), 2L, 2L)

      if (design %in% c('controlled', 'external')) {
        yb2 <- ybar2_sim[i, ]
        S2i <- matrix(c(S2_11[i], S2_12[i], S2_12[i], S2_22[i]), 2L, 2L)
      } else {
        yb2 <- NULL; S2i <- NULL
      }

      # Call pbayespostpred2cont to get region probability vector
      Pr_R <- pbayespostpred2cont(
        prob = prob, design = design, prior = prior,
        theta.TV1 = theta.TV1, theta.MAV1 = theta.MAV1,
        theta.TV2 = theta.TV2, theta.MAV2 = theta.MAV2,
        theta.NULL1 = theta.NULL1, theta.NULL2 = theta.NULL2,
        n1 = n1, n2 = n2,
        ybar1 = yb1, S1 = S1i,
        ybar2 = yb2, S2 = S2i,
        m1 = m1, m2 = m2,
        kappa01 = kappa01, nu01 = nu01, mu01 = mu01, Lambda01 = Lambda01,
        kappa02 = kappa02, nu02 = nu02, mu02 = mu02, Lambda02 = Lambda02,
        r = r,
        ne1 = ne1, ne2 = ne2, alpha01e = alpha01e, alpha02e = alpha02e,
        ybar_e1 = ybar_e1, ybar_e2 = ybar_e2, Se1 = Se1, Se2 = Se2,
        nMC = nMC, method = method
      )

      PrGo_vec[i]   <- sum(Pr_R[GoRegions])
      PrNoGo_vec[i] <- sum(Pr_R[NoGoRegions])
    }

    # Classify each replicate into Go, NoGo, Miss
    ind_Go   <- (PrGo_vec >= gamma1) & (PrNoGo_vec <  gamma2)
    ind_NoGo <- (PrGo_vec <  gamma1) & (PrNoGo_vec >= gamma2)
    ind_Miss <- (PrGo_vec >= gamma1) & (PrNoGo_vec >= gamma2)

    result_mat[s, 1L] <- mean(ind_Go)
    result_mat[s, 2L] <- mean(ind_NoGo)
    result_mat[s, 3L] <- mean(ind_Miss)
  }

  # ---------------------------------------------------------------------------
  # Section 4: Assemble output
  # ---------------------------------------------------------------------------

  # Suppress floating-point noise near zero
  result_mat[result_mat < .Machine$double.eps ^ 0.25] <- 0

  if (error_if_Miss && any(result_mat[, 3L] > 0))
    stop("Positive Miss probability detected. Please re-consider the chosen thresholds.")

  GrayProb <- if (Gray_inc_Miss) {
    1 - result_mat[, 1L] - result_mat[, 2L]
  } else {
    1 - rowSums(result_mat)
  }

  if (design == 'uncontrolled') {
    results <- data.frame(
      mu1_ep1 = mu1[, 1L],
      mu1_ep2 = mu1[, 2L],
      Go      = result_mat[, 1L],
      Gray    = GrayProb,
      NoGo    = result_mat[, 2L]
    )
  } else {
    results <- data.frame(
      mu1_ep1 = mu1[, 1L],
      mu1_ep2 = mu1[, 2L],
      mu2_ep1 = mu2[, 1L],
      mu2_ep2 = mu2[, 2L],
      Go      = result_mat[, 1L],
      Gray    = GrayProb,
      NoGo    = result_mat[, 2L]
    )
  }

  if (!error_if_Miss && !Gray_inc_Miss)
    results$Miss <- result_mat[, 3L]

  # Attach metadata as attributes
  attr(results, 'prob')          <- prob
  attr(results, 'design')        <- design
  attr(results, 'prior')         <- prior
  attr(results, 'nsim')          <- nsim
  attr(results, 'nMC')           <- nMC
  attr(results, 'method')        <- method
  attr(results, 'GoRegions')     <- GoRegions
  attr(results, 'NoGoRegions')   <- NoGoRegions
  attr(results, 'gamma1')        <- gamma1
  attr(results, 'gamma2')        <- gamma2
  attr(results, 'theta.TV1')     <- theta.TV1
  attr(results, 'theta.MAV1')    <- theta.MAV1
  attr(results, 'theta.TV2')     <- theta.TV2
  attr(results, 'theta.MAV2')    <- theta.MAV2
  attr(results, 'theta.NULL1')   <- theta.NULL1
  attr(results, 'theta.NULL2')   <- theta.NULL2
  attr(results, 'n1')            <- n1
  attr(results, 'n2')            <- n2
  attr(results, 'm1')            <- m1
  attr(results, 'm2')            <- m2
  attr(results, 'Sigma1')        <- Sigma1
  attr(results, 'Sigma2')        <- Sigma2
  attr(results, 'kappa01')       <- kappa01
  attr(results, 'nu01')          <- nu01
  attr(results, 'mu01')          <- mu01
  attr(results, 'Lambda01')      <- Lambda01
  attr(results, 'kappa02')       <- kappa02
  attr(results, 'nu02')          <- nu02
  attr(results, 'mu02')          <- mu02
  attr(results, 'Lambda02')      <- Lambda02
  attr(results, 'r')             <- r
  attr(results, 'ne1')           <- ne1
  attr(results, 'ne2')           <- ne2
  attr(results, 'alpha01e')      <- alpha01e
  attr(results, 'alpha02e')      <- alpha02e
  attr(results, 'ybar_e1')       <- ybar_e1
  attr(results, 'ybar_e2')       <- ybar_e2
  attr(results, 'Se1')           <- Se1
  attr(results, 'Se2')           <- Se2
  attr(results, 'error_if_Miss') <- error_if_Miss
  attr(results, 'Gray_inc_Miss') <- Gray_inc_Miss
  attr(results, 'seed')          <- seed

  class(results) <- c('pbayesdecisionprob2cont', 'data.frame')

  return(results)
}

# ==============================================================================

#' Print Method for pbayesdecisionprob2cont Objects
#'
#' Displays a formatted summary of Go/NoGo/Gray decision probabilities for
#' two-continuous-endpoint results returned by \code{\link{pbayesdecisionprob2cont}}.
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

  fmt <- function(v) {
    if (is.null(v))    return('NULL')
    if (is.matrix(v))  return(paste0('[', paste(v, collapse = ', '), ']'))
    if (length(v) > 1) return(paste0('(', paste(v, collapse = ', '), ')'))
    as.character(v)
  }

  prob          <- attr(x, 'prob')
  design        <- attr(x, 'design')
  prior         <- attr(x, 'prior')
  nsim          <- attr(x, 'nsim')
  nMC           <- attr(x, 'nMC')
  GoRegions     <- attr(x, 'GoRegions')
  NoGoRegions   <- attr(x, 'NoGoRegions')
  gamma1        <- attr(x, 'gamma1')
  gamma2        <- attr(x, 'gamma2')
  n1            <- attr(x, 'n1')
  n2            <- attr(x, 'n2')
  m1            <- attr(x, 'm1')
  m2            <- attr(x, 'm2')
  Sigma1        <- attr(x, 'Sigma1')
  Sigma2        <- attr(x, 'Sigma2')
  kappa01       <- attr(x, 'kappa01')
  nu01          <- attr(x, 'nu01')
  mu01          <- attr(x, 'mu01')
  Lambda01      <- attr(x, 'Lambda01')
  kappa02       <- attr(x, 'kappa02')
  nu02          <- attr(x, 'nu02')
  mu02          <- attr(x, 'mu02')
  Lambda02      <- attr(x, 'Lambda02')
  r             <- attr(x, 'r')
  ne1           <- attr(x, 'ne1')
  ne2           <- attr(x, 'ne2')
  alpha01e      <- attr(x, 'alpha01e')
  alpha02e      <- attr(x, 'alpha02e')
  ybar_e1       <- attr(x, 'ybar_e1')
  ybar_e2       <- attr(x, 'ybar_e2')
  error_if_Miss <- attr(x, 'error_if_Miss')
  Gray_inc_Miss <- attr(x, 'Gray_inc_Miss')
  seed          <- attr(x, 'seed')
  method        <- attr(x, 'method')

  if (prob == 'posterior') {
    thresh_str <- sprintf('TV1 = %s, MAV1 = %s, TV2 = %s, MAV2 = %s',
                          fmt(attr(x, 'theta.TV1')), fmt(attr(x, 'theta.MAV1')),
                          fmt(attr(x, 'theta.TV2')), fmt(attr(x, 'theta.MAV2')))
  } else {
    thresh_str <- sprintf('NULL1 = %s, NULL2 = %s',
                          fmt(attr(x, 'theta.NULL1')), fmt(attr(x, 'theta.NULL2')))
  }

  prior_label <- if (prior == 'vague') 'vague' else 'N-Inv-Wishart'

  cat('Go/NoGo/Gray Decision Probabilities (Two Continuous Endpoints)\n')
  cat(strrep('-', 65), '\n')
  cat(sprintf('  Probability type : %s\n',   prob))
  cat(sprintf('  Design           : %s\n',   design))
  cat(sprintf('  Prior            : %s\n',   prior_label))
  cat(sprintf('  Simulations      : nsim = %s, nMC = %s\n', fmt(nsim), fmt(nMC)))
  cat(sprintf('  Method           : %s\n',   fmt(method)))
  cat(sprintf('  Seed             : %s\n',   fmt(seed)))
  cat(sprintf('  Threshold(s)     : %s\n',   thresh_str))
  cat(sprintf('  Go  threshold    : gamma1 = %s\n', fmt(gamma1)))
  cat(sprintf('  NoGo threshold   : gamma2 = %s\n', fmt(gamma2)))
  cat(sprintf('  Go  regions      : {%s}\n', paste(GoRegions,   collapse = ', ')))
  cat(sprintf('  NoGo regions     : {%s}\n', paste(NoGoRegions, collapse = ', ')))
  cat(sprintf('  Sample size      : n1 = %s, n2 = %s\n', fmt(n1), fmt(n2)))

  if (prior == 'N-Inv-Wishart') {
    cat(sprintf('  Prior Grp1 (NIW) : kappa01 = %s, nu01 = %s\n',
                fmt(kappa01), fmt(nu01)))
    cat(sprintf('                     mu01 = %s, Lambda01 = %s\n',
                fmt(mu01), fmt(Lambda01)))
    if (design %in% c('controlled', 'external')) {
      cat(sprintf('  Prior Grp2 (NIW) : kappa02 = %s, nu02 = %s\n',
                  fmt(kappa02), fmt(nu02)))
      cat(sprintf('                     mu02 = %s, Lambda02 = %s\n',
                  fmt(mu02), fmt(Lambda02)))
    }
  }

  if (design == 'uncontrolled')
    cat(sprintf('  Hyp. control     : mu02 = %s, r = %s\n', fmt(mu02), fmt(r)))

  if (prob == 'predictive')
    cat(sprintf('  Future trial     : m1 = %s, m2 = %s\n', fmt(m1), fmt(m2)))

  if (design == 'external') {
    cat(sprintf('  External data    : ne1 = %s, ne2 = %s\n', fmt(ne1), fmt(ne2)))
    cat(sprintf('                     alpha01e = %s, alpha02e = %s\n',
                fmt(alpha01e), fmt(alpha02e)))
    if (!is.null(ybar_e1))
      cat(sprintf('                     ybar_e1 = %s\n', fmt(ybar_e1)))
    if (!is.null(ybar_e2))
      cat(sprintf('                     ybar_e2 = %s\n', fmt(ybar_e2)))
  }

  cat(sprintf('  Miss handling    : error_if_Miss = %s, Gray_inc_Miss = %s\n',
              fmt(error_if_Miss), fmt(Gray_inc_Miss)))
  cat(strrep('-', 65), '\n')

  scenario_cols <- c('mu1_ep1', 'mu1_ep2', 'mu2_ep1', 'mu2_ep2')
  prob_cols     <- names(x)[!names(x) %in% scenario_cols]

  x_print <- x
  x_print[prob_cols] <- lapply(x[prob_cols], function(col) {
    formatC(col, digits = digits, format = 'f')
  })

  print.data.frame(x_print, row.names = FALSE, quote = FALSE)
  cat(strrep('-', 65), '\n')

  invisible(x)
}
