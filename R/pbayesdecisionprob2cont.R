#' Go/NoGo/Gray Decision Probabilities for Two Continuous Endpoints
#'
#' Estimates Go, NoGo, and Gray operating characteristics for a two-endpoint
#' continuous Bayesian Go/NoGo decision framework by Monte Carlo simulation.
#' For each simulated dataset, region probabilities are computed by calling
#' \code{\link{pbayespostpred2cont}} in vectorised mode, consistent with the
#' single-endpoint analogue \code{\link{pbayesdecisionprob1cont}}.
#'
#' @param nsim A positive integer. Number of simulated datasets per scenario.
#' @param prob Character scalar, either \code{'posterior'} or
#'        \code{'predictive'}.
#' @param design Character scalar: \code{'controlled'}, \code{'uncontrolled'},
#'        or \code{'external'}.
#' @param prior Character scalar: \code{'vague'} or \code{'N-Inv-Wishart'}.
#' @param GoRegions Integer vector. Region indices (1 to 9 for posterior,
#'        1 to 4 for predictive) that define a Go decision.
#' @param NoGoRegions Integer vector. Region indices that define a NoGo
#'        decision. Must be disjoint from \code{GoRegions}.
#' @param gamma1 Numeric scalar in (0, 1). Go threshold: a Go decision is made
#'        when \eqn{PrGo \ge \gamma_1} and \eqn{PrNoGo < \gamma_2}.
#' @param gamma2 Numeric scalar in (0, 1). NoGo threshold: a NoGo decision is
#'        made when \eqn{PrNoGo \ge \gamma_2} and \eqn{PrGo < \gamma_1}.
#' @param theta.TV1 Numeric scalar. Target value for Endpoint 1
#'        (required when \code{prob = 'posterior'}).
#' @param theta.MAV1 Numeric scalar. Minimum acceptable value for Endpoint 1
#'        (required when \code{prob = 'posterior'}).
#' @param theta.TV2 Numeric scalar. Target value for Endpoint 2
#'        (required when \code{prob = 'posterior'}).
#' @param theta.MAV2 Numeric scalar. Minimum acceptable value for Endpoint 2
#'        (required when \code{prob = 'posterior'}).
#' @param theta.NULL1 Numeric scalar. Null threshold for Endpoint 1
#'        (required when \code{prob = 'predictive'}).
#' @param theta.NULL2 Numeric scalar. Null threshold for Endpoint 2
#'        (required when \code{prob = 'predictive'}).
#' @param n1 Positive integer. PoC sample size for the treatment arm.
#' @param n2 Positive integer or \code{NULL}. PoC sample size for the control
#'        arm (not used when \code{design = 'uncontrolled'}).
#' @param m1 Positive integer or \code{NULL}. Pivotal study sample size for
#'        the treatment arm (required when \code{prob = 'predictive'}).
#' @param m2 Positive integer or \code{NULL}. Pivotal study sample size for
#'        the control arm (required when \code{prob = 'predictive'} and
#'        \code{design != 'uncontrolled'}).
#' @param mu1 Numeric matrix with \eqn{K} rows and 2 columns. Each row
#'        specifies the true treatment mean vector for one scenario.
#' @param Sigma1 A 2x2 positive-definite numeric matrix. True covariance
#'        matrix for the treatment arm (shared across scenarios).
#' @param mu2 Numeric matrix with \eqn{K} rows and 2 columns or a length-2
#'        vector. True control mean vector(s). For
#'        \code{design = 'uncontrolled'}, only a single row or vector is
#'        used (hypothetical control mean, not used in simulation).
#' @param Sigma2 A 2x2 positive-definite numeric matrix. True covariance
#'        matrix for the control arm. Set to a diagonal matrix or any valid
#'        matrix for \code{design = 'uncontrolled'} (not used in simulation).
#' @param kappa01 Positive numeric scalar. NIW prior concentration for the
#'        treatment arm (required when \code{prior = 'N-Inv-Wishart'}).
#' @param nu01 Numeric scalar \eqn{> 3}. NIW prior degrees of freedom for
#'        the treatment arm (required when \code{prior = 'N-Inv-Wishart'}).
#' @param mu01 Length-2 numeric vector. NIW prior mean for the treatment arm
#'        (required when \code{prior = 'N-Inv-Wishart'}).
#' @param Lambda01 A 2x2 positive-definite numeric matrix. NIW prior scale
#'        matrix for the treatment arm (required when
#'        \code{prior = 'N-Inv-Wishart'}).
#' @param kappa02 Positive numeric scalar. NIW prior concentration for the
#'        control arm (required for controlled/external when
#'        \code{prior = 'N-Inv-Wishart'}).
#' @param nu02 Numeric scalar \eqn{> 3}. NIW prior degrees of freedom for
#'        the control arm.
#' @param mu02 Length-2 numeric vector. NIW prior mean for the control arm
#'        or hypothetical control location
#'        (required when \code{design = 'uncontrolled'}).
#' @param Lambda02 A 2x2 positive-definite numeric matrix. NIW prior scale
#'        matrix for the control arm.
#' @param r Positive numeric scalar. Variance scaling factor for the
#'        hypothetical control (required when
#'        \code{design = 'uncontrolled'}).
#' @param ne1 Positive integer or \code{NULL}. External treatment sample size
#'        (used when \code{design = 'external'}).
#' @param ne2 Positive integer or \code{NULL}. External control sample size
#'        (used when \code{design = 'external'}).
#' @param alpha01e Numeric scalar in \code{(0, 1]}. Power prior weight for
#'        external treatment data.
#' @param alpha02e Numeric scalar in \code{(0, 1]}. Power prior weight for
#'        external control data.
#' @param ybar_e1 Length-2 numeric vector. External treatment sample mean.
#' @param ybar_e2 Length-2 numeric vector. External control sample mean.
#' @param Se1 A 2x2 numeric matrix. External treatment sum-of-squares matrix.
#' @param Se2 A 2x2 numeric matrix. External control sum-of-squares matrix.
#' @param nMC Positive integer. Number of Monte Carlo draws used by
#'        \code{\link{pbayespostpred2cont}} to compute region probabilities
#'        for each simulated dataset. Default \code{1000L}.
#' @param method Character scalar: \code{'MC'} (default) or \code{'MM'}
#'        (Moment-Matching via \code{mvtnorm::pmvt}).
#' @param error_if_Miss Logical scalar. If \code{TRUE} (default), the function
#'        stops when any scenario yields a positive Miss probability.
#' @param Gray_inc_Miss Logical scalar. If \code{TRUE}, Miss probability is
#'        folded into Gray. If \code{FALSE} (default), Miss is reported
#'        separately. Active only when \code{error_if_Miss = FALSE}.
#' @param seed A numeric value for reproducible random number generation.
#'
#' @return A \code{data.frame} with one row per scenario and columns
#'   \code{mu1_ep1}, \code{mu1_ep2} (and \code{mu2_ep1}, \code{mu2_ep2} for
#'   controlled/external), \code{Go}, \code{Gray}, \code{NoGo}, and optionally
#'   \code{Miss}.  The object has S3 class
#'   \code{c("pbayesdecisionprob2cont", "data.frame")} with all input
#'   parameters attached as attributes.
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
                                    prob,
                                    design,
                                    prior,
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
      n1 != floor(n1) || n1 < 2L)
    stop("'n1' must be a single integer >= 2")
  n1 <- as.integer(n1)

  if (design %in% c('controlled', 'external')) {
    if (is.null(n2) || !is.numeric(n2) || length(n2) != 1L || is.na(n2) ||
        n2 != floor(n2) || n2 < 2L)
      stop("'n2' must be a single integer >= 2 for controlled/external designs")
    n2 <- as.integer(n2)
  }

  # mu1 / mu2: accept matrix or vector; coerce to matrix
  if (is.numeric(mu1) && !is.matrix(mu1)) mu1 <- matrix(mu1, nrow = 1L)
  if (!is.matrix(mu1) || ncol(mu1) != 2L || !is.numeric(mu1))
    stop("'mu1' must be a numeric matrix with 2 columns")
  n_scen <- nrow(mu1)

  if (is.numeric(mu2) && !is.matrix(mu2)) mu2 <- matrix(mu2, nrow = 1L)
  if (!is.matrix(mu2) || ncol(mu2) != 2L || !is.numeric(mu2))
    stop("'mu2' must be a numeric matrix with 2 columns")
  if (design %in% c('controlled', 'external') && nrow(mu2) != n_scen)
    stop("'mu2' must have the same number of rows as 'mu1'")

  for (nm in c('Sigma1', 'Sigma2')) {
    val <- get(nm)
    if (!is.matrix(val) || !is.numeric(val) || nrow(val) != 2L || ncol(val) != 2L)
      stop(paste0("'", nm, "' must be a 2x2 numeric matrix"))
  }

  if (prior == 'N-Inv-Wishart') {
    for (nm in c('kappa01', 'nu01')) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L || is.na(val) ||
          val <= 0)
        stop(paste0("'", nm, "' must be a single positive numeric for NIW prior"))
    }
    if (is.null(mu01) || !is.numeric(mu01) || length(mu01) != 2L)
      stop("'mu01' must be a length-2 numeric vector for NIW prior")
    if (is.null(Lambda01) || !is.matrix(Lambda01) || nrow(Lambda01) != 2L)
      stop("'Lambda01' must be a 2x2 numeric matrix for NIW prior")
    if (design %in% c('controlled', 'external')) {
      for (nm in c('kappa02', 'nu02')) {
        val <- get(nm)
        if (is.null(val) || !is.numeric(val) || length(val) != 1L ||
            is.na(val) || val <= 0)
          stop(paste0("'", nm, "' must be a single positive numeric for NIW prior"))
      }
      if (is.null(mu02) || !is.numeric(mu02) || length(mu02) != 2L)
        stop("'mu02' must be a length-2 numeric vector for NIW prior")
      if (is.null(Lambda02) || !is.matrix(Lambda02) || nrow(Lambda02) != 2L)
        stop("'Lambda02' must be a 2x2 numeric matrix for NIW prior")
    }
  }

  if (design == 'uncontrolled') {
    if (is.null(r) || !is.numeric(r) || length(r) != 1L || is.na(r) || r <= 0)
      stop("'r' must be a single positive numeric when design = 'uncontrolled'")
    if (is.null(mu02) || !is.numeric(mu02) || length(mu02) != 2L)
      stop("'mu02' must be a length-2 numeric vector when design = 'uncontrolled'")
  }

  if (design == 'external') {
    has_ext1 <- !is.null(ne1) && !is.null(alpha01e) &&
      !is.null(ybar_e1) && !is.null(Se1)
    has_ext2 <- !is.null(ne2) && !is.null(alpha02e) &&
      !is.null(ybar_e2) && !is.null(Se2)
    if (!has_ext1 && !has_ext2)
      stop(paste0("For design = 'external', at least one complete set of ",
                  "external data must be provided"))
  }

  if (!is.numeric(nMC) || length(nMC) != 1L || is.na(nMC) ||
      nMC != floor(nMC) || nMC < 1L)
    stop("'nMC' must be a single positive integer")
  nMC <- as.integer(nMC)

  if (!is.character(method) || length(method) != 1L ||
      !method %in% c('MC', 'MM'))
    stop("'method' must be either 'MC' or 'MM'")

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

  R_Sigma1 <- chol(Sigma1)
  Z1_raw   <- matrix(rnorm(nsim * n1 * 2L),
                     nrow = nsim * n1, ncol = 2L) %*% R_Sigma1
  block1   <- rep(seq_len(nsim), each = n1)

  Z1_colsums      <- apply(Z1_raw, 2L, function(col) tapply(col, block1, sum))
  Z1_colmeans_rep <- Z1_colsums[block1, ] / n1
  Z1_centered     <- Z1_raw - Z1_colmeans_rep

  # S1 stored as 3 unique elements per replicate (symmetric 2x2)
  S1_11 <- tapply(Z1_centered[, 1L] ^ 2,                block1, sum)
  S1_12 <- tapply(Z1_centered[, 1L] * Z1_centered[, 2L], block1, sum)
  S1_22 <- tapply(Z1_centered[, 2L] ^ 2,                block1, sum)

  if (design %in% c('controlled', 'external')) {
    R_Sigma2 <- chol(Sigma2)
    Z2_raw   <- matrix(rnorm(nsim * n2 * 2L),
                       nrow = nsim * n2, ncol = 2L) %*% R_Sigma2
    block2   <- rep(seq_len(nsim), each = n2)

    Z2_colsums      <- apply(Z2_raw, 2L, function(col) tapply(col, block2, sum))
    Z2_colmeans_rep <- Z2_colsums[block2, ] / n2
    Z2_centered     <- Z2_raw - Z2_colmeans_rep

    S2_11 <- tapply(Z2_centered[, 1L] ^ 2,                block2, sum)
    S2_12 <- tapply(Z2_centered[, 1L] * Z2_centered[, 2L], block2, sum)
    S2_22 <- tapply(Z2_centered[, 2L] ^ 2,                block2, sum)
  }

  # ---------------------------------------------------------------------------
  # Section 3: Scenario loop
  #
  # For each scenario s:
  #   1. Construct nsim sufficient statistics by shifting the pre-generated
  #      residuals by the scenario mean.
  #   2. Call pbayespostpred2cont once in vectorised mode (ybar1 as matrix,
  #      S1 as list) to obtain an nsim x n_regions matrix of Pr_R values.
  #   3. Compute PrGo and PrNoGo for each replicate and classify.
  # ---------------------------------------------------------------------------

  result_mat <- matrix(0, nrow = n_scen, ncol = 3L)

  for (s in seq_len(n_scen)) {

    # Shift residual column sums by scenario mean (vectorised over nsim)
    ybar1_sim <- sweep(Z1_colsums / n1, 2L, mu1[s, ], '+')
    # ybar1_sim: nsim x 2 matrix

    # Build S1 list for this scenario (scenario-invariant values)
    S1_list <- vector('list', nsim)
    for (i in seq_len(nsim)) {
      S1_list[[i]] <- matrix(c(S1_11[i], S1_12[i], S1_12[i], S1_22[i]),
                             nrow = 2L, ncol = 2L)
    }

    if (design %in% c('controlled', 'external')) {
      ybar2_sim <- sweep(Z2_colsums / n2, 2L, mu2[s, ], '+')
      S2_list <- vector('list', nsim)
      for (i in seq_len(nsim)) {
        S2_list[[i]] <- matrix(c(S2_11[i], S2_12[i], S2_12[i], S2_22[i]),
                               nrow = 2L, ncol = 2L)
      }
    } else {
      ybar2_sim <- NULL
      S2_list   <- NULL
    }

    # Vectorised call: returns nsim x n_regions matrix
    Pr_R_mat <- pbayespostpred2cont(
      prob        = prob,
      design      = design,
      prior       = prior,
      theta.TV1   = theta.TV1,   theta.MAV1  = theta.MAV1,
      theta.TV2   = theta.TV2,   theta.MAV2  = theta.MAV2,
      theta.NULL1 = theta.NULL1, theta.NULL2 = theta.NULL2,
      n1 = n1, n2 = n2,
      ybar1 = ybar1_sim, S1 = S1_list,
      ybar2 = ybar2_sim, S2 = S2_list,
      m1 = m1, m2 = m2,
      kappa01  = kappa01,  nu01     = nu01,
      mu01     = mu01,     Lambda01 = Lambda01,
      kappa02  = kappa02,  nu02     = nu02,
      mu02     = mu02,     Lambda02 = Lambda02,
      r        = r,
      ne1 = ne1, ne2 = ne2,
      alpha01e = alpha01e, alpha02e = alpha02e,
      ybar_e1  = ybar_e1,  ybar_e2  = ybar_e2,
      Se1      = Se1,      Se2      = Se2,
      nMC    = nMC,
      method = method
    )
    # Pr_R_mat: nsim x n_regions (column names R1...R9 or R1...R4)

    # Go/NoGo probability per replicate
    PrGo_vec   <- rowSums(Pr_R_mat[, GoRegions,   drop = FALSE])
    PrNoGo_vec <- rowSums(Pr_R_mat[, NoGoRegions, drop = FALSE])

    # Classify into Go, NoGo, Miss
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

  # Suppress floating-point noise in output probability columns
  prob_cols <- c('Go', 'Gray', 'NoGo', 'Miss')
  prob_cols <- prob_cols[prob_cols %in% names(results)]
  results[prob_cols] <- lapply(results[prob_cols], function(col) {
    ifelse(col < .Machine$double.eps ^ 0.25, 0, col)
  })

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
