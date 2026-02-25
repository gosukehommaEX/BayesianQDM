#' Find Optimal Go/NoGo Thresholds for Two Continuous Endpoints
#'
#' Computes the optimal Go threshold \eqn{\gamma_1} and NoGo threshold
#' \eqn{\gamma_2} for two continuous endpoints by searching over a
#' two-dimensional grid of candidate value pairs.  The search follows a
#' two-stage simulate-then-sweep strategy: \code{nsim} bivariate datasets
#' are generated from the specified scenario parameters, region probabilities
#' \eqn{\hat{g}_{Go,i}} and \eqn{\hat{g}_{NoGo,i}} are computed once for
#' each replicate via \code{\link{pbayespostpred2cont}}, and operating
#' characteristics under each \eqn{(\gamma_1, \gamma_2)} pair are obtained
#' by a fast proportion-based evaluation without additional simulation.
#'
#' @param nsim A positive integer giving the number of Monte Carlo
#'        datasets to simulate.  Default is \code{10000L}.
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or
#'        \code{'external'}.
#' @param prior A character string specifying the prior distribution.
#'        Must be \code{'vague'} or \code{'N-Inv-Wishart'}.
#' @param GoRegions An integer vector of region indices (subset of
#'        \code{1:9}) that constitute the Go region.  The 9 regions are
#'        defined by the cross-classification of the two treatment effects
#'        \eqn{(\theta_1, \theta_2)} relative to \code{(TV1, MAV1)} and
#'        \code{(TV2, MAV2)}: region \eqn{k = (z_1 - 1) \times 3 + z_2}
#'        where \eqn{z_i = 1} if \eqn{\theta_i > TV_i}, \eqn{z_i = 2} if
#'        \eqn{MAV_i < \theta_i \le TV_i}, \eqn{z_i = 3} if
#'        \eqn{\theta_i \le MAV_i}.
#' @param NoGoRegions An integer vector of region indices (subset of
#'        \code{1:9}) that constitute the NoGo region.
#'        Must be disjoint from \code{GoRegions}.
#' @param mu1 A length-2 numeric vector giving the true bivariate mean
#'        for group 1 (treatment arm) for the scenario used to search
#'        for the optimal thresholds.
#' @param Sigma1 A 2x2 positive-definite numeric matrix giving the
#'        true within-arm covariance in the treatment arm.
#' @param mu2 A length-2 numeric vector giving the true bivariate mean
#'        for group 2 (control arm).  Required for
#'        \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param Sigma2 A 2x2 positive-definite numeric matrix giving the
#'        true within-arm covariance in the control arm.  Required for
#'        \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param target_go A numeric scalar in \code{(0, 1)} giving the target
#'        value for \eqn{\Pr(\mathrm{Go})} used to select the optimal
#'        \eqn{\gamma_1}.  For each candidate \eqn{\gamma_1}, the
#'        worst-case \eqn{\Pr(\mathrm{Go})} over all \eqn{\gamma_2} in
#'        \code{gamma2_grid} is compared against \code{target_go} using
#'        the operator specified by \code{crit_go}.
#' @param target_nogo A numeric scalar in \code{(0, 1)} giving the
#'        target value for \eqn{\Pr(\mathrm{NoGo})} used to select the
#'        optimal \eqn{\gamma_2}.  For each candidate \eqn{\gamma_2},
#'        the worst-case \eqn{\Pr(\mathrm{NoGo})} over all
#'        \eqn{\gamma_1} in \code{gamma1_grid} is compared against
#'        \code{target_nogo} using the operator specified by
#'        \code{crit_nogo}.
#' @param crit_go A character string specifying the comparison operator
#'        applied to \eqn{\Pr(\mathrm{Go})} when searching for
#'        \eqn{\gamma_1}.  Must be one of \code{"<"}, \code{"<="},
#'        \code{">"}, or \code{">="}.  Default is \code{"<"}.
#' @param crit_nogo A character string specifying the comparison
#'        operator applied to \eqn{\Pr(\mathrm{NoGo})} when searching
#'        for \eqn{\gamma_2}.  Must be one of \code{"<"}, \code{"<="},
#'        \code{">"}, or \code{">="}.  Default is \code{"<"}.
#' @param sel_go A character string specifying whether to select the
#'        \code{"smallest"} or \code{"largest"} value in
#'        \code{gamma1_grid} satisfying the \code{crit_go} criterion.
#'        Default is \code{"smallest"}.
#' @param sel_nogo A character string specifying whether to select the
#'        \code{"smallest"} or \code{"largest"} value in
#'        \code{gamma2_grid} satisfying the \code{crit_nogo} criterion.
#'        Default is \code{"largest"}.
#' @param n1 A positive integer giving the number of patients in group 1
#'        (treatment arm) in the PoC trial.
#' @param n2 A positive integer giving the number of patients in group 2
#'        (control arm) in the PoC trial.  Set to \code{NULL} for
#'        \code{design = 'uncontrolled'}.
#' @param theta.TV1 A numeric scalar giving the TV threshold for
#'        Endpoint 1.  Required when \code{prob = 'posterior'};
#'        otherwise set to \code{NULL}.
#' @param theta.MAV1 A numeric scalar giving the MAV threshold for
#'        Endpoint 1.  Required when \code{prob = 'posterior'};
#'        otherwise set to \code{NULL}.
#' @param theta.TV2 A numeric scalar giving the TV threshold for
#'        Endpoint 2.  Required when \code{prob = 'posterior'};
#'        otherwise set to \code{NULL}.
#' @param theta.MAV2 A numeric scalar giving the MAV threshold for
#'        Endpoint 2.  Required when \code{prob = 'posterior'};
#'        otherwise set to \code{NULL}.
#' @param theta.NULL1 A numeric scalar giving the null hypothesis
#'        threshold for Endpoint 1.  Required when
#'        \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param theta.NULL2 A numeric scalar giving the null hypothesis
#'        threshold for Endpoint 2.  Required when
#'        \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param m1 A positive integer giving the future sample size for
#'        group 1.  Required when \code{prob = 'predictive'};
#'        set to \code{NULL} otherwise.
#' @param m2 A positive integer giving the future sample size for
#'        group 2.  Required when \code{prob = 'predictive'};
#'        set to \code{NULL} otherwise.
#' @param kappa01 Positive numeric scalar.  NIW prior hyperparameter
#'        \eqn{\kappa_{01}} for group 1.  Required when
#'        \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param nu01 Positive numeric scalar.  NIW prior degrees of freedom
#'        \eqn{\nu_{01}} for group 1.  Required when
#'        \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param mu01 Length-2 numeric vector.  NIW prior mean \eqn{\mu_{01}}
#'        for group 1.  Required when \code{prior = 'N-Inv-Wishart'};
#'        otherwise \code{NULL}.
#' @param Lambda01 A 2x2 positive-definite numeric matrix.  NIW prior
#'        scale matrix \eqn{\Lambda_{01}} for group 1.  Required when
#'        \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param kappa02 Positive numeric scalar; see \code{kappa01}.
#'        For group 2.
#' @param nu02 Positive numeric scalar; see \code{nu01}.
#'        For group 2.
#' @param mu02 Length-2 numeric vector; see \code{mu01}.  For group 2.
#'        May be required for the vague prior uncontrolled design;
#'        see \code{\link{pbayesdecisionprob2cont}}.
#' @param Lambda02 A 2x2 matrix; see \code{Lambda01}.
#'        For group 2.
#' @param r A positive numeric scalar giving the power prior weight for
#'        group 2 when \code{design = 'uncontrolled'} and
#'        \code{prior = 'vague'}.  Otherwise \code{NULL}.
#' @param ne1 A positive integer giving the external treatment sample
#'        size.  Required when \code{design = 'external'} and external
#'        treatment data are used; otherwise \code{NULL}.
#' @param ne2 A positive integer giving the external control sample
#'        size.  Required when \code{design = 'external'} and external
#'        control data are used; otherwise \code{NULL}.
#' @param alpha01e A numeric scalar in \code{(0, 1]} giving the power
#'        prior weight for external treatment data.  Required when
#'        external treatment data are used; otherwise \code{NULL}.
#' @param alpha02e A numeric scalar in \code{(0, 1]} giving the power
#'        prior weight for external control data.  Required when
#'        external control data are used; otherwise \code{NULL}.
#' @param ybar_e1 A length-2 numeric vector.  External treatment sample
#'        mean.  Required when external treatment data are used;
#'        otherwise \code{NULL}.
#' @param ybar_e2 A length-2 numeric vector.  External control sample
#'        mean.  Required when external control data are used;
#'        otherwise \code{NULL}.
#' @param Se1 A 2x2 numeric matrix.  External treatment
#'        sum-of-squares matrix.  Required when external treatment data
#'        are used; otherwise \code{NULL}.
#' @param Se2 A 2x2 numeric matrix.  External control sum-of-squares
#'        matrix.  Required when external control data are used;
#'        otherwise \code{NULL}.
#' @param nMC A positive integer giving the number of Monte Carlo draws
#'        passed to \code{\link{pbayespostpred2cont}} when
#'        \code{method = 'MC'}.  Set to \code{NULL} when
#'        \code{method = 'MM'}.  Default is \code{NULL}.
#' @param method A character string specifying the computation method
#'        passed to \code{\link{pbayespostpred2cont}}.  Must be
#'        \code{'MC'} (default) or \code{'MM'}.
#' @param gamma1_grid A numeric vector of candidate Go threshold values
#'        in \code{(0, 1)} to search over.  Defaults to
#'        \code{seq(0.01, 0.99, by = 0.01)}.
#' @param gamma2_grid A numeric vector of candidate NoGo threshold
#'        values in \code{(0, 1)} to search over.  Defaults to
#'        \code{seq(0.01, 0.99, by = 0.01)}.
#' @param seed An integer value for reproducible random number
#'        generation.  Default is \code{NULL} (no seed set).
#'
#' @return A list of class \code{getgamma2cont} with the following
#'   elements:
#' \describe{
#'   \item{gamma1}{Optimal Go threshold selected from
#'         \code{gamma1_grid} according to \code{crit_go} and
#'         \code{sel_go}.  \code{NA} if no value satisfies the
#'         criterion.}
#'   \item{gamma2}{Optimal NoGo threshold selected from
#'         \code{gamma2_grid} according to \code{crit_nogo} and
#'         \code{sel_nogo}.  \code{NA} if no value satisfies the
#'         criterion.}
#'   \item{PrGo_at_gamma}{Pr(Go) evaluated at
#'         \code{(gamma1, gamma2)}.  \code{NA} if either threshold
#'         is \code{NA}.}
#'   \item{PrNoGo_at_gamma}{Pr(NoGo) evaluated at
#'         \code{(gamma1, gamma2)}.  \code{NA} if either threshold
#'         is \code{NA}.}
#'   \item{gamma1_grid}{The candidate Go threshold grid used.}
#'   \item{gamma2_grid}{The candidate NoGo threshold grid used.}
#'   \item{PrGo_grid}{Numeric matrix of dimensions
#'         \code{length(gamma1_grid)} x \code{length(gamma2_grid)}
#'         giving Pr(Go) at each grid point.}
#'   \item{PrNoGo_grid}{Numeric matrix of the same dimensions giving
#'         Pr(NoGo) at each grid point.}
#' }
#'
#' @details
#' The function uses a two-stage simulate-then-sweep strategy:
#'
#' \strong{Stage 1 (simulation and precomputation)}: \code{nsim}
#' bivariate datasets are generated by drawing from
#' \eqn{N_2(\mu_1, \Sigma_1)} (and \eqn{N_2(\mu_2, \Sigma_2)} for
#' controlled/external designs).  Sufficient statistics
#' \eqn{(\bar{y}_{1,i}, S_{1,i})} are computed for each replicate, and
#' \code{\link{pbayespostpred2cont}} is called once in vectorised mode
#' to return an \eqn{nsim \times 9} matrix of region probabilities.
#' The probabilities are summed over \code{GoRegions} and
#' \code{NoGoRegions} to obtain \eqn{\hat{g}_{Go,i}} and
#' \eqn{\hat{g}_{NoGo,i}} for each replicate, independent of the
#' decision thresholds \eqn{(\gamma_1, \gamma_2)}.
#'
#' \strong{Stage 2 (gamma sweep)}: For each pair
#' \eqn{(\gamma_1, \gamma_2)} in the two-dimensional grid, the
#' operating characteristics are computed as:
#' \deqn{\Pr(\mathrm{Go}) = \frac{1}{n_{sim}} \sum_{i=1}^{n_{sim}}
#'   \mathbf{1}\!\left[\hat{g}_{Go,i} \ge \gamma_1,\;
#'   \hat{g}_{NoGo,i} < \gamma_2\right]}
#' \deqn{\Pr(\mathrm{NoGo}) = \frac{1}{n_{sim}} \sum_{i=1}^{n_{sim}}
#'   \mathbf{1}\!\left[\hat{g}_{NoGo,i} \ge \gamma_2,\;
#'   \hat{g}_{Go,i} < \gamma_1\right]}
#' No additional Monte Carlo draws are required for this step.
#'
#' \strong{Stage 3 (optimal threshold selection)}: For each candidate
#' \eqn{\gamma_1}, the maximum \eqn{\Pr(\mathrm{Go})} over all
#' \eqn{\gamma_2} in \code{gamma2_grid} is computed; the optimal
#' \eqn{\gamma_1} is then chosen as the \code{sel_go} value satisfying
#' the \code{crit_go} comparison against \code{target_go}.  Analogously,
#' the optimal \eqn{\gamma_2} is the \code{sel_nogo} value for which
#' the maximum \eqn{\Pr(\mathrm{NoGo})} over all \eqn{\gamma_1}
#' satisfies \code{crit_nogo} against \code{target_nogo}.
#'
#' @examples
#' # Example 1: Controlled design, posterior probability, vague prior
#' # Scenario: mu1 - mu2 = c(5, 1) at MAV boundary (intermediate effect)
#' # gamma1: smallest gamma1 s.t. max_{gamma2} Pr(Go)   < 0.05
#' # gamma2: largest  gamma2 s.t. max_{gamma1} Pr(NoGo) < 0.20
#' \dontrun{
#' Sigma <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' getgamma2cont(
#'   nsim = 1000L, prob = 'posterior', design = 'controlled',
#'   prior = 'vague',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   mu1 = c(-5.0, 0.0), Sigma1 = Sigma,
#'   mu2 = c(-10.0, -1.0), Sigma2 = Sigma,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n1 = 30L, n2 = 30L,
#'   theta.TV1 = 10.0, theta.MAV1 = 5.0,
#'   theta.TV2 = 2.0,  theta.MAV2 = 1.0,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   m1 = NULL, m2 = NULL,
#'   kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
#'   kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
#'   r = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
#'   ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
#'   nMC = 500L, method = 'MC',
#'   gamma1_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma2_grid = seq(0.01, 0.99, by = 0.01),
#'   seed = 1L
#' )
#' }
#'
#' # Example 2: Controlled design, posterior probability, N-Inv-Wishart prior
#' \dontrun{
#' Sigma  <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' Lambda <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' getgamma2cont(
#'   nsim = 1000L, prob = 'posterior', design = 'controlled',
#'   prior = 'N-Inv-Wishart',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   mu1 = c(-5.0, 0.0), Sigma1 = Sigma,
#'   mu2 = c(-10.0, -1.0), Sigma2 = Sigma,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n1 = 30L, n2 = 30L,
#'   theta.TV1 = 10.0, theta.MAV1 = 5.0,
#'   theta.TV2 = 2.0,  theta.MAV2 = 1.0,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   m1 = NULL, m2 = NULL,
#'   kappa01 = 0.1, nu01 = 4.0, mu01 = c(0.0, 1.0),  Lambda01 = Lambda,
#'   kappa02 = 0.1, nu02 = 4.0, mu02 = c(-10.0, -1.0), Lambda02 = Lambda,
#'   r = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
#'   ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
#'   nMC = 500L, method = 'MC',
#'   gamma1_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma2_grid = seq(0.01, 0.99, by = 0.01),
#'   seed = 1L
#' )
#' }
#'
#' # Example 3: Uncontrolled design, posterior probability, vague prior
#' \dontrun{
#' Sigma <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' getgamma2cont(
#'   nsim = 1000L, prob = 'posterior', design = 'uncontrolled',
#'   prior = 'vague',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   mu1 = c(-5.0, 0.0), Sigma1 = Sigma,
#'   mu2 = NULL, Sigma2 = NULL,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n1 = 30L, n2 = NULL,
#'   theta.TV1 = 10.0, theta.MAV1 = 5.0,
#'   theta.TV2 = 2.0,  theta.MAV2 = 1.0,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   m1 = NULL, m2 = NULL,
#'   kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
#'   kappa02 = NULL, nu02 = NULL, mu02 = c(-10.0, -1.0), Lambda02 = NULL,
#'   r = 1.0,
#'   ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
#'   ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
#'   nMC = NULL, method = 'MM',
#'   gamma1_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma2_grid = seq(0.01, 0.99, by = 0.01),
#'   seed = 1L
#' )
#' }
#'
#' # Example 4: Controlled design, predictive probability, vague prior
#' # Note: prob = 'predictive' returns 4 regions; GoRegions = 1L (both
#' # endpoints exceed NULL), NoGoRegions = 4L (both endpoints below NULL)
#' \dontrun{
#' Sigma <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' getgamma2cont(
#'   nsim = 1000L, prob = 'predictive', design = 'controlled',
#'   prior = 'vague',
#'   GoRegions = 1L, NoGoRegions = 4L,
#'   mu1 = c(-5.0, 0.0), Sigma1 = Sigma,
#'   mu2 = c(-10.0, -1.0), Sigma2 = Sigma,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n1 = 30L, n2 = 30L,
#'   theta.TV1 = NULL, theta.MAV1 = NULL,
#'   theta.TV2 = NULL, theta.MAV2 = NULL,
#'   theta.NULL1 = 5.0, theta.NULL2 = 1.0,
#'   m1 = 100L, m2 = 100L,
#'   kappa01 = NULL, nu01 = NULL, mu01 = NULL, Lambda01 = NULL,
#'   kappa02 = NULL, nu02 = NULL, mu02 = NULL, Lambda02 = NULL,
#'   r = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01e = NULL, alpha02e = NULL,
#'   ybar_e1 = NULL, ybar_e2 = NULL, Se1 = NULL, Se2 = NULL,
#'   nMC = 500L, method = 'MC',
#'   gamma1_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma2_grid = seq(0.01, 0.99, by = 0.01),
#'   seed = 1L
#' )
#' }
#'
#' @importFrom stats rnorm
#' @export
getgamma2cont <- function(nsim        = 10000L,
                          prob        = 'posterior',
                          design      = 'controlled',
                          prior       = 'vague',
                          GoRegions, NoGoRegions,
                          mu1, Sigma1,
                          mu2    = NULL, Sigma2 = NULL,
                          target_go, target_nogo,
                          crit_go  = '<',        crit_nogo  = '<',
                          sel_go   = 'smallest', sel_nogo   = 'largest',
                          n1, n2 = NULL,
                          theta.TV1   = NULL, theta.MAV1  = NULL,
                          theta.TV2   = NULL, theta.MAV2  = NULL,
                          theta.NULL1 = NULL, theta.NULL2 = NULL,
                          m1 = NULL, m2 = NULL,
                          kappa01  = NULL, nu01     = NULL,
                          mu01     = NULL, Lambda01 = NULL,
                          kappa02  = NULL, nu02     = NULL,
                          mu02     = NULL, Lambda02 = NULL,
                          r        = NULL,
                          ne1 = NULL, ne2 = NULL,
                          alpha01e = NULL, alpha02e = NULL,
                          ybar_e1  = NULL, ybar_e2  = NULL,
                          Se1      = NULL, Se2      = NULL,
                          nMC    = NULL,
                          method = 'MC',
                          gamma1_grid = seq(0.01, 0.99, by = 0.01),
                          gamma2_grid = seq(0.01, 0.99, by = 0.01),
                          seed = NULL) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) ||
      nsim != floor(nsim) || nsim < 1L) {
    stop("'nsim' must be a single positive integer")
  }
  nsim <- as.integer(nsim)

  if (!is.character(prob) || length(prob) != 1L ||
      !prob %in% c('posterior', 'predictive')) {
    stop("'prob' must be either 'posterior' or 'predictive'")
  }

  if (!is.character(design) || length(design) != 1L ||
      !design %in% c('controlled', 'uncontrolled', 'external')) {
    stop("'design' must be 'controlled', 'uncontrolled', or 'external'")
  }

  if (!is.character(prior) || length(prior) != 1L ||
      !prior %in% c('vague', 'N-Inv-Wishart')) {
    stop("'prior' must be 'vague' or 'N-Inv-Wishart'")
  }

  if (!is.integer(GoRegions) || length(GoRegions) < 1L ||
      any(is.na(GoRegions)) || any(GoRegions < 1L) || any(GoRegions > 9L)) {
    stop("'GoRegions' must be an integer vector with all values in 1:9")
  }
  if (!is.integer(NoGoRegions) || length(NoGoRegions) < 1L ||
      any(is.na(NoGoRegions)) || any(NoGoRegions < 1L) ||
      any(NoGoRegions > 9L)) {
    stop("'NoGoRegions' must be an integer vector with all values in 1:9")
  }
  if (length(intersect(GoRegions, NoGoRegions)) > 0L) {
    stop("'GoRegions' and 'NoGoRegions' must be disjoint")
  }

  if (!is.numeric(mu1) || length(mu1) != 2L || any(is.na(mu1))) {
    stop("'mu1' must be a length-2 numeric vector")
  }
  if (!is.matrix(Sigma1) || !all(dim(Sigma1) == c(2L, 2L)) ||
      !is.numeric(Sigma1) || any(is.na(Sigma1))) {
    stop("'Sigma1' must be a 2x2 numeric matrix")
  }

  if (design != 'uncontrolled') {
    if (is.null(mu2) || !is.numeric(mu2) || length(mu2) != 2L ||
        any(is.na(mu2))) {
      stop("'mu2' must be a length-2 numeric vector for controlled or external design")
    }
    if (is.null(Sigma2) || !is.matrix(Sigma2) ||
        !all(dim(Sigma2) == c(2L, 2L)) || !is.numeric(Sigma2) ||
        any(is.na(Sigma2))) {
      stop("'Sigma2' must be a 2x2 numeric matrix for controlled or external design")
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
    if (is.null(n2) || !is.numeric(n2) || length(n2) != 1L || is.na(n2) ||
        n2 != floor(n2) || n2 < 1L) {
      stop("'n2' must be a single positive integer for controlled or external design")
    }
  }

  if (prob == 'posterior') {
    for (nm in c("theta.TV1", "theta.MAV1", "theta.TV2", "theta.MAV2")) {
      if (is.null(get(nm))) {
        stop(paste0("'", nm, "' must be non-NULL when prob = 'posterior'"))
      }
    }
    if (theta.TV1 <= theta.MAV1) {
      stop("'theta.TV1' must be strictly greater than 'theta.MAV1'")
    }
    if (theta.TV2 <= theta.MAV2) {
      stop("'theta.TV2' must be strictly greater than 'theta.MAV2'")
    }
  } else {
    for (nm in c("theta.NULL1", "theta.NULL2")) {
      if (is.null(get(nm))) {
        stop(paste0("'", nm, "' must be non-NULL when prob = 'predictive'"))
      }
    }
    if (is.null(m1) || !is.numeric(m1) || length(m1) != 1L || is.na(m1) ||
        m1 != floor(m1) || m1 < 1L) {
      stop("'m1' must be a single positive integer when prob = 'predictive'")
    }
    if (is.null(m2) || !is.numeric(m2) || length(m2) != 1L || is.na(m2) ||
        m2 != floor(m2) || m2 < 1L) {
      stop("'m2' must be a single positive integer when prob = 'predictive'")
    }
  }

  if (!is.character(method) || length(method) != 1L ||
      !method %in% c('MC', 'MM')) {
    stop("'method' must be either 'MC' or 'MM'")
  }
  if (method == 'MC') {
    if (is.null(nMC) || !is.numeric(nMC) || length(nMC) != 1L ||
        is.na(nMC) || nMC != floor(nMC) || nMC < 1L) {
      stop("'nMC' must be a single positive integer when method = 'MC'")
    }
  }

  for (gname in c("gamma1_grid", "gamma2_grid")) {
    gval <- get(gname)
    if (!is.numeric(gval) || length(gval) < 1L ||
        any(is.na(gval)) || any(gval <= 0) || any(gval >= 1)) {
      stop(paste0("'", gname, "' must be a numeric vector with all values in (0, 1)"))
    }
  }
  gamma1_grid <- sort(unique(gamma1_grid))
  gamma2_grid <- sort(unique(gamma2_grid))
  ng1         <- length(gamma1_grid)
  ng2         <- length(gamma2_grid)

  # ---------------------------------------------------------------------------
  # Stage 1: Simulate nsim datasets and precompute PrGo_vec, PrNoGo_vec
  #          via pbayespostpred2cont (vectorised call)
  # ---------------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  n1 <- as.integer(n1)

  # Generate nsim * n1 bivariate normal residuals (scenario-invariant)
  R_Sigma1   <- chol(Sigma1)
  Z1_raw     <- matrix(rnorm(nsim * n1 * 2L),
                       nrow = nsim * n1, ncol = 2L) %*% R_Sigma1
  block1     <- rep(seq_len(nsim), each = n1)
  Z1_colsums <- apply(Z1_raw, 2L, function(col) tapply(col, block1, sum))
  # Z1_colsums: nsim x 2 matrix of column sums

  # Scatter matrices for treatment arm
  Z1_colmeans_rep <- Z1_colsums[block1, ] / n1
  Z1_centered     <- Z1_raw - Z1_colmeans_rep
  S1_11 <- tapply(Z1_centered[, 1L] ^ 2,                 block1, sum)
  S1_12 <- tapply(Z1_centered[, 1L] * Z1_centered[, 2L], block1, sum)
  S1_22 <- tapply(Z1_centered[, 2L] ^ 2,                 block1, sum)

  # Shift residuals by the (single) scenario mean: nsim x 2
  ybar1_sim <- sweep(Z1_colsums / n1, 2L, mu1, '+')

  # Build S1 list for vectorised call
  S1_list <- vector('list', nsim)
  for (i in seq_len(nsim)) {
    S1_list[[i]] <- matrix(c(S1_11[i], S1_12[i], S1_12[i], S1_22[i]),
                           nrow = 2L, ncol = 2L)
  }

  if (design %in% c('controlled', 'external')) {
    n2 <- as.integer(n2)
    R_Sigma2   <- chol(Sigma2)
    Z2_raw     <- matrix(rnorm(nsim * n2 * 2L),
                         nrow = nsim * n2, ncol = 2L) %*% R_Sigma2
    block2     <- rep(seq_len(nsim), each = n2)
    Z2_colsums <- apply(Z2_raw, 2L, function(col) tapply(col, block2, sum))
    Z2_colmeans_rep <- Z2_colsums[block2, ] / n2
    Z2_centered     <- Z2_raw - Z2_colmeans_rep
    S2_11 <- tapply(Z2_centered[, 1L] ^ 2,                 block2, sum)
    S2_12 <- tapply(Z2_centered[, 1L] * Z2_centered[, 2L], block2, sum)
    S2_22 <- tapply(Z2_centered[, 2L] ^ 2,                 block2, sum)
    ybar2_sim <- sweep(Z2_colsums / n2, 2L, mu2, '+')
    S2_list   <- vector('list', nsim)
    for (i in seq_len(nsim)) {
      S2_list[[i]] <- matrix(c(S2_11[i], S2_12[i], S2_12[i], S2_22[i]),
                             nrow = 2L, ncol = 2L)
    }
  } else {
    ybar2_sim <- NULL
    S2_list   <- NULL
  }

  # Vectorised call to pbayespostpred2cont: returns nsim x n_regions matrix
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
  # Pr_R_mat: nsim x 9 matrix

  # Go/NoGo probability per replicate (length nsim vectors)
  PrGo_vec   <- rowSums(Pr_R_mat[, GoRegions,   drop = FALSE])
  PrNoGo_vec <- rowSums(Pr_R_mat[, NoGoRegions, drop = FALSE])

  # ---------------------------------------------------------------------------
  # Stage 2: Sweep (gamma1_grid x gamma2_grid)
  #
  # For each (g1, g2):
  #   Pr(Go)   = mean( I(PrGo_vec >= g1  AND PrNoGo_vec <  g2) )
  #   Pr(NoGo) = mean( I(PrNoGo_vec >= g2 AND PrGo_vec  <  g1) )
  # ---------------------------------------------------------------------------
  PrGo_grid   <- matrix(NA_real_, nrow = ng1, ncol = ng2)
  PrNoGo_grid <- matrix(NA_real_, nrow = ng1, ncol = ng2)

  for (k1 in seq_len(ng1)) {
    g1 <- gamma1_grid[k1]
    for (k2 in seq_len(ng2)) {
      g2 <- gamma2_grid[k2]

      go_mask   <- (PrGo_vec   >= g1) & (PrNoGo_vec <  g2)
      nogo_mask <- (PrNoGo_vec >= g2) & (PrGo_vec   <  g1)

      PrGo_grid[k1, k2]   <- mean(go_mask)
      PrNoGo_grid[k1, k2] <- mean(nogo_mask)
    }
  }

  # ---------------------------------------------------------------------------
  # Stage 3: Select optimal (gamma1, gamma2)
  #
  # gamma1: worst-case (max) Pr(Go) over gamma2_grid for each gamma1 candidate
  # gamma2: worst-case (max) Pr(NoGo) over gamma1_grid for each gamma2 candidate
  # ---------------------------------------------------------------------------
  .compare <- function(x, op, val) {
    switch(op,
           '<'  = x <  val,
           '<=' = x <= val,
           '>'  = x >  val,
           '>=' = x >= val)
  }

  .select_idx <- function(mask, sel) {
    idx <- which(mask)
    if (length(idx) == 0L) return(NA_integer_)
    if (sel == 'smallest') min(idx) else max(idx)
  }

  max_PrGo_per_g1   <- apply(PrGo_grid,   1L, max, na.rm = TRUE)
  opt1 <- .select_idx(.compare(max_PrGo_per_g1,   crit_go,   target_go),   sel_go)

  max_PrNoGo_per_g2 <- apply(PrNoGo_grid, 2L, max, na.rm = TRUE)
  opt2 <- .select_idx(.compare(max_PrNoGo_per_g2, crit_nogo, target_nogo), sel_nogo)

  if (is.na(opt1) || is.na(opt2)) {
    gamma1          <- if (!is.na(opt1)) gamma1_grid[opt1] else NA_real_
    gamma2          <- if (!is.na(opt2)) gamma2_grid[opt2] else NA_real_
    PrGo_at_gamma   <- NA_real_
    PrNoGo_at_gamma <- NA_real_
  } else {
    gamma1          <- gamma1_grid[opt1]
    gamma2          <- gamma2_grid[opt2]
    PrGo_at_gamma   <- PrGo_grid[opt1, opt2]
    PrNoGo_at_gamma <- PrNoGo_grid[opt1, opt2]
  }

  # ---------------------------------------------------------------------------
  # Build and return result
  # ---------------------------------------------------------------------------
  result <- list(
    gamma1          = gamma1,
    gamma2          = gamma2,
    PrGo_at_gamma   = PrGo_at_gamma,
    PrNoGo_at_gamma = PrNoGo_at_gamma,
    gamma1_grid     = gamma1_grid,
    gamma2_grid     = gamma2_grid,
    PrGo_grid       = PrGo_grid,
    PrNoGo_grid     = PrNoGo_grid
  )

  class(result) <- 'getgamma2cont'
  return(result)
}
