#' Go/NoGo/Gray Decision Probabilities for a Clinical Trial with Two Binary
#' Endpoints
#'
#' Evaluates operating characteristics (Go, NoGo, Gray probabilities) for
#' clinical trials with two binary endpoints under the Bayesian framework.
#' The function supports controlled, uncontrolled, and external-control designs,
#' and uses both posterior probability and posterior predictive probability
#' criteria.
#'
#' Computation proceeds in two stages for efficiency.  In Stage 1, posterior
#' or predictive probabilities are precomputed for every possible multinomial
#' count combination \eqn{(x_t, x_c)}.  In Stage 2, operating characteristics
#' for each scenario are obtained by weighting the Stage 1 decisions by their
#' multinomial probabilities under the true parameters, avoiding repeated Monte
#' Carlo sampling per scenario.
#'
#' @param prob A character string specifying the probability type for
#'        decision-making.  Must be \code{'posterior'} or
#'        \code{'predictive'}.
#' @param design A character string specifying the trial design.  Must be
#'        \code{'controlled'}, \code{'uncontrolled'}, or \code{'external'}.
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
#' @param gamma1 A numeric scalar in \code{(0, 1)} giving the minimum
#'        posterior/predictive probability required for a Go decision.
#' @param gamma2 A numeric scalar in \code{(0, 1)} giving the minimum
#'        posterior/predictive probability required for a NoGo decision.
#'        Unlike single-endpoint designs, \code{gamma2} may be greater than,
#'        equal to, or less than \code{gamma1}, because the Go and NoGo
#'        regions are structurally asymmetric (e.g., Go = R1 only vs
#'        NoGo = R9 only) and their calibrated thresholds are independent.
#' @param pi_t1 A numeric vector of true treatment response probabilities
#'        for Endpoint 1.  Each element must be in \code{(0, 1)}.
#' @param pi_t2 A numeric vector of true treatment response probabilities
#'        for Endpoint 2.  Must have the same length as \code{pi_t1}.
#' @param rho_t A numeric vector of true within-treatment-arm correlations
#'        between Endpoint 1 and Endpoint 2.  Must have the same length as
#'        \code{pi_t1}.  Each element must be within the feasible range
#'        for the corresponding \code{(pi_t1, pi_t2)} pair (checked via
#'        \code{\link{getjointbin}}).
#' @param pi_c1 A numeric vector of true control response probabilities for
#'        Endpoint 1.  For \code{design = 'uncontrolled'}, this parameter
#'        is not used in probability calculations but must still be supplied;
#'        it is retained in the output for reference.  Must have the same
#'        length as \code{pi_t1}.
#' @param pi_c2 A numeric vector of true control response probabilities for
#'        Endpoint 2.  For \code{design = 'uncontrolled'}, see note for
#'        \code{pi_c1}.  Must have the same length as \code{pi_t1}.
#' @param rho_c A numeric vector of true within-control-arm correlations.
#'        For \code{design = 'uncontrolled'}, not used in probability
#'        calculations.  Must have the same length as \code{pi_t1}.
#' @param n1 A positive integer giving the number of patients in group 1
#'        (treatment) in the proof-of-concept (PoC) trial.
#' @param n2 A positive integer giving the number of patients in group 2
#'        (control) in the PoC trial.
#' @param a1_00 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (0,0) response pattern in group 1.
#' @param a1_01 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (0,1) response pattern in group 1.
#' @param a1_10 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (1,0) response pattern in group 1.
#' @param a1_11 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (1,1) response pattern in group 1.
#' @param a2_00 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (0,0) response pattern in group 2.  For
#'        \code{design = 'uncontrolled'}, serves as a hyperparameter of
#'        the hypothetical control distribution.
#' @param a2_01 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (0,1) response pattern in group 2.
#' @param a2_10 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (1,0) response pattern in group 2.
#' @param a2_11 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (1,1) response pattern in group 2.
#' @param m1 A positive integer giving the future sample size for group 1.
#'        Required when \code{prob = 'predictive'}; otherwise set to
#'        \code{NULL}.
#' @param m2 A positive integer giving the future sample size for group 2.
#'        Required when \code{prob = 'predictive'}; otherwise set to
#'        \code{NULL}.
#' @param theta.TV1 A numeric scalar giving the TV threshold for Endpoint 1.
#'        Required when \code{prob = 'posterior'}; otherwise set to
#'        \code{NULL}.
#' @param theta.MAV1 A numeric scalar giving the MAV threshold for Endpoint
#'        1.  Required when \code{prob = 'posterior'}; otherwise set to
#'        \code{NULL}.
#' @param theta.TV2 A numeric scalar giving the TV threshold for Endpoint 2.
#'        Required when \code{prob = 'posterior'}; otherwise set to
#'        \code{NULL}.
#' @param theta.MAV2 A numeric scalar giving the MAV threshold for Endpoint
#'        2.  Required when \code{prob = 'posterior'}; otherwise set to
#'        \code{NULL}.
#' @param theta.NULL1 A numeric scalar giving the null hypothesis threshold
#'        for Endpoint 1.  Required when \code{prob = 'predictive'};
#'        otherwise set to \code{NULL}.
#' @param theta.NULL2 A numeric scalar giving the null hypothesis threshold
#'        for Endpoint 2.  Required when \code{prob = 'predictive'};
#'        otherwise set to \code{NULL}.
#' @param z00 A non-negative integer giving the hypothetical control count
#'        for pattern (0,0).  Required when \code{design = 'uncontrolled'};
#'        otherwise set to \code{NULL}.
#' @param z01 A non-negative integer giving the hypothetical control count
#'        for pattern (0,1).  Required when \code{design = 'uncontrolled'};
#'        otherwise set to \code{NULL}.
#' @param z10 A non-negative integer giving the hypothetical control count
#'        for pattern (1,0).  Required when \code{design = 'uncontrolled'};
#'        otherwise set to \code{NULL}.
#' @param z11 A non-negative integer giving the hypothetical control count
#'        for pattern (1,1).  Required when \code{design = 'uncontrolled'};
#'        otherwise set to \code{NULL}.
#' @param xe1_00 A non-negative integer giving the external group 1 count
#'        for pattern (0,0).  Required when \code{design = 'external'} and
#'        external treatment data are used; otherwise \code{NULL}.
#' @param xe1_01 A non-negative integer; see \code{xe1_00}.
#' @param xe1_10 A non-negative integer; see \code{xe1_00}.
#' @param xe1_11 A non-negative integer; see \code{xe1_00}.
#' @param xe2_00 A non-negative integer giving the external group 2 count
#'        for pattern (0,0).  Required when \code{design = 'external'} and
#'        external control data are used; otherwise \code{NULL}.
#' @param xe2_01 A non-negative integer; see \code{xe2_00}.
#' @param xe2_10 A non-negative integer; see \code{xe2_00}.
#' @param xe2_11 A non-negative integer; see \code{xe2_00}.
#' @param ae1 A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for group 1.  Required when external treatment data are
#'        used; otherwise \code{NULL}.
#' @param ae2 A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for group 2.  Required when external control data are
#'        used; otherwise \code{NULL}.
#' @param nMC A positive integer giving the number of Monte Carlo draws
#'        used inside \code{\link{pbayespostpred2bin}} for each count combination.
#'        Default is \code{10000}.
#' @param error_if_Miss A logical scalar; if \code{TRUE} (default), the
#'        function stops with an error if the Miss probability is positive,
#'        prompting reconsideration of the thresholds.
#' @param Gray_inc_Miss A logical scalar; if \code{TRUE}, the Miss
#'        probability is added to Gray.  If \code{FALSE} (default), Miss is
#'        reported separately.  Active only when
#'        \code{error_if_Miss = FALSE}.
#'
#' @return A data frame with one row per scenario and columns:
#' \describe{
#'   \item{pi_t1}{True treatment response probability for Endpoint 1.}
#'   \item{pi_t2}{True treatment response probability for Endpoint 2.}
#'   \item{rho_t}{True within-arm correlation in the treatment arm.}
#'   \item{pi_c1}{True control response probability for Endpoint 1
#'                (omitted when \code{design = 'uncontrolled'}).}
#'   \item{pi_c2}{True control response probability for Endpoint 2
#'                (omitted when \code{design = 'uncontrolled'}).}
#'   \item{rho_c}{True within-arm correlation in the control arm
#'                (omitted when \code{design = 'uncontrolled'}).}
#'   \item{Go}{Probability of making a Go decision.}
#'   \item{Gray}{Probability of making a Gray (inconclusive) decision.}
#'   \item{NoGo}{Probability of making a NoGo decision.}
#'   \item{Miss}{Probability where Go and NoGo criteria are simultaneously
#'               met.  Included when \code{error_if_Miss = FALSE} and
#'               \code{Gray_inc_Miss = FALSE}.}
#' }
#' The returned object has S3 class \code{pbayesdecisionprob2bin} with an associated
#' \code{print} method.
#'
#' @details
#' \strong{Two-stage computation.}
#'
#' Stage 1 (precomputation): All possible multinomial count combinations
#' \eqn{(x_t, x_c)} are enumerated using \code{\link{allmultinom}}.
#' For each combination, \code{\link{pbayespostpred2bin}} computes the Go and NoGo
#' posterior/predictive probabilities via Monte Carlo, yielding matrices
#' \code{PrGo[i, j]} and \code{PrNoGo[i, j]}.  This stage is independent of
#' the true scenario parameters.
#'
#' Stage 2 (scenario evaluation): For each scenario
#' \eqn{(\pi_{t1}, \pi_{t2}, \rho_t, \pi_{c1}, \pi_{c2}, \rho_c)},
#' \code{\link{getjointbin}} converts the marginal rates and correlation
#' into the four-cell probability vector \eqn{p_k}.  The multinomial
#' probability mass function weights the Stage 1 decisions:
#' \deqn{\Pr(\mathrm{Go} \mid \Theta) = \sum_{i,j}
#'       \mathbf{1}(\mathrm{Go}_{ij}) \cdot
#'       P(x_{t,i} \mid n_1, p_t) \cdot P(x_{c,j} \mid n_2, p_c),}
#' where \eqn{\mathbf{1}(\mathrm{Go}_{ij})} is the Go indicator from Stage 1.
#' Stage 2 requires no additional Monte Carlo and is fast even for large
#' scenario grids.
#'
#' \strong{Decision categories.}
#' \itemize{
#'   \item \strong{Go}: sum of Go region probabilities \eqn{\ge \gamma_1}
#'         AND sum of NoGo region probabilities \eqn{< \gamma_2}.
#'   \item \strong{NoGo}: sum of Go region probabilities \eqn{< \gamma_1}
#'         AND sum of NoGo region probabilities \eqn{\ge \gamma_2}.
#'   \item \strong{Miss}: both Go and NoGo criteria satisfied simultaneously.
#'   \item \strong{Gray}: neither Go nor NoGo criterion satisfied.
#' }
#'
#' @examples
#' # Example 1: Posterior probability, controlled design (rho = 0)
#' # Corresponds to Table 3(a) of the reference with independent endpoints.
#' pbayesdecisionprob2bin(
#'   prob        = 'posterior',
#'   design      = 'controlled',
#'   GoRegions   = 1L,
#'   NoGoRegions = 9L,
#'   gamma1      = 0.27,
#'   gamma2      = 0.36,
#'   pi_t1       = c(0.15, 0.25, 0.30, 0.40),
#'   pi_t2       = c(0.20, 0.30, 0.35, 0.45),
#'   rho_t       = rep(0.0, 4),
#'   pi_c1       = rep(0.15, 4),
#'   pi_c2       = rep(0.20, 4),
#'   rho_c       = rep(0.0, 4),
#'   n1 = 10, n2 = 10,
#'   a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
#'   a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
#'   m1 = NULL, m2 = NULL,
#'   theta.TV1   = 0.15, theta.MAV1 = 0.10,
#'   theta.TV2   = 0.15, theta.MAV2 = 0.10,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 100,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 2: Posterior probability, controlled design (rho > 0)
#' # Explores the impact of positive within-arm endpoint correlation.
#' pbayesdecisionprob2bin(
#'   prob        = 'posterior',
#'   design      = 'controlled',
#'   GoRegions   = 1L,
#'   NoGoRegions = 9L,
#'   gamma1      = 0.52,
#'   gamma2      = 0.52,
#'   pi_t1       = c(0.15, 0.25, 0.30, 0.40),
#'   pi_t2       = c(0.20, 0.30, 0.35, 0.45),
#'   rho_t       = rep(0.6, 4),
#'   pi_c1       = rep(0.15, 4),
#'   pi_c2       = rep(0.20, 4),
#'   rho_c       = rep(0.6, 4),
#'   n1 = 10, n2 = 10,
#'   a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
#'   a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
#'   m1 = NULL, m2 = NULL,
#'   theta.TV1   = 0.15, theta.MAV1 = 0.10,
#'   theta.TV2   = 0.15, theta.MAV2 = 0.10,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 100,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 3: Posterior probability, uncontrolled design
#' # Hypothetical control specified via prior (a2_*) and pseudo-counts (z*).
#' pbayesdecisionprob2bin(
#'   prob        = 'posterior',
#'   design      = 'uncontrolled',
#'   GoRegions   = 1L,
#'   NoGoRegions = 9L,
#'   gamma1      = 0.27,
#'   gamma2      = 0.36,
#'   pi_t1       = c(0.15, 0.30, 0.40),
#'   pi_t2       = c(0.20, 0.35, 0.45),
#'   rho_t       = rep(0.0, 3),
#'   pi_c1       = rep(0.15, 3),
#'   pi_c2       = rep(0.20, 3),
#'   rho_c       = rep(0.0, 3),
#'   n1 = 10, n2 = 10,
#'   a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
#'   a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
#'   m1 = NULL, m2 = NULL,
#'   theta.TV1   = 0.15, theta.MAV1 = 0.10,
#'   theta.TV2   = 0.15, theta.MAV2 = 0.10,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   z00 = 2L, z01 = 1L, z10 = 2L, z11 = 1L,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 100,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 4: Posterior probability, external control design
#' # External control data incorporated via power prior (ae2 = 0.5).
#' pbayesdecisionprob2bin(
#'   prob        = 'posterior',
#'   design      = 'external',
#'   GoRegions   = 1L,
#'   NoGoRegions = 9L,
#'   gamma1      = 0.27,
#'   gamma2      = 0.36,
#'   pi_t1       = c(0.15, 0.25, 0.30, 0.40),
#'   pi_t2       = c(0.20, 0.30, 0.35, 0.45),
#'   rho_t       = rep(0.0, 4),
#'   pi_c1       = rep(0.15, 4),
#'   pi_c2       = rep(0.20, 4),
#'   rho_c       = rep(0.0, 4),
#'   n1 = 10, n2 = 10,
#'   a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
#'   a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
#'   m1 = NULL, m2 = NULL,
#'   theta.TV1   = 0.15, theta.MAV1 = 0.10,
#'   theta.TV2   = 0.15, theta.MAV2 = 0.10,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = 4L,  xe2_01 = 2L,  xe2_10 = 3L,  xe2_11 = 1L,
#'   ae1 = NULL, ae2 = 0.5,
#'   nMC = 100,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 5: Predictive probability, controlled design
#' # Future trial has m1 = m2 = 30 patients per arm.
#' pbayesdecisionprob2bin(
#'   prob        = 'predictive',
#'   design      = 'controlled',
#'   GoRegions   = 1L,
#'   NoGoRegions = 4L,
#'   gamma1      = 0.60,
#'   gamma2      = 0.80,
#'   pi_t1       = c(0.15, 0.25, 0.30, 0.40),
#'   pi_t2       = c(0.20, 0.30, 0.35, 0.45),
#'   rho_t       = rep(0.0, 4),
#'   pi_c1       = rep(0.15, 4),
#'   pi_c2       = rep(0.20, 4),
#'   rho_c       = rep(0.0, 4),
#'   n1 = 10, n2 = 10,
#'   a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
#'   a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
#'   m1 = 30L, m2 = 30L,
#'   theta.TV1   = NULL, theta.MAV1 = NULL,
#'   theta.TV2   = NULL, theta.MAV2 = NULL,
#'   theta.NULL1 = 0.10, theta.NULL2 = 0.10,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 100,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 6: Predictive probability, uncontrolled design
#' # Hypothetical control specified via pseudo-counts; future trial m1 = m2 = 30.
#' pbayesdecisionprob2bin(
#'   prob        = 'predictive',
#'   design      = 'uncontrolled',
#'   GoRegions   = 1L,
#'   NoGoRegions = 4L,
#'   gamma1      = 0.60,
#'   gamma2      = 0.80,
#'   pi_t1       = c(0.15, 0.30, 0.40),
#'   pi_t2       = c(0.20, 0.35, 0.45),
#'   rho_t       = rep(0.0, 3),
#'   pi_c1       = rep(0.15, 3),
#'   pi_c2       = rep(0.20, 3),
#'   rho_c       = rep(0.0, 3),
#'   n1 = 10, n2 = 10,
#'   a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
#'   a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
#'   m1 = 30L, m2 = 30L,
#'   theta.TV1   = NULL, theta.MAV1 = NULL,
#'   theta.TV2   = NULL, theta.MAV2 = NULL,
#'   theta.NULL1 = 0.10, theta.NULL2 = 0.10,
#'   z00 = 2L, z01 = 1L, z10 = 2L, z11 = 1L,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 100,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 7: Predictive probability, external control design
#' # External treatment data incorporated via power prior (ae1 = 0.5).
#' pbayesdecisionprob2bin(
#'   prob        = 'predictive',
#'   design      = 'external',
#'   GoRegions   = 1L,
#'   NoGoRegions = 4L,
#'   gamma1      = 0.60,
#'   gamma2      = 0.80,
#'   pi_t1       = c(0.15, 0.25, 0.30, 0.40),
#'   pi_t2       = c(0.20, 0.30, 0.35, 0.45),
#'   rho_t       = rep(0.0, 4),
#'   pi_c1       = rep(0.15, 4),
#'   pi_c2       = rep(0.20, 4),
#'   rho_c       = rep(0.0, 4),
#'   n1 = 10, n2 = 10,
#'   a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
#'   a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
#'   m1 = 30L, m2 = 30L,
#'   theta.TV1   = NULL, theta.MAV1 = NULL,
#'   theta.TV2   = NULL, theta.MAV2 = NULL,
#'   theta.NULL1 = 0.10, theta.NULL2 = 0.10,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe1_00 = 3L, xe1_01 = 2L, xe1_10 = 3L, xe1_11 = 2L,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = 0.5, ae2 = NULL,
#'   nMC = 100,
#'   error_if_Miss = TRUE, Gray_inc_Miss = FALSE
#' )
#'
#' @importFrom stats dmultinom
#' @export
pbayesdecisionprob2bin <- function(prob        = 'posterior',
                                   design      = 'controlled',
                                   GoRegions,
                                   NoGoRegions,
                                   gamma1,
                                   gamma2,
                                   pi_t1, pi_t2, rho_t,
                                   pi_c1, pi_c2, rho_c,
                                   n1, n2,
                                   a1_00, a1_01, a1_10, a1_11,
                                   a2_00, a2_01, a2_10, a2_11,
                                   m1          = NULL,
                                   m2          = NULL,
                                   theta.TV1   = NULL, theta.MAV1  = NULL,
                                   theta.TV2   = NULL, theta.MAV2  = NULL,
                                   theta.NULL1 = NULL, theta.NULL2 = NULL,
                                   z00 = NULL, z01 = NULL,
                                   z10 = NULL, z11 = NULL,
                                   xe1_00 = NULL, xe1_01 = NULL,
                                   xe1_10 = NULL, xe1_11 = NULL,
                                   xe2_00 = NULL, xe2_01 = NULL,
                                   xe2_10 = NULL, xe2_11 = NULL,
                                   ae1         = NULL,
                                   ae2         = NULL,
                                   nMC         = 10000L,
                                   error_if_Miss = TRUE,
                                   Gray_inc_Miss = FALSE) {

  # ---------------------------------------------------------------------------
  # Section 1: Input validation
  # ---------------------------------------------------------------------------

  if (!is.character(prob) || length(prob) != 1L ||
      !prob %in% c('posterior', 'predictive')) {
    stop("'prob' must be either 'posterior' or 'predictive'")
  }

  if (!is.character(design) || length(design) != 1L ||
      !design %in% c('controlled', 'uncontrolled', 'external')) {
    stop("'design' must be 'controlled', 'uncontrolled', or 'external'")
  }

  # Determine valid region range based on prob type
  max_region <- if (prob == 'posterior') 9L else 4L

  if (!is.numeric(GoRegions) || length(GoRegions) < 1L ||
      any(is.na(GoRegions)) || any(GoRegions != floor(GoRegions)) ||
      any(GoRegions < 1L) || any(GoRegions > max_region)) {
    stop(sprintf("'GoRegions' must be an integer vector with values in 1:%d",
                 max_region))
  }

  if (!is.numeric(NoGoRegions) || length(NoGoRegions) < 1L ||
      any(is.na(NoGoRegions)) || any(NoGoRegions != floor(NoGoRegions)) ||
      any(NoGoRegions < 1L) || any(NoGoRegions > max_region)) {
    stop(sprintf("'NoGoRegions' must be an integer vector with values in 1:%d",
                 max_region))
  }

  if (length(intersect(GoRegions, NoGoRegions)) > 0L) {
    stop("'GoRegions' and 'NoGoRegions' must be disjoint")
  }

  if (!is.numeric(gamma1) || length(gamma1) != 1L || is.na(gamma1) ||
      gamma1 <= 0 || gamma1 >= 1) {
    stop("'gamma1' must be a single numeric value in (0, 1)")
  }

  if (!is.numeric(gamma2) || length(gamma2) != 1L || is.na(gamma2) ||
      gamma2 <= 0 || gamma2 >= 1) {
    stop("'gamma2' must be a single numeric value in (0, 1)")
  }

  # Note: gamma2 >= gamma1 is allowed for two-endpoint designs because
  # GoRegions and NoGoRegions are asymmetric and their thresholds are
  # calibrated independently (see Table 3 of the reference).

  # --- Scenario vectors ---
  pi_t1 <- as.numeric(pi_t1)
  pi_t2 <- as.numeric(pi_t2)
  rho_t  <- as.numeric(rho_t)
  pi_c1 <- as.numeric(pi_c1)
  pi_c2 <- as.numeric(pi_c2)
  rho_c  <- as.numeric(rho_c)

  n_scen <- length(pi_t1)

  if (any(is.na(pi_t1)) || any(pi_t1 <= 0) || any(pi_t1 >= 1))
    stop("All elements of 'pi_t1' must be in (0, 1)")
  if (length(pi_t2) != n_scen || any(is.na(pi_t2)) ||
      any(pi_t2 <= 0) || any(pi_t2 >= 1))
    stop("'pi_t2' must have the same length as 'pi_t1' with all elements in (0, 1)")
  if (length(rho_t) != n_scen || any(is.na(rho_t)))
    stop("'rho_t' must have the same length as 'pi_t1'")
  if (length(pi_c1) != n_scen || any(is.na(pi_c1)) ||
      any(pi_c1 <= 0) || any(pi_c1 >= 1))
    stop("'pi_c1' must have the same length as 'pi_t1' with all elements in (0, 1)")
  if (length(pi_c2) != n_scen || any(is.na(pi_c2)) ||
      any(pi_c2 <= 0) || any(pi_c2 >= 1))
    stop("'pi_c2' must have the same length as 'pi_t1' with all elements in (0, 1)")
  if (length(rho_c) != n_scen || any(is.na(rho_c)))
    stop("'rho_c' must have the same length as 'pi_t1'")

  # --- Sample sizes ---
  for (nm in c("n1", "n2")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val != floor(val) || val < 1L)
      stop(paste0("'", nm, "' must be a single positive integer"))
  }
  n1 <- as.integer(n1)
  n2 <- as.integer(n2)

  # --- Dirichlet prior parameters ---
  for (nm in c("a1_00", "a1_01", "a1_10", "a1_11",
               "a2_00", "a2_01", "a2_10", "a2_11")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0)
      stop(paste0("'", nm, "' must be a single positive numeric value"))
  }

  # --- Threshold parameters ---
  if (prob == 'posterior') {
    for (nm in c("theta.TV1", "theta.MAV1", "theta.TV2", "theta.MAV2")) {
      val <- get(nm)
      if (is.null(val))
        stop(paste0("'", nm, "' must be non-NULL when prob = 'posterior'"))
      if (!is.numeric(val) || length(val) != 1L || is.na(val))
        stop(paste0("'", nm, "' must be a single numeric value"))
    }
    if (theta.TV1 <= theta.MAV1)
      stop("'theta.TV1' must be strictly greater than 'theta.MAV1'")
    if (theta.TV2 <= theta.MAV2)
      stop("'theta.TV2' must be strictly greater than 'theta.MAV2'")
  } else {
    for (nm in c("theta.NULL1", "theta.NULL2")) {
      val <- get(nm)
      if (is.null(val))
        stop(paste0("'", nm, "' must be non-NULL when prob = 'predictive'"))
      if (!is.numeric(val) || length(val) != 1L || is.na(val))
        stop(paste0("'", nm, "' must be a single numeric value"))
    }
    # For pbayespostpred2bin, predictive uses TV = MAV = NULL threshold
    theta.TV1  <- theta.NULL1;  theta.MAV1 <- theta.NULL1
    theta.TV2  <- theta.NULL2;  theta.MAV2 <- theta.NULL2
  }

  # --- Future sample sizes ---
  if (prob == 'predictive') {
    if (is.null(m1) || is.null(m2))
      stop("'m1' and 'm2' must be non-NULL when prob = 'predictive'")
    for (nm in c("m1", "m2")) {
      val <- get(nm)
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 1L)
        stop(paste0("'", nm, "' must be a single positive integer"))
    }
    m1 <- as.integer(m1); m2 <- as.integer(m2)
  }

  # --- Uncontrolled design: hypothetical control counts ---
  if (design == 'uncontrolled') {
    for (nm in c("z00", "z01", "z10", "z11")) {
      val <- get(nm)
      if (is.null(val))
        stop(paste0("'", nm, "' must be non-NULL when design = 'uncontrolled'"))
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 0L)
        stop(paste0("'", nm, "' must be a single non-negative integer"))
    }
  }

  # --- External design ---
  if (design == 'external') {
    has_ext1 <- !is.null(xe1_00) && !is.null(xe1_01) &&
      !is.null(xe1_10) && !is.null(xe1_11) && !is.null(ae1)
    has_ext2 <- !is.null(xe2_00) && !is.null(xe2_01) &&
      !is.null(xe2_10) && !is.null(xe2_11) && !is.null(ae2)
    if (!has_ext1 && !has_ext2)
      stop(paste0("For design = 'external', at least one complete set of ",
                  "external data must be provided"))
  }

  if (!is.logical(error_if_Miss) || length(error_if_Miss) != 1L ||
      is.na(error_if_Miss))
    stop("'error_if_Miss' must be a single logical value")

  if (!is.logical(Gray_inc_Miss) || length(Gray_inc_Miss) != 1L ||
      is.na(Gray_inc_Miss))
    stop("'Gray_inc_Miss' must be a single logical value")

  # ---------------------------------------------------------------------------
  # Section 2: Stage 1 -- Precompute decision probabilities for all count
  #            combinations using fully vectorised batch Gamma sampling.
  #
  # Key idea: Dirichlet posterior Dir(alpha + x) can be sampled via
  #   Gamma(alpha_k + x_k, 1) normalised by their row sum.
  # We split this into:
  #   Y_k ~ Gamma(alpha_k, 1)   [prior component, shared across all x]
  #   Z_k ~ Gamma(x_k,     1)   [data  component, specific to each x]
  # and generate all Y and Z draws in a single rgamma() call each,
  # then add them to form the unnormalised posterior samples.
  # This reduces function-call overhead from n_t * n_c calls to O(1) calls.
  # ---------------------------------------------------------------------------

  # Enumerate all possible count vectors for each arm
  counts_t <- allmultinom(n1)   # (n_t x 4) integer matrix
  counts_c <- allmultinom(n2)   # (n_c x 4) integer matrix
  n_t      <- nrow(counts_t)
  n_c      <- nrow(counts_c)

  # --- Build posterior Dirichlet base parameters (prior + external data) ---
  xe1_w <- if (!is.null(ae1) && design == 'external') ae1 else 0
  xe2_w <- if (!is.null(ae2) && design == 'external') ae2 else 0

  alpha1_base <- c(
    a1_00 + xe1_w * ifelse(!is.null(xe1_00), xe1_00, 0),
    a1_01 + xe1_w * ifelse(!is.null(xe1_01), xe1_01, 0),
    a1_10 + xe1_w * ifelse(!is.null(xe1_10), xe1_10, 0),
    a1_11 + xe1_w * ifelse(!is.null(xe1_11), xe1_11, 0)
  )

  if (design == 'uncontrolled') {
    # Hypothetical control: fixed Dir(a2 + z), no data component varies
    alpha2_fixed <- c(a2_00 + z00, a2_01 + z01, a2_10 + z10, a2_11 + z11)
  } else {
    alpha2_base <- c(
      a2_00 + xe2_w * ifelse(!is.null(xe2_00), xe2_00, 0),
      a2_01 + xe2_w * ifelse(!is.null(xe2_01), xe2_01, 0),
      a2_10 + xe2_w * ifelse(!is.null(xe2_10), xe2_10, 0),
      a2_11 + xe2_w * ifelse(!is.null(xe2_11), xe2_11, 0)
    )
  }

  # --- Batch Gamma sampling: prior component (nMC x 4) ---
  # Y1[u, k] ~ Gamma(alpha1_base[k], 1),  shared across all x_t
  # Y2[u, k] ~ Gamma(alpha2_base[k], 1),  shared across all x_c
  Y1 <- matrix(rgamma(nMC * 4L, shape = rep(alpha1_base, each = nMC)),
               nrow = nMC, ncol = 4L)

  if (design == 'uncontrolled') {
    # Control is fixed: sample once from Dir(alpha2_fixed)
    G2_fixed <- matrix(rgamma(nMC * 4L,
                              shape = rep(alpha2_fixed, each = nMC)),
                       nrow = nMC, ncol = 4L)
    S2_fixed  <- rowSums(G2_fixed)
    p2_fixed  <- G2_fixed / S2_fixed   # (nMC x 4) Dirichlet samples
  } else {
    Y2 <- matrix(rgamma(nMC * 4L, shape = rep(alpha2_base, each = nMC)),
                 nrow = nMC, ncol = 4L)
  }

  # --- For predictive: pre-generate uniform random draws for future data ---
  # (Sequential binomial approach; uniform variates reused across count combos
  #  is NOT valid because each combo has different p parameters.
  #  We therefore keep the sequential-binomial loop but operate on the
  #  nMC-length p vectors produced per combo from the batch Gamma above.)

  # --- Precompute Gamma data components for all count combinations ---
  # Z1[i, u, k] ~ Gamma(counts_t[i, k] + 1, 1) would be expensive;
  # instead we use the additive property:
  #   Gamma(alpha + x, 1) = Gamma(alpha, 1) + Gamma(x, 1)  [in distribution]
  # Z1_arr[i, u, k]: nMC draws of Gamma(counts_t[i,k], 1) for each i.
  # We generate all at once: shape vector of length n_t * nMC * 4.
  # To avoid Gamma(0, 1) (undefined), we clamp shapes to a small epsilon
  # and note that Gamma(0,1) = 0 a.s., so we handle zeros separately.

  # Vectorised Gamma draws for treatment arm data component:
  # shape matrix: (n_t x 4), replicated nMC times -> (n_t * nMC x 4)
  shapes_t <- counts_t + 0L   # keep as integer; Gamma(0,1) handled below

  # Draw Gamma(shape, 1) for all (i, k) pairs and all nMC replicates at once.
  # Layout: for each i, nMC consecutive rows -> total n_t * nMC rows, 4 cols.
  # shape for row (i-1)*nMC + u, col k  =  counts_t[i, k]
  shape_vec_t <- as.numeric(rep(t(shapes_t), each = nMC))  # length n_t*nMC*4
  # Replace 0 shapes with 1 temporarily; zero-count columns contribute 0
  zero_mask_t <- shape_vec_t == 0
  shape_vec_t[zero_mask_t] <- 1
  raw_t <- rgamma(n_t * nMC * 4L, shape = shape_vec_t, rate = 1)
  raw_t[zero_mask_t] <- 0   # Gamma(0,1) = 0 a.s.
  # Arrange into (n_t * nMC) x 4 matrix
  Z1_mat <- matrix(raw_t, nrow = n_t * nMC, ncol = 4L)

  if (design != 'uncontrolled') {
    shapes_c    <- counts_c + 0L
    shape_vec_c <- as.numeric(rep(t(shapes_c), each = nMC))
    zero_mask_c <- shape_vec_c == 0
    shape_vec_c[zero_mask_c] <- 1
    raw_c <- rgamma(n_c * nMC * 4L, shape = shape_vec_c, rate = 1)
    raw_c[zero_mask_c] <- 0
    Z2_mat <- matrix(raw_c, nrow = n_c * nMC, ncol = 4L)
  }

  # --- Helper: compute PrGo and PrNoGo from a set of (p1, p2) samples ---
  # p1, p2: (nMC x 4) matrices of Dirichlet samples
  .compute_PrGoNoGo <- function(p1, p2) {

    # Marginal response rates and treatment effects
    theta1 <- (p1[, 3L] + p1[, 4L]) - (p2[, 3L] + p2[, 4L])
    theta2 <- (p1[, 2L] + p1[, 4L]) - (p2[, 2L] + p2[, 4L])

    if (prob == 'posterior') {
      # Region index: row-major, Endpoint 1 slowest (3 x 3 grid -> 9 regions)
      r1 <- 3L - as.integer(theta1 > theta.MAV1) -
        as.integer(theta1 > theta.TV1)
      r2 <- 3L - as.integer(theta2 > theta.MAV2) -
        as.integer(theta2 > theta.TV2)
      region <- (r1 - 1L) * 3L + r2
      Pr_R   <- tabulate(region, nbins = 9L) / nMC

    } else {
      # Predictive: simulate future multinomial data via sequential binomials
      x1f <- matrix(0L, nrow = nMC, ncol = 4L)
      x2f <- matrix(0L, nrow = nMC, ncol = 4L)
      rem1 <- rep(m1, nMC);  rem2 <- rep(m2, nMC)
      u1   <- rep(0,  nMC);  u2   <- rep(0,  nMC)

      for (jj in seq_len(3L)) {
        d1 <- pmax(1 - u1, 0);  d2 <- pmax(1 - u2, 0)
        q1 <- pmin(pmax(ifelse(d1 > 0, p1[, jj] / d1, 0), 0), 1)
        q2 <- pmin(pmax(ifelse(d2 > 0, p2[, jj] / d2, 0), 0), 1)
        b1 <- rbinom(nMC, rem1, q1);  b2 <- rbinom(nMC, rem2, q2)
        x1f[, jj] <- b1;  x2f[, jj] <- b2
        rem1 <- rem1 - b1;  rem2 <- rem2 - b2
        u1   <- u1 + p1[, jj];  u2 <- u2 + p2[, jj]
      }
      x1f[, 4L] <- rem1;  x2f[, 4L] <- rem2

      theta1 <- (x1f[, 3L] + x1f[, 4L]) / m1 -
        (x2f[, 3L] + x2f[, 4L]) / m2
      theta2 <- (x1f[, 2L] + x1f[, 4L]) / m1 -
        (x2f[, 2L] + x2f[, 4L]) / m2

      r1     <- 2L - as.integer(theta1 > theta.TV1)
      r2     <- 2L - as.integer(theta2 > theta.TV2)
      region <- (r1 - 1L) * 2L + r2
      Pr_R   <- tabulate(region, nbins = 4L) / nMC
    }

    list(PrGo   = sum(Pr_R[GoRegions]),
         PrNoGo = sum(Pr_R[NoGoRegions]))
  }

  # --- Main Stage 1 loop: iterate over treatment count combinations ---
  # For the controlled/external case, we eliminate the inner j-loop by
  # computing region memberships for all n_c control combos simultaneously.
  #
  # For each treatment combo i:
  #   p1[u, k] = (Y1[u,k] + Z1[i,u,k]) / rowSum  -- (nMC x 4)
  #
  # theta1_t[u] = p1[u,3] + p1[u,4]  (treatment marginal for Endpoint 1)
  # theta2_t[u] = p1[u,2] + p1[u,4]  (treatment marginal for Endpoint 2)
  #
  # For each control combo j:
  #   theta1_c[u] = p2[u,3] + p2[u,4]
  #   theta2_c[u] = p2[u,2] + p2[u,4]
  #
  # region_ij = f(theta1_t[u] - theta1_c[u], theta2_t[u] - theta2_c[u])
  #
  # Key vectorisation: pre-compute normalised p2 for ALL j at once:
  #   p2_all: (n_c * nMC x 4)  -- all control combos stacked
  # Then for fixed i, broadcast theta_t (length nMC) against
  # theta_c (n_c blocks of nMC) to get all n_c results in one pass.

  PrGo_mat   <- matrix(0, nrow = n_t, ncol = n_c)
  PrNoGo_mat <- matrix(0, nrow = n_t, ncol = n_c)

  if (design == 'uncontrolled') {

    # --- Uncontrolled: single fixed control distribution ---
    for (i in seq_len(n_t)) {
      idx1 <- ((i - 1L) * nMC + 1L):(i * nMC)
      G1   <- Y1 + Z1_mat[idx1, , drop = FALSE]
      p1   <- G1 / rowSums(G1)

      res <- .compute_PrGoNoGo(p1, p2_fixed)
      # All columns share the same value (Stage 2 will use only column 1)
      PrGo_mat[i, ]   <- res$PrGo
      PrNoGo_mat[i, ] <- res$PrNoGo
    }

  } else {

    # --- Controlled / External: vectorise over all j for each i ---

    # Pre-normalise all control combos: p2_all is (n_c * nMC) x 4
    G2_all <- Y2[rep(seq_len(nMC), times = n_c), , drop = FALSE] +
      Z2_mat
    p2_all <- G2_all / rowSums(G2_all)

    # Pre-compute control marginals for all combos: each is length n_c * nMC
    # theta_c1[j-block, u] = p2_all[(j-1)*nMC+u, 3] + p2_all[(j-1)*nMC+u, 4]
    tc1_all <- p2_all[, 3L] + p2_all[, 4L]  # length n_c * nMC
    tc2_all <- p2_all[, 2L] + p2_all[, 4L]  # length n_c * nMC

    # Pre-compute control marginals reshaped as (nMC x n_c) matrices.
    # tc1_mat[u, j] = marginal for Endpoint 1 in control combo j, MC draw u
    # tc2_mat[u, j] = marginal for Endpoint 2 in control combo j, MC draw u
    # This avoids rep() inside the i-loop (case B optimisation).
    tc1_mat <- matrix(tc1_all, nrow = nMC, ncol = n_c)
    tc2_mat <- matrix(tc2_all, nrow = nMC, ncol = n_c)

    if (prob == 'predictive') {
      # Pre-compute future multinomial draws for ALL control combos at once.
      # Layout of p2_all: (n_c * nMC) x 4, blocks of nMC rows per combo j.
      # Sequential binomial for all n_c * nMC rows simultaneously (case C).
      x2f_all <- matrix(0L, nrow = n_c * nMC, ncol = 4L)
      rem2_all <- rep(m2, n_c * nMC)
      u2_all   <- rep(0,  n_c * nMC)
      for (jj in seq_len(3L)) {
        d2_all  <- pmax(1 - u2_all, 0)
        q2_all  <- pmin(pmax(
          ifelse(d2_all > 0, p2_all[, jj] / d2_all, 0), 0), 1)
        b2_all  <- rbinom(n_c * nMC, rem2_all, q2_all)
        x2f_all[, jj] <- b2_all
        rem2_all <- rem2_all - b2_all
        u2_all   <- u2_all + p2_all[, jj]
      }
      x2f_all[, 4L] <- rem2_all

      # Endpoint marginals for future control data: (nMC x n_c) matrices
      # fut_tc1_mat[u, j] = (x2f_all[(j-1)*nMC+u, 3] + x2f_all[..., 4]) / m2
      fut_tc1_mat <- matrix(
        (x2f_all[, 3L] + x2f_all[, 4L]) / m2, nrow = nMC, ncol = n_c)
      fut_tc2_mat <- matrix(
        (x2f_all[, 2L] + x2f_all[, 4L]) / m2, nrow = nMC, ncol = n_c)
    }

    for (i in seq_len(n_t)) {
      idx1 <- ((i - 1L) * nMC + 1L):(i * nMC)
      G1   <- Y1 + Z1_mat[idx1, , drop = FALSE]
      p1   <- G1 / rowSums(G1)

      # Treatment marginals: (nMC x 1) vectors
      tt1 <- p1[, 3L] + p1[, 4L]
      tt2 <- p1[, 2L] + p1[, 4L]

      if (prob == 'posterior') {

        # th1_mat[u, j] = tt1[u] - tc1_mat[u, j]  (nMC x n_c, no rep())
        th1_mat <- tt1 - tc1_mat   # R broadcasts (nMC) against (nMC x n_c)
        th2_mat <- tt2 - tc2_mat

        # Region classification: (nMC x n_c) integer matrices
        r1_mat <- 3L - (th1_mat > theta.MAV1) - (th1_mat > theta.TV1)
        r2_mat <- 3L - (th2_mat > theta.MAV2) - (th2_mat > theta.TV2)
        reg_mat <- (r1_mat - 1L) * 3L + r2_mat   # values in 1..9

        # Go/NoGo indicator matrices (nMC x n_c), then column means = PrGo[i,]
        # matrix() restores dimensions lost by %in% (which returns a plain vector)
        PrGo_mat[i, ]   <- colMeans(
          matrix(reg_mat %in% GoRegions,   nrow = nMC, ncol = n_c))
        PrNoGo_mat[i, ] <- colMeans(
          matrix(reg_mat %in% NoGoRegions, nrow = nMC, ncol = n_c))

      } else {

        # Predictive: simulate future treatment data for this combo i
        x1f <- matrix(0L, nrow = nMC, ncol = 4L)
        rem1 <- rep(m1, nMC)
        u1   <- rep(0,  nMC)
        for (jj in seq_len(3L)) {
          d1 <- pmax(1 - u1, 0)
          q1 <- pmin(pmax(ifelse(d1 > 0, p1[, jj] / d1, 0), 0), 1)
          b1 <- rbinom(nMC, rem1, q1)
          x1f[, jj] <- b1
          rem1 <- rem1 - b1
          u1   <- u1 + p1[, jj]
        }
        x1f[, 4L] <- rem1

        # Future treatment marginals: (nMC x 1) vectors
        fut_tt1 <- (x1f[, 3L] + x1f[, 4L]) / m1
        fut_tt2 <- (x1f[, 2L] + x1f[, 4L]) / m1

        # Difference matrices (nMC x n_c) using pre-computed control futures
        fth1_mat <- fut_tt1 - fut_tc1_mat
        fth2_mat <- fut_tt2 - fut_tc2_mat

        r1_mat  <- 2L - (fth1_mat > theta.TV1)
        r2_mat  <- 2L - (fth2_mat > theta.TV2)
        reg_mat <- (r1_mat - 1L) * 2L + r2_mat

        PrGo_mat[i, ]   <- colMeans(
          matrix(reg_mat %in% GoRegions,   nrow = nMC, ncol = n_c))
        PrNoGo_mat[i, ] <- colMeans(
          matrix(reg_mat %in% NoGoRegions, nrow = nMC, ncol = n_c))
      }
    }
  }

  # Decision indicator matrices
  ind_Go   <- (PrGo_mat   >= gamma1) & (PrNoGo_mat <  gamma2)
  ind_NoGo <- (PrGo_mat   <  gamma1) & (PrNoGo_mat >= gamma2)
  ind_Miss <- (PrGo_mat   >= gamma1) & (PrNoGo_mat >= gamma2)

  # ---------------------------------------------------------------------------
  # Section 3: Stage 2 -- Compute operating characteristics per scenario
  #            by weighting Stage 1 decisions with multinomial probabilities
  # ---------------------------------------------------------------------------

  # Result matrix: columns = Go, NoGo, Miss
  result_mat <- matrix(0, nrow = n_scen, ncol = 3L)

  for (s in seq_len(n_scen)) {

    # Convert (pi, rho) to four-cell probability vector for each arm
    p_t <- getjointbin(pi1 = pi_t1[s], pi2 = pi_t2[s], rho = rho_t[s])

    if (design == 'uncontrolled') {
      # Control arm weights do not depend on (pi_c, rho_c): use unit weight
      # The Stage 1 precomputation already fixed the hypothetical control
      # distribution via z*, so x_c iterates over a single row (all zeros)
      # conceptually -- but pbayespostpred2bin draws from Dir(a2 + z) regardless
      # of x_c.  We therefore fix x_c to zeros and iterate over x_t only.
      # Multinomial weight for control is degenerate (point mass at x_c=0000
      # is not meaningful); instead collapse the j loop to j=1 with weight 1.
      # Stage 1 already stored results for all (i,j) but only the i loop
      # matters; the j=1 column (all zeros) is the reference.
      # Correct approach: treat Stage 1 x_c index as dummy and sum over j=1.
      w_t <- apply(counts_t, 1L, function(x)
        dmultinom(x, size = n1, prob = p_t))

      result_mat[s, 1L] <- sum(ind_Go[,   1L] * w_t)
      result_mat[s, 2L] <- sum(ind_NoGo[, 1L] * w_t)
      result_mat[s, 3L] <- sum(ind_Miss[, 1L] * w_t)

    } else {
      # Controlled or external: both arms contribute multinomial weights
      p_c <- getjointbin(pi1 = pi_c1[s], pi2 = pi_c2[s], rho = rho_c[s])

      w_t <- apply(counts_t, 1L, function(x)
        dmultinom(x, size = n1, prob = p_t))
      w_c <- apply(counts_c, 1L, function(x)
        dmultinom(x, size = n2, prob = p_c))

      # Outer product of weights: w[i,j] = w_t[i] * w_c[j]
      w_mat <- outer(w_t, w_c)

      result_mat[s, 1L] <- sum(ind_Go   * w_mat)
      result_mat[s, 2L] <- sum(ind_NoGo * w_mat)
      result_mat[s, 3L] <- sum(ind_Miss * w_mat)
    }
  }

  # ---------------------------------------------------------------------------
  # Section 4: Assemble output
  # ---------------------------------------------------------------------------

  # Suppress floating-point rounding artefacts near zero before Miss check
  result_mat[result_mat < .Machine$double.eps ^ 0.25] <- 0

  # Check for positive Miss probability (after zeroing out numerical noise)
  if (error_if_Miss && any(result_mat[, 3L] > 0)) {
    stop("Positive Miss probability detected. Please re-consider the chosen thresholds.")
  }

  # Gray probability = complement of Go + NoGo (+ Miss if not included)
  if (Gray_inc_Miss) {
    GrayProb <- 1 - result_mat[, 1L] - result_mat[, 2L]
  } else {
    GrayProb <- 1 - rowSums(result_mat)
  }

  # Build the results data frame
  if (design == 'uncontrolled') {
    results <- data.frame(
      pi_t1 = pi_t1,
      pi_t2 = pi_t2,
      rho_t  = rho_t,
      Go    = result_mat[, 1L],
      Gray  = GrayProb,
      NoGo  = result_mat[, 2L]
    )
  } else {
    results <- data.frame(
      pi_t1 = pi_t1,
      pi_t2 = pi_t2,
      rho_t  = rho_t,
      pi_c1 = pi_c1,
      pi_c2 = pi_c2,
      rho_c  = rho_c,
      Go    = result_mat[, 1L],
      Gray  = GrayProb,
      NoGo  = result_mat[, 2L]
    )
  }

  if (!error_if_Miss && !Gray_inc_Miss) {
    results$Miss <- result_mat[, 3L]
  }

  # Attach metadata as attributes for use in print()
  attr(results, 'prob')          <- prob
  attr(results, 'design')        <- design
  attr(results, 'GoRegions')     <- GoRegions
  attr(results, 'NoGoRegions')   <- NoGoRegions
  attr(results, 'gamma1')        <- gamma1
  attr(results, 'gamma2')        <- gamma2
  attr(results, 'n1')            <- n1
  attr(results, 'n2')            <- n2
  attr(results, 'a1_00')         <- a1_00
  attr(results, 'a1_01')         <- a1_01
  attr(results, 'a1_10')         <- a1_10
  attr(results, 'a1_11')         <- a1_11
  attr(results, 'a2_00')         <- a2_00
  attr(results, 'a2_01')         <- a2_01
  attr(results, 'a2_10')         <- a2_10
  attr(results, 'a2_11')         <- a2_11
  attr(results, 'm1')            <- m1
  attr(results, 'm2')            <- m2
  attr(results, 'theta.TV1')     <- theta.TV1
  attr(results, 'theta.MAV1')    <- theta.MAV1
  attr(results, 'theta.TV2')     <- theta.TV2
  attr(results, 'theta.MAV2')    <- theta.MAV2
  attr(results, 'theta.NULL1')   <- theta.NULL1
  attr(results, 'theta.NULL2')   <- theta.NULL2
  attr(results, 'z00')           <- z00
  attr(results, 'z01')           <- z01
  attr(results, 'z10')           <- z10
  attr(results, 'z11')           <- z11
  attr(results, 'ae1')           <- ae1
  attr(results, 'ae2')           <- ae2
  attr(results, 'nMC')           <- nMC
  attr(results, 'error_if_Miss') <- error_if_Miss
  attr(results, 'Gray_inc_Miss') <- Gray_inc_Miss

  class(results) <- c('pbayesdecisionprob2bin', 'data.frame')

  return(results)
}

# ==============================================================================

#' Print Method for pbayesdecisionprob2bin Objects
#'
#' Displays a formatted summary of Go/NoGo/Gray decision probabilities for
#' two-binary-endpoint results returned by \code{\link{pbayesdecisionprob2bin}}.
#'
#' @param x An object of class \code{pbayesdecisionprob2bin}.
#' @param digits A positive integer specifying the number of decimal places
#'        for probability values.  Default is 4.
#' @param ... Further arguments passed to or from other methods (ignored).
#'
#' @return Invisibly returns \code{x}.
#'
#' @export
print.pbayesdecisionprob2bin <- function(x, digits = 4, ...) {

  # Helper: format a value as string (NULL -> "NULL")
  fmt <- function(v) if (is.null(v)) 'NULL' else as.character(v)

  # Extract metadata
  prob          <- attr(x, 'prob')
  design        <- attr(x, 'design')
  GoRegions     <- attr(x, 'GoRegions')
  NoGoRegions   <- attr(x, 'NoGoRegions')
  gamma1        <- attr(x, 'gamma1')
  gamma2        <- attr(x, 'gamma2')
  n1            <- attr(x, 'n1')
  n2            <- attr(x, 'n2')
  a1_00         <- attr(x, 'a1_00')
  a1_01         <- attr(x, 'a1_01')
  a1_10         <- attr(x, 'a1_10')
  a1_11         <- attr(x, 'a1_11')
  a2_00         <- attr(x, 'a2_00')
  a2_01         <- attr(x, 'a2_01')
  a2_10         <- attr(x, 'a2_10')
  a2_11         <- attr(x, 'a2_11')
  m1            <- attr(x, 'm1')
  m2            <- attr(x, 'm2')
  theta.TV1     <- attr(x, 'theta.TV1')
  theta.MAV1    <- attr(x, 'theta.MAV1')
  theta.TV2     <- attr(x, 'theta.TV2')
  theta.MAV2    <- attr(x, 'theta.MAV2')
  theta.NULL1   <- attr(x, 'theta.NULL1')
  theta.NULL2   <- attr(x, 'theta.NULL2')
  z00           <- attr(x, 'z00')
  z01           <- attr(x, 'z01')
  z10           <- attr(x, 'z10')
  z11           <- attr(x, 'z11')
  ae1           <- attr(x, 'ae1')
  ae2           <- attr(x, 'ae2')
  nMC           <- attr(x, 'nMC')
  error_if_Miss <- attr(x, 'error_if_Miss')
  Gray_inc_Miss <- attr(x, 'Gray_inc_Miss')

  # Build threshold string
  if (prob == 'posterior') {
    thresh_str <- sprintf(
      'TV1 = %s, MAV1 = %s, TV2 = %s, MAV2 = %s',
      fmt(theta.TV1), fmt(theta.MAV1), fmt(theta.TV2), fmt(theta.MAV2)
    )
  } else {
    thresh_str <- sprintf(
      'NULL1 = %s, NULL2 = %s',
      fmt(theta.NULL1), fmt(theta.NULL2)
    )
  }

  # Build Dirichlet prior string
  prior1_str <- sprintf('(%s, %s, %s, %s)',
                        fmt(a1_00), fmt(a1_01), fmt(a1_10), fmt(a1_11))
  prior2_str <- sprintf('(%s, %s, %s, %s)',
                        fmt(a2_00), fmt(a2_01), fmt(a2_10), fmt(a2_11))

  # Print header
  cat('Go/NoGo/Gray Decision Probabilities (Two Binary Endpoints)\n')
  cat(strrep('-', 65), '\n')
  cat(sprintf('  Probability type : %s\n', prob))
  cat(sprintf('  Design           : %s\n', design))
  cat(sprintf('  Threshold(s)     : %s\n', thresh_str))
  cat(sprintf('  Go  threshold    : gamma1 = %s\n', fmt(gamma1)))
  cat(sprintf('  NoGo threshold   : gamma2 = %s\n', fmt(gamma2)))
  cat(sprintf('  Go  regions      : {%s}\n',
              paste(GoRegions,   collapse = ', ')))
  cat(sprintf('  NoGo regions     : {%s}\n',
              paste(NoGoRegions, collapse = ', ')))
  cat(sprintf('  Sample size      : n1 = %s, n2 = %s\n', fmt(n1), fmt(n2)))
  cat(sprintf('  Prior Grp1 (Dir) : a = %s\n', prior1_str))
  cat(sprintf('  Prior Grp2 (Dir) : a = %s\n', prior2_str))

  if (design == 'uncontrolled') {
    cat(sprintf(
      '  Hyp. control     : z = (%s, %s, %s, %s) [z00, z01, z10, z11]\n',
      fmt(z00), fmt(z01), fmt(z10), fmt(z11)
    ))
  }
  if (prob == 'predictive') {
    cat(sprintf('  Future trial     : m1 = %s, m2 = %s\n', fmt(m1), fmt(m2)))
  }
  if (!is.null(ae1) || !is.null(ae2)) {
    cat(sprintf('  Power prior      : ae1 = %s, ae2 = %s\n',
                fmt(ae1), fmt(ae2)))
  }
  cat(sprintf('  MC samples (nMC) : %s\n', fmt(nMC)))
  cat(sprintf('  Miss handling    : error_if_Miss = %s, Gray_inc_Miss = %s\n',
              fmt(error_if_Miss), fmt(Gray_inc_Miss)))
  cat(strrep('-', 65), '\n')

  # Format probability columns only (not scenario columns)
  scenario_cols <- c('pi_t1', 'pi_t2', 'rho_t', 'pi_c1', 'pi_c2', 'rho_c')
  prob_cols     <- names(x)[!names(x) %in% scenario_cols]

  x_print <- x
  x_print[prob_cols] <- lapply(x[prob_cols], function(col) {
    formatC(col, digits = digits, format = 'f')
  })

  print.data.frame(x_print, row.names = FALSE, quote = FALSE)
  cat(strrep('-', 65), '\n')

  invisible(x)
}
