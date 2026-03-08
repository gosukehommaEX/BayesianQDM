#' Go/NoGo/Gray Decision Probabilities for a Clinical Trial with Two Binary
#' Endpoints
#'
#' Evaluates operating characteristics (Go, NoGo, Gray probabilities) for
#' clinical trials with two binary endpoints under the Bayesian framework.
#' The function supports controlled, uncontrolled, and external designs,
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
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or
#'        \code{'external'}.
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
#' @param n_t A positive integer giving the number of patients in the
#'        treatment group in the proof-of-concept (PoC) trial.
#' @param n_c A positive integer giving the number of patients in the
#'        control group in the PoC trial. For \code{design = 'uncontrolled'},
#'        this is the hypothetical control sample size (required for
#'        consistency with other designs).
#' @param a_t_00 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(0, 0)} response pattern in the treatment
#'        group.
#' @param a_t_01 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(0, 1)} response pattern in the treatment
#'        group.
#' @param a_t_10 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(1, 0)} response pattern in the treatment
#'        group.
#' @param a_t_11 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(1, 1)} response pattern in the treatment
#'        group.
#' @param a_c_00 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(0, 0)} response pattern in the control
#'        group. For \code{design = 'uncontrolled'}, serves as a
#'        hyperparameter of the hypothetical control distribution.
#' @param a_c_01 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(0, 1)} response pattern in the control
#'        group. For \code{design = 'uncontrolled'}, serves as a
#'        hyperparameter of the hypothetical control distribution.
#' @param a_c_10 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(1, 0)} response pattern in the control
#'        group. For \code{design = 'uncontrolled'}, serves as a
#'        hyperparameter of the hypothetical control distribution.
#' @param a_c_11 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the \code{(1, 1)} response pattern in the control
#'        group. For \code{design = 'uncontrolled'}, serves as a
#'        hyperparameter of the hypothetical control distribution.
#' @param m_t A positive integer giving the number of patients in the
#'        treatment group for the future trial. Required when
#'        \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param m_c A positive integer giving the number of patients in the
#'        control group for the future trial. Required when
#'        \code{prob = 'predictive'}; otherwise set to \code{NULL}.
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
#' @param z00 A non-negative integer giving the hypothetical control count
#'        for pattern \code{(0, 0)}. Required when
#'        \code{design = 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param z01 A non-negative integer giving the hypothetical control count
#'        for pattern \code{(0, 1)}. Required when
#'        \code{design = 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param z10 A non-negative integer giving the hypothetical control count
#'        for pattern \code{(1, 0)}. Required when
#'        \code{design = 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param z11 A non-negative integer giving the hypothetical control count
#'        for pattern \code{(1, 1)}. Required when
#'        \code{design = 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param xe_t_00 A non-negative integer giving the external treatment group
#'        count for pattern \code{(0, 0)}. Required when
#'        \code{design = 'external'} and external treatment data are used;
#'        otherwise \code{NULL}.
#' @param xe_t_01 A non-negative integer giving the external treatment group
#'        count for pattern \code{(0, 1)}. Required for external treatment
#'        data; otherwise \code{NULL}.
#' @param xe_t_10 A non-negative integer giving the external treatment group
#'        count for pattern \code{(1, 0)}. Required for external treatment
#'        data; otherwise \code{NULL}.
#' @param xe_t_11 A non-negative integer giving the external treatment group
#'        count for pattern \code{(1, 1)}. Required for external treatment
#'        data; otherwise \code{NULL}.
#' @param xe_c_00 A non-negative integer giving the external control group
#'        count for pattern \code{(0, 0)}. Required when
#'        \code{design = 'external'} and external control data are used;
#'        otherwise \code{NULL}.
#' @param xe_c_01 A non-negative integer giving the external control group
#'        count for pattern \code{(0, 1)}. Required for external control
#'        data; otherwise \code{NULL}.
#' @param xe_c_10 A non-negative integer giving the external control group
#'        count for pattern \code{(1, 0)}. Required for external control
#'        data; otherwise \code{NULL}.
#' @param xe_c_11 A non-negative integer giving the external control group
#'        count for pattern \code{(1, 1)}. Required for external control
#'        data; otherwise \code{NULL}.
#' @param alpha0e_t A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for the treatment group.  Required when external treatment
#'        data are used; otherwise \code{NULL}.
#' @param alpha0e_c A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for the control group.  Required when external control
#'        data are used; otherwise \code{NULL}.
#' @param nsim A positive integer giving the number of PoC count vectors
#'        sampled via \code{rmultinom} per arm per scenario when
#'        \code{CalcMethod = 'MC'} (outer loop).  The sampled vectors are
#'        deduplicated into \eqn{K_t} and \eqn{K_c} unique vectors
#'        (\eqn{K_t, K_c \ll} \code{nsim}); Dirichlet sampling is then
#'        performed only for these unique vectors using \code{nMC} draws each.
#'        Ignored when \code{CalcMethod = 'Exact'}.  Default is \code{10000}.
#' @param nMC A positive integer giving the number of Dirichlet draws used to
#'        evaluate the decision probability for each count combination in
#'        Stage 1.  Used by both \code{CalcMethod = 'Exact'} and
#'        \code{CalcMethod = 'MC'} (inner loop).  Default is \code{10000}.
#' @param CalcMethod A character string specifying the computation method.
#'        Must be \code{'Exact'} (default) or \code{'MC'}.
#'        \code{'Exact'} uses full enumeration of all possible multinomial
#'        count combinations (two-stage approach described in Details).
#'        \code{'MC'} samples \code{nsim} \eqn{(x_t, x_c)} pairs via
#'        \code{rmultinom}, deduplicates them into \eqn{K} unique pairs
#'        (\eqn{K \ll} \code{nsim}), calls
#'        \code{\link{pbayespostpred2bin}} for each unique pair to obtain
#'        Go/NoGo probabilities, and weights the decisions by the observed
#'        pair frequencies.  This reuses the same validated probability logic
#'        as the Exact method, avoiding any duplication of computation.
#'        The \code{'MC'} method trades some Monte Carlo variance for
#'        substantially reduced computation time when \code{n_t} and/or
#'        \code{n_c} are large.
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
#' The returned object has S3 class \code{pbayesdecisionprob2bin} with
#' an associated \code{print} method.
#'
#' @details
#' \strong{Two-stage computation.}
#'
#' Stage 1 (precomputation): All possible multinomial count combinations
#' \eqn{x_t} are enumerated using \code{\link{allmultinom}}.
#' For \code{design = 'controlled'} or \code{'external'}, all combinations
#' \eqn{x_c} are also enumerated and the Go/NoGo probability matrix
#' \code{PrGo[i, j]} has dimensions \eqn{n_t \times n_c}.
#' For \code{design = 'uncontrolled'}, the control distribution is fixed as
#' \eqn{\mathrm{Dir}(\alpha_{2,**} + z_{**})} and does not depend on any
#' observed control counts; the probability matrix therefore has dimensions
#' \eqn{n_t \times 1}, eliminating the \eqn{\binom{n_2+3}{3}}-row control
#' enumeration and the associated Gamma sampling.
#' This stage is independent of the true scenario parameters.
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
#'
#' # Example 1: Posterior probability, controlled design (rho > 0)
#' # Explores the impact of positive within-arm endpoint correlation.
#' pbayesdecisionprob2bin(
#'   prob        = 'posterior',
#'   design      = 'controlled',
#'   GoRegions   = 1L,
#'   NoGoRegions = 9L,
#'   gamma_go    = 0.52,
#'   gamma_nogo  = 0.52,
#'   pi_t1       = c(0.15, 0.25, 0.30, 0.40),
#'   pi_t2       = c(0.20, 0.30, 0.35, 0.45),
#'   rho_t       = rep(0.6, 4),
#'   pi_c1       = rep(0.15, 4),
#'   pi_c2       = rep(0.20, 4),
#'   rho_c       = rep(0.6, 4),
#'   n_t = 10, n_c = 10,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   m_t = NULL, m_c = NULL,
#'   theta_TV1   = 0.15, theta_MAV1 = 0.10,
#'   theta_TV2   = 0.15, theta_MAV2 = 0.10,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL,
#'   nMC = 100,
#'   error_if_Miss = FALSE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 2: Posterior probability, uncontrolled design
#' # Hypothetical control specified via pseudo-counts.
#' pbayesdecisionprob2bin(
#'   prob        = 'posterior',
#'   design      = 'uncontrolled',
#'   GoRegions   = 1L,
#'   NoGoRegions = 9L,
#'   gamma_go    = 0.27,
#'   gamma_nogo  = 0.36,
#'   pi_t1       = c(0.15, 0.30, 0.40),
#'   pi_t2       = c(0.20, 0.35, 0.45),
#'   rho_t       = rep(0.2, 3),
#'   pi_c1       = NULL,
#'   pi_c2       = NULL,
#'   rho_c       = NULL,
#'   n_t = 10, n_c = NULL,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   m_t = NULL, m_c = NULL,
#'   theta_TV1   = 0.15, theta_MAV1 = 0.10,
#'   theta_TV2   = 0.15, theta_MAV2 = 0.10,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   z00 = 2L, z01 = 1L, z10 = 2L, z11 = 1L,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL,
#'   nMC = 100,
#'   error_if_Miss = FALSE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 3: Posterior probability, external control design
#' # External control data incorporated via power prior (alpha0e_c = 0.5).
#' pbayesdecisionprob2bin(
#'   prob        = 'posterior',
#'   design      = 'external',
#'   GoRegions   = 1L,
#'   NoGoRegions = 9L,
#'   gamma_go    = 0.27,
#'   gamma_nogo  = 0.36,
#'   pi_t1       = c(0.15, 0.25, 0.30, 0.40),
#'   pi_t2       = c(0.20, 0.30, 0.35, 0.45),
#'   rho_t       = rep(0.0, 4),
#'   pi_c1       = rep(0.15, 4),
#'   pi_c2       = rep(0.20, 4),
#'   rho_c       = rep(0.0, 4),
#'   n_t = 10, n_c = 10,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   m_t = NULL, m_c = NULL,
#'   theta_TV1   = 0.15, theta_MAV1 = 0.10,
#'   theta_TV2   = 0.15, theta_MAV2 = 0.10,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = 4L,  xe_c_01 = 2L,  xe_c_10 = 3L,  xe_c_11 = 1L,
#'   alpha0e_t = NULL, alpha0e_c = 0.5,
#'   nMC = 100,
#'   error_if_Miss = FALSE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 4: Predictive probability, controlled design
#' # Future trial has m_t = m_c = 30 patients per arm.
#' pbayesdecisionprob2bin(
#'   prob        = 'predictive',
#'   design      = 'controlled',
#'   GoRegions   = 1L,
#'   NoGoRegions = 4L,
#'   gamma_go    = 0.60,
#'   gamma_nogo  = 0.80,
#'   pi_t1       = c(0.15, 0.25, 0.30, 0.40),
#'   pi_t2       = c(0.20, 0.30, 0.35, 0.45),
#'   rho_t       = rep(0.3, 4),
#'   pi_c1       = rep(0.15, 4),
#'   pi_c2       = rep(0.20, 4),
#'   rho_c       = rep(0.3, 4),
#'   n_t = 10, n_c = 10,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   m_t = 30L, m_c = 30L,
#'   theta_TV1   = NULL, theta_MAV1 = NULL,
#'   theta_TV2   = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 0.10, theta_NULL2 = 0.10,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL,
#'   nMC = 100,
#'   error_if_Miss = FALSE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 5: Predictive probability, uncontrolled design
#' # Hypothetical control specified via pseudo-counts; future trial m_t = m_c = 30.
#' pbayesdecisionprob2bin(
#'   prob        = 'predictive',
#'   design      = 'uncontrolled',
#'   GoRegions   = 1L,
#'   NoGoRegions = 4L,
#'   gamma_go    = 0.60,
#'   gamma_nogo  = 0.80,
#'   pi_t1       = c(0.15, 0.30, 0.40),
#'   pi_t2       = c(0.20, 0.35, 0.45),
#'   rho_t       = rep(0.7, 3),
#'   pi_c1       = rep(0.15, 3),
#'   pi_c2       = rep(0.20, 3),
#'   rho_c       = rep(0.7, 3),
#'   n_t = 10, n_c = 10,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   m_t = 30L, m_c = 30L,
#'   theta_TV1   = NULL, theta_MAV1 = NULL,
#'   theta_TV2   = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 0.10, theta_NULL2 = 0.10,
#'   z00 = 2L, z01 = 1L, z10 = 2L, z11 = 1L,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL,
#'   nMC = 100,
#'   error_if_Miss = FALSE, Gray_inc_Miss = FALSE
#' )
#'
#' # Example 6: Predictive probability, external treatment design
#' # External treatment data incorporated via power prior (alpha0e_t = 0.5).
#' pbayesdecisionprob2bin(
#'   prob        = 'predictive',
#'   design      = 'external',
#'   GoRegions   = 1L,
#'   NoGoRegions = 4L,
#'   gamma_go    = 0.60,
#'   gamma_nogo  = 0.80,
#'   pi_t1       = c(0.15, 0.25, 0.30, 0.40),
#'   pi_t2       = c(0.20, 0.30, 0.35, 0.45),
#'   rho_t       = rep(0.5, 4),
#'   pi_c1       = rep(0.15, 4),
#'   pi_c2       = rep(0.20, 4),
#'   rho_c       = rep(0.5, 4),
#'   n_t = 10, n_c = 10,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   m_t = 30L, m_c = 30L,
#'   theta_TV1   = NULL, theta_MAV1 = NULL,
#'   theta_TV2   = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 0.10, theta_NULL2 = 0.10,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe_t_00 = 3L, xe_t_01 = 2L, xe_t_10 = 3L, xe_t_11 = 2L,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   alpha0e_t = 0.5, alpha0e_c = NULL,
#'   nMC = 100,
#'   error_if_Miss = FALSE, Gray_inc_Miss = FALSE
#' )
#'
#' @importFrom stats dmultinom rmultinom rgamma rbinom
#' @export
pbayesdecisionprob2bin <- function(nsim        = 10000L,
                                   prob        = 'posterior',
                                   design      = 'controlled',
                                   GoRegions,
                                   NoGoRegions,
                                   gamma_go,
                                   gamma_nogo,
                                   pi_t1, pi_t2, rho_t,
                                   pi_c1 = NULL, pi_c2 = NULL, rho_c = NULL,
                                   n_t, n_c = NULL,
                                   a_t_00, a_t_01, a_t_10, a_t_11,
                                   a_c_00, a_c_01, a_c_10, a_c_11,
                                   m_t         = NULL,
                                   m_c         = NULL,
                                   theta_TV1   = NULL, theta_MAV1  = NULL,
                                   theta_TV2   = NULL, theta_MAV2  = NULL,
                                   theta_NULL1 = NULL, theta_NULL2 = NULL,
                                   z00 = NULL, z01 = NULL,
                                   z10 = NULL, z11 = NULL,
                                   xe_t_00 = NULL, xe_t_01 = NULL,
                                   xe_t_10 = NULL, xe_t_11 = NULL,
                                   xe_c_00 = NULL, xe_c_01 = NULL,
                                   xe_c_10 = NULL, xe_c_11 = NULL,
                                   alpha0e_t     = NULL,
                                   alpha0e_c     = NULL,
                                   nMC           = 10000L,
                                   CalcMethod    = 'Exact',
                                   error_if_Miss = TRUE,
                                   Gray_inc_Miss = FALSE) {

  # ---------------------------------------------------------------------------
  # Section 1: Input validation
  # ---------------------------------------------------------------------------

  if (!is.character(CalcMethod) || length(CalcMethod) != 1L ||
      !CalcMethod %in% c('Exact', 'MC'))
    stop("'CalcMethod' must be either 'Exact' or 'MC'")

  if (CalcMethod == 'MC') {
    if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) ||
        nsim != floor(nsim) || nsim < 1L)
      stop("'nsim' must be a single positive integer")
    nsim <- as.integer(nsim)
  }

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

  if (!is.numeric(gamma_go) || length(gamma_go) != 1L || is.na(gamma_go) ||
      gamma_go <= 0 || gamma_go >= 1) {
    stop("'gamma_go' must be a single numeric value in (0, 1)")
  }

  if (!is.numeric(gamma_nogo) || length(gamma_nogo) != 1L || is.na(gamma_nogo) ||
      gamma_nogo <= 0 || gamma_nogo >= 1) {
    stop("'gamma_nogo' must be a single numeric value in (0, 1)")
  }

  # No ordering constraint on gamma_go and gamma_nogo is imposed; their
  # combination determines the frequency of Miss outcomes.

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

  if (design != 'uncontrolled') {
    pi_c1 <- as.numeric(pi_c1)
    pi_c2 <- as.numeric(pi_c2)
    rho_c  <- as.numeric(rho_c)
    if (length(pi_c1) != n_scen || any(is.na(pi_c1)) ||
        any(pi_c1 <= 0) || any(pi_c1 >= 1))
      stop("'pi_c1' must have the same length as 'pi_t1' with all elements in (0, 1)")
    if (length(pi_c2) != n_scen || any(is.na(pi_c2)) ||
        any(pi_c2 <= 0) || any(pi_c2 >= 1))
      stop("'pi_c2' must have the same length as 'pi_t1' with all elements in (0, 1)")
    if (length(rho_c) != n_scen || any(is.na(rho_c)))
      stop("'rho_c' must have the same length as 'pi_t1'")
  }

  # --- Sample sizes ---
  if (!is.numeric(n_t) || length(n_t) != 1L || is.na(n_t) ||
      n_t != floor(n_t) || n_t < 1L)
    stop("'n_t' must be a single positive integer")
  n_t <- as.integer(n_t)

  if (design != 'uncontrolled') {
    if (!is.numeric(n_c) || length(n_c) != 1L || is.na(n_c) ||
        n_c != floor(n_c) || n_c < 1L)
      stop("'n_c' must be a single positive integer")
    n_c <- as.integer(n_c)
  }

  # --- Dirichlet prior parameters ---
  for (nm in c("a_t_00", "a_t_01", "a_t_10", "a_t_11",
               "a_c_00", "a_c_01", "a_c_10", "a_c_11")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0)
      stop(paste0("'", nm, "' must be a single positive numeric value"))
  }

  # --- Threshold parameters ---
  if (prob == 'posterior') {
    for (nm in c("theta_TV1", "theta_MAV1", "theta_TV2", "theta_MAV2")) {
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
    for (nm in c("theta_NULL1", "theta_NULL2")) {
      val <- get(nm)
      if (is.null(val))
        stop(paste0("'", nm, "' must be non-NULL when prob = 'predictive'"))
      if (!is.numeric(val) || length(val) != 1L || is.na(val))
        stop(paste0("'", nm, "' must be a single numeric value"))
    }
    # For pbayespostpred2bin, predictive uses TV = MAV = NULL threshold
    theta_TV1  <- theta_NULL1;  theta_MAV1 <- theta_NULL1
    theta_TV2  <- theta_NULL2;  theta_MAV2 <- theta_NULL2
  }

  # --- Future sample sizes ---
  if (prob == 'predictive') {
    if (is.null(m_t) || is.null(m_c))
      stop("'m_t' and 'm_c' must be non-NULL when prob = 'predictive'")
    for (nm in c("m_t", "m_c")) {
      val <- get(nm)
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 1L)
        stop(paste0("'", nm, "' must be a single positive integer"))
    }
    m_t <- as.integer(m_t); m_c <- as.integer(m_c)
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
    has_ext_t <- !is.null(xe_t_00) && !is.null(xe_t_01) &&
      !is.null(xe_t_10) && !is.null(xe_t_11) && !is.null(alpha0e_t)
    has_ext_c <- !is.null(xe_c_00) && !is.null(xe_c_01) &&
      !is.null(xe_c_10) && !is.null(xe_c_11) && !is.null(alpha0e_c)
    if (!has_ext_t && !has_ext_c)
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
  # Section 2 (MC): Unique-pair weighted Monte Carlo per scenario
  #
  # For each scenario s:
  #   1. Sample nsim count vectors via rmultinom for each arm.
  #   2. Deduplicate (x_t, x_c) pairs jointly; extract K unique pairs and
  #      their observed frequencies w (where K << nsim^2).
  #   3. Call pbayespostpred2bin for each unique pair to obtain PrGo/PrNoGo,
  #      reusing the same validated logic as the Exact method.
  #   4. Accumulate Go/NoGo/Miss probabilities weighted by w.
  # ---------------------------------------------------------------------------

  if (CalcMethod == 'MC') {

    # ---------------------------------------------------------------------------
    # MC method: sample nsim (x_t, x_c) pairs, deduplicate as joint pairs,
    # call pbayespostpred2bin for each unique pair, then weight by frequency.
    #
    # Using pbayespostpred2bin directly guarantees the same computation as the
    # Exact method's Stage 1, and avoids duplicating the posterior/predictive
    # probability logic.
    # ---------------------------------------------------------------------------

    result_mat <- matrix(0, nrow = n_scen, ncol = 3L)

    for (s in seq_len(n_scen)) {

      # Sample nsim treatment count vectors
      p_t    <- getjointbin(pi1 = pi_t1[s], pi2 = pi_t2[s], rho = rho_t[s])
      xt_raw <- t(rmultinom(nsim, size = n_t, prob = p_t))   # nsim x 4

      if (design == 'uncontrolled') {

        # Deduplicate treatment counts only
        xt_keys   <- apply(xt_raw, 1L, paste, collapse = "_")
        uniq_keys <- unique(xt_keys)
        w_pairs   <- tabulate(match(xt_keys, uniq_keys)) / nsim
        K         <- length(uniq_keys)
        xt_uniq   <- xt_raw[match(uniq_keys, xt_keys), , drop = FALSE]

        PrGo_vec   <- numeric(K)
        PrNoGo_vec <- numeric(K)

        for (i in seq_len(K)) {
          Pr_R <- pbayespostpred2bin(
            prob        = prob,        design      = design,
            theta_TV1   = theta_TV1,   theta_MAV1  = theta_MAV1,
            theta_TV2   = theta_TV2,   theta_MAV2  = theta_MAV2,
            theta_NULL1 = theta_NULL1, theta_NULL2 = theta_NULL2,
            x_t_00 = xt_uniq[i, 1L], x_t_01 = xt_uniq[i, 2L],
            x_t_10 = xt_uniq[i, 3L], x_t_11 = xt_uniq[i, 4L],
            x_c_00 = NULL, x_c_01 = NULL, x_c_10 = NULL, x_c_11 = NULL,
            a_t_00 = a_t_00, a_t_01 = a_t_01, a_t_10 = a_t_10, a_t_11 = a_t_11,
            a_c_00 = a_c_00, a_c_01 = a_c_01, a_c_10 = a_c_10, a_c_11 = a_c_11,
            m_t = m_t, m_c = m_c,
            z00 = z00, z01 = z01, z10 = z10, z11 = z11,
            xe_t_00 = xe_t_00, xe_t_01 = xe_t_01,
            xe_t_10 = xe_t_10, xe_t_11 = xe_t_11,
            xe_c_00 = xe_c_00, xe_c_01 = xe_c_01,
            xe_c_10 = xe_c_10, xe_c_11 = xe_c_11,
            alpha0e_t = alpha0e_t, alpha0e_c = alpha0e_c,
            nMC = nMC
          )
          PrGo_vec[i]   <- sum(Pr_R[GoRegions])
          PrNoGo_vec[i] <- sum(Pr_R[NoGoRegions])
        }

      } else {

        # Sample nsim control count vectors
        p_c    <- getjointbin(pi1 = pi_c1[s], pi2 = pi_c2[s], rho = rho_c[s])
        xc_raw <- t(rmultinom(nsim, size = n_c, prob = p_c))   # nsim x 4

        # Deduplicate (x_t, x_c) pairs jointly
        pair_keys <- paste(
          apply(xt_raw, 1L, paste, collapse = "_"),
          apply(xc_raw, 1L, paste, collapse = "_"),
          sep = "|"
        )
        uniq_keys <- unique(pair_keys)
        w_pairs   <- tabulate(match(pair_keys, uniq_keys)) / nsim
        K         <- length(uniq_keys)

        uniq_idx <- match(uniq_keys, pair_keys)
        xt_uniq  <- xt_raw[uniq_idx, , drop = FALSE]
        xc_uniq  <- xc_raw[uniq_idx, , drop = FALSE]

        PrGo_vec   <- numeric(K)
        PrNoGo_vec <- numeric(K)

        for (i in seq_len(K)) {
          Pr_R <- pbayespostpred2bin(
            prob        = prob,        design      = design,
            theta_TV1   = theta_TV1,   theta_MAV1  = theta_MAV1,
            theta_TV2   = theta_TV2,   theta_MAV2  = theta_MAV2,
            theta_NULL1 = theta_NULL1, theta_NULL2 = theta_NULL2,
            x_t_00 = xt_uniq[i, 1L], x_t_01 = xt_uniq[i, 2L],
            x_t_10 = xt_uniq[i, 3L], x_t_11 = xt_uniq[i, 4L],
            x_c_00 = xc_uniq[i, 1L], x_c_01 = xc_uniq[i, 2L],
            x_c_10 = xc_uniq[i, 3L], x_c_11 = xc_uniq[i, 4L],
            a_t_00 = a_t_00, a_t_01 = a_t_01, a_t_10 = a_t_10, a_t_11 = a_t_11,
            a_c_00 = a_c_00, a_c_01 = a_c_01, a_c_10 = a_c_10, a_c_11 = a_c_11,
            m_t = m_t, m_c = m_c,
            z00 = z00, z01 = z01, z10 = z10, z11 = z11,
            xe_t_00 = xe_t_00, xe_t_01 = xe_t_01,
            xe_t_10 = xe_t_10, xe_t_11 = xe_t_11,
            xe_c_00 = xe_c_00, xe_c_01 = xe_c_01,
            xe_c_10 = xe_c_10, xe_c_11 = xe_c_11,
            alpha0e_t = alpha0e_t, alpha0e_c = alpha0e_c,
            nMC = nMC
          )
          PrGo_vec[i]   <- sum(Pr_R[GoRegions])
          PrNoGo_vec[i] <- sum(Pr_R[NoGoRegions])
        }
      }

      # --- Decision indicators (Go, NoGo, Miss are mutually exclusive; Gray is the complement) ---
      # Accumulate Go/NoGo/Miss using frequency weights
      ind_Go   <- (PrGo_vec   >= gamma_go) & (PrNoGo_vec <  gamma_nogo)
      ind_NoGo <- (PrGo_vec   <  gamma_go) & (PrNoGo_vec >= gamma_nogo)
      ind_Miss <- (PrGo_vec   >= gamma_go) & (PrNoGo_vec >= gamma_nogo)

      result_mat[s, 1L] <- sum(ind_Go   * w_pairs)
      result_mat[s, 2L] <- sum(ind_NoGo * w_pairs)
      result_mat[s, 3L] <- sum(ind_Miss * w_pairs)
    }

  } else {

    # ---------------------------------------------------------------------------
    # Section 2 (Exact): Stage 1 -- Precompute decision probabilities for all
    #            count combinations using fully vectorised batch Gamma sampling.
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
    counts_t <- allmultinom(n_t)   # (n_t x 4) integer matrix
    n_row_t      <- nrow(counts_t)

    # For uncontrolled design, the control distribution is fixed (Dir(a2 + z))
    # and does not vary with observed control counts.  The Stage 1 probability
    # matrix therefore needs only a single column (n_c = 1), avoiding the
    # generation of the full allmultinom(n_c) table (C(n_c+3,3) rows).
    if (design == 'uncontrolled') {
      n_row_c <- 1L
    } else {
      counts_c <- allmultinom(n_c)   # (n_c x 4) integer matrix
      n_row_c      <- nrow(counts_c)
    }

    # --- Build posterior Dirichlet base parameters (prior + external data) ---
    xe_t_w <- if (!is.null(alpha0e_t) && design == 'external') alpha0e_t else 0
    xe_c_w <- if (!is.null(alpha0e_c) && design == 'external') alpha0e_c else 0

    alpha_t_base <- c(
      a_t_00 + xe_t_w * ifelse(!is.null(xe_t_00), xe_t_00, 0),
      a_t_01 + xe_t_w * ifelse(!is.null(xe_t_01), xe_t_01, 0),
      a_t_10 + xe_t_w * ifelse(!is.null(xe_t_10), xe_t_10, 0),
      a_t_11 + xe_t_w * ifelse(!is.null(xe_t_11), xe_t_11, 0)
    )

    if (design == 'uncontrolled') {
      # Hypothetical control: fixed Dir(a2 + z), no data component varies
      alpha_c_fixed <- c(a_c_00 + z00, a_c_01 + z01, a_c_10 + z10, a_c_11 + z11)
    } else {
      alpha_c_base <- c(
        a_c_00 + xe_c_w * ifelse(!is.null(xe_c_00), xe_c_00, 0),
        a_c_01 + xe_c_w * ifelse(!is.null(xe_c_01), xe_c_01, 0),
        a_c_10 + xe_c_w * ifelse(!is.null(xe_c_10), xe_c_10, 0),
        a_c_11 + xe_c_w * ifelse(!is.null(xe_c_11), xe_c_11, 0)
      )
    }

    # --- Batch Gamma sampling: prior component (nMC x 4) ---
    # Y_t[u, k] ~ Gamma(alpha_t_base[k], 1),  shared across all x_t
    # Y_c[u, k] ~ Gamma(alpha_c_base[k], 1),  shared across all x_c
    Y_t <- matrix(rgamma(nMC * 4L, shape = rep(alpha_t_base, each = nMC)),
                  nrow = nMC, ncol = 4L)

    if (design == 'uncontrolled') {
      # Control is fixed: sample once from Dir(alpha_c_fixed)
      G_c_fixed <- matrix(rgamma(nMC * 4L,
                                 shape = rep(alpha_c_fixed, each = nMC)),
                          nrow = nMC, ncol = 4L)
      S_c_fixed  <- rowSums(G_c_fixed)
      p_c_fixed  <- G_c_fixed / S_c_fixed   # (nMC x 4) Dirichlet samples
    } else {
      Y_c <- matrix(rgamma(nMC * 4L, shape = rep(alpha_c_base, each = nMC)),
                    nrow = nMC, ncol = 4L)
    }

    # --- For predictive: pre-generate uniform random draws for future data ---
    # (Sequential binomial approach; uniform variates reused across count combos
    #  is NOT valid because each combo has different p parameters.
    #  We therefore keep the sequential-binomial loop but operate on the
    #  nMC-length p vectors produced per combo from the batch Gamma above.)

    # --- Precompute Gamma data components for all count combinations ---
    # Z_t[i, u, k] ~ Gamma(counts_t[i, k] + 1, 1) would be expensive;
    # instead we use the additive property:
    #   Gamma(alpha + x, 1) = Gamma(alpha, 1) + Gamma(x, 1)  [in distribution]
    # Z_t_arr[i, u, k]: nMC draws of Gamma(counts_t[i,k], 1) for each i.
    # We generate all at once: shape vector of length n_row_t * nMC * 4.
    # To avoid Gamma(0, 1) (undefined), we clamp shapes to a small epsilon
    # and note that Gamma(0,1) = 0 a.s., so we handle zeros separately.

    # Vectorised Gamma draws for treatment arm data component.
    # Replicate each row of counts_t nMC times: row block i holds counts_t[i, ]
    # for nMC consecutive rows, giving the correct (n_row_t * nMC x 4) shape matrix.
    shapes_t_rep <- counts_t[rep(seq_len(n_row_t), each = nMC), , drop = FALSE]
    zero_mask_t  <- shapes_t_rep == 0L
    shapes_t_rep[zero_mask_t] <- 1L
    raw_t <- rgamma(n_row_t * nMC * 4L, shape = shapes_t_rep, rate = 1)
    raw_t[zero_mask_t] <- 0   # Gamma(0,1) = 0 a.s.
    Z_t_mat <- matrix(raw_t, nrow = n_row_t * nMC, ncol = 4L)

    if (design != 'uncontrolled') {
      shapes_c_rep <- counts_c[rep(seq_len(n_row_c), each = nMC), , drop = FALSE]
      zero_mask_c  <- shapes_c_rep == 0L
      shapes_c_rep[zero_mask_c] <- 1L
      raw_c <- rgamma(n_row_c * nMC * 4L, shape = shapes_c_rep, rate = 1)
      raw_c[zero_mask_c] <- 0
      Z_c_mat <- matrix(raw_c, nrow = n_row_c * nMC, ncol = 4L)
    }

    # --- Helper: compute PrGo and PrNoGo from a set of (p_t, p_c) samples ---
    # p_t, p_c: (nMC x 4) matrices of Dirichlet samples
    .compute_PrGoNoGo <- function(p_t, p_c) {

      # Marginal response rates and treatment effects
      theta1 <- (p_t[, 3L] + p_t[, 4L]) - (p_c[, 3L] + p_c[, 4L])
      theta2 <- (p_t[, 2L] + p_t[, 4L]) - (p_c[, 2L] + p_c[, 4L])

      if (prob == 'posterior') {
        # Region index: row-major, Endpoint 1 slowest (3 x 3 grid -> 9 regions)
        r1 <- 3L - as.integer(theta1 > theta_MAV1) -
          as.integer(theta1 > theta_TV1)
        r2 <- 3L - as.integer(theta2 > theta_MAV2) -
          as.integer(theta2 > theta_TV2)
        region <- (r1 - 1L) * 3L + r2
        Pr_R   <- tabulate(region, nbins = 9L) / nMC

      } else {
        # Predictive: simulate future multinomial data via sequential binomials
        x1f <- matrix(0L, nrow = nMC, ncol = 4L)
        x2f <- matrix(0L, nrow = nMC, ncol = 4L)
        rem1 <- rep(m_t, nMC);  rem2 <- rep(m_c, nMC)
        u1   <- rep(0,  nMC);  u2   <- rep(0,  nMC)

        for (jj in seq_len(3L)) {
          d1 <- pmax(1 - u1, 0);  d2 <- pmax(1 - u2, 0)
          q1 <- pmin(pmax(ifelse(d1 > 0, p_t[, jj] / d1, 0), 0), 1)
          q2 <- pmin(pmax(ifelse(d2 > 0, p_c[, jj] / d2, 0), 0), 1)
          b1 <- rbinom(nMC, rem1, q1);  b2 <- rbinom(nMC, rem2, q2)
          x1f[, jj] <- b1;  x2f[, jj] <- b2
          rem1 <- rem1 - b1;  rem2 <- rem2 - b2
          u1   <- u1 + p_t[, jj];  u2 <- u2 + p_c[, jj]
        }
        x1f[, 4L] <- rem1;  x2f[, 4L] <- rem2

        theta1 <- (x1f[, 3L] + x1f[, 4L]) / m_t -
          (x2f[, 3L] + x2f[, 4L]) / m_c
        theta2 <- (x1f[, 2L] + x1f[, 4L]) / m_t -
          (x2f[, 2L] + x2f[, 4L]) / m_c

        r1     <- 2L - as.integer(theta1 > theta_TV1)
        r2     <- 2L - as.integer(theta2 > theta_TV2)
        region <- (r1 - 1L) * 2L + r2
        Pr_R   <- tabulate(region, nbins = 4L) / nMC
      }

      list(PrGo   = sum(Pr_R[GoRegions]),
           PrNoGo = sum(Pr_R[NoGoRegions]))
    }

    # --- Main Stage 1 loop: iterate over treatment count combinations ---
    # For the controlled/external case, we eliminate the inner j-loop by
    # computing region memberships for all n_row_c control combos simultaneously.
    #
    # For each treatment combo i:
    #   p_t[u, k] = (Y_t[u,k] + Z_t[i,u,k]) / rowSum  -- (nMC x 4)
    #
    # theta1_t[u] = p_t[u,3] + p_t[u,4]  (treatment marginal for Endpoint 1)
    # theta2_t[u] = p_t[u,2] + p_t[u,4]  (treatment marginal for Endpoint 2)
    #
    # For each control combo j:
    #   theta1_c[u] = p_c[u,3] + p_c[u,4]
    #   theta2_c[u] = p_c[u,2] + p_c[u,4]
    #
    # region_ij = f(theta1_t[u] - theta1_c[u], theta2_t[u] - theta2_c[u])
    #
    # Key vectorisation: pre-compute normalised p_c for ALL j at once:
    #   p_c_all: (n_row_c * nMC x 4)  -- all control combos stacked
    # Then for fixed i, broadcast theta_t (length nMC) against
    # theta_c (n_row_c blocks of nMC) to get all n_row_c results in one pass.

    PrGo_mat   <- matrix(0, nrow = n_row_t, ncol = n_row_c)
    PrNoGo_mat <- matrix(0, nrow = n_row_t, ncol = n_row_c)

    if (design == 'uncontrolled') {

      # --- Uncontrolled: single fixed control distribution ---
      for (i in seq_len(n_row_t)) {
        idx1 <- ((i - 1L) * nMC + 1L):(i * nMC)
        G1   <- Y_t + Z_t_mat[idx1, , drop = FALSE]
        p_t   <- G1 / rowSums(G1)

        res <- .compute_PrGoNoGo(p_t, p_c_fixed)
        # n_row_c = 1 for uncontrolled, so only column 1 exists
        PrGo_mat[i, 1L]   <- res$PrGo
        PrNoGo_mat[i, 1L] <- res$PrNoGo
      }

    } else {

      # --- Controlled / External: vectorise over all j for each i ---

      # Pre-normalise all control combos: p_c_all is (n_row_c * nMC) x 4
      G_c_all <- Y_c[rep(seq_len(nMC), times = n_row_c), , drop = FALSE] + Z_c_mat
      p_c_all <- G_c_all / rowSums(G_c_all)

      # Pre-compute control marginals for all combos: each is length n_row_c * nMC
      # theta_c1[j-block, u] = p_c_all[(j-1)*nMC+u, 3] + p_c_all[(j-1)*nMC+u, 4]
      tc1_all <- p_c_all[, 3L] + p_c_all[, 4L]  # length n_row_c * nMC
      tc2_all <- p_c_all[, 2L] + p_c_all[, 4L]  # length n_row_c * nMC

      # Pre-compute control marginals reshaped as (nMC x n_row_c) matrices.
      # tc1_mat[u, j] = marginal for Endpoint 1 in control combo j, MC draw u
      # tc2_mat[u, j] = marginal for Endpoint 2 in control combo j, MC draw u
      # This avoids rep() inside the i-loop (case B optimisation).
      tc1_mat <- matrix(tc1_all, nrow = nMC, ncol = n_row_c)
      tc2_mat <- matrix(tc2_all, nrow = nMC, ncol = n_row_c)

      if (prob == 'predictive') {
        # Pre-compute future multinomial draws for ALL control combos at once.
        # Layout of p_c_all: (n_row_c * nMC) x 4, blocks of nMC rows per combo j.
        # Sequential binomial for all n_row_c * nMC rows simultaneously (case C).
        x2f_all <- matrix(0L, nrow = n_row_c * nMC, ncol = 4L)
        rem2_all <- rep(m_c, n_row_c * nMC)
        u2_all   <- rep(0,  n_row_c * nMC)
        for (jj in seq_len(3L)) {
          d2_all  <- pmax(1 - u2_all, 0)
          q2_all  <- pmin(pmax(
            ifelse(d2_all > 0, p_c_all[, jj] / d2_all, 0), 0), 1)
          b2_all  <- rbinom(n_row_c * nMC, rem2_all, q2_all)
          x2f_all[, jj] <- b2_all
          rem2_all <- rem2_all - b2_all
          u2_all   <- u2_all + p_c_all[, jj]
        }
        x2f_all[, 4L] <- rem2_all

        # Endpoint marginals for future control data: (nMC x n_row_c) matrices
        # fut_tc1_mat[u, j] = (x2f_all[(j-1)*nMC+u, 3] + x2f_all[..., 4]) / m_c
        fut_tc1_mat <- matrix(
          (x2f_all[, 3L] + x2f_all[, 4L]) / m_c, nrow = nMC, ncol = n_row_c)
        fut_tc2_mat <- matrix(
          (x2f_all[, 2L] + x2f_all[, 4L]) / m_c, nrow = nMC, ncol = n_row_c)
      }

      for (i in seq_len(n_row_t)) {
        idx1 <- ((i - 1L) * nMC + 1L):(i * nMC)
        G1   <- Y_t + Z_t_mat[idx1, , drop = FALSE]
        p_t   <- G1 / rowSums(G1)

        # Treatment marginals: (nMC x 1) vectors
        tt1 <- p_t[, 3L] + p_t[, 4L]
        tt2 <- p_t[, 2L] + p_t[, 4L]

        if (prob == 'posterior') {

          # th1_mat[u, j] = tt1[u] - tc1_mat[u, j]  (nMC x n_row_c, no rep())
          th1_mat <- tt1 - tc1_mat   # R broadcasts (nMC) against (nMC x n_row_c)
          th2_mat <- tt2 - tc2_mat

          # Region classification: (nMC x n_row_c) integer matrices
          r1_mat <- 3L - (th1_mat > theta_MAV1) - (th1_mat > theta_TV1)
          r2_mat <- 3L - (th2_mat > theta_MAV2) - (th2_mat > theta_TV2)
          reg_mat <- (r1_mat - 1L) * 3L + r2_mat   # values in 1..9

          # Go/NoGo indicator matrices (nMC x n_row_c), then column means = PrGo[i,]
          # matrix() restores dimensions lost by %in% (which returns a plain vector)
          PrGo_mat[i, ]   <- colMeans(
            matrix(reg_mat %in% GoRegions,   nrow = nMC, ncol = n_row_c))
          PrNoGo_mat[i, ] <- colMeans(
            matrix(reg_mat %in% NoGoRegions, nrow = nMC, ncol = n_row_c))

        } else {

          # Predictive: simulate future treatment data for this combo i
          x1f <- matrix(0L, nrow = nMC, ncol = 4L)
          rem1 <- rep(m_t, nMC)
          u1   <- rep(0,  nMC)
          for (jj in seq_len(3L)) {
            d1 <- pmax(1 - u1, 0)
            q1 <- pmin(pmax(ifelse(d1 > 0, p_t[, jj] / d1, 0), 0), 1)
            b1 <- rbinom(nMC, rem1, q1)
            x1f[, jj] <- b1
            rem1 <- rem1 - b1
            u1   <- u1 + p_t[, jj]
          }
          x1f[, 4L] <- rem1

          # Future treatment marginals: (nMC x 1) vectors
          fut_tt1 <- (x1f[, 3L] + x1f[, 4L]) / m_t
          fut_tt2 <- (x1f[, 2L] + x1f[, 4L]) / m_t

          # Difference matrices (nMC x n_row_c) using pre-computed control futures
          fth1_mat <- fut_tt1 - fut_tc1_mat
          fth2_mat <- fut_tt2 - fut_tc2_mat

          r1_mat  <- 2L - (fth1_mat > theta_TV1)
          r2_mat  <- 2L - (fth2_mat > theta_TV2)
          reg_mat <- (r1_mat - 1L) * 2L + r2_mat

          PrGo_mat[i, ]   <- colMeans(
            matrix(reg_mat %in% GoRegions,   nrow = nMC, ncol = n_row_c))
          PrNoGo_mat[i, ] <- colMeans(
            matrix(reg_mat %in% NoGoRegions, nrow = nMC, ncol = n_row_c))
        }
      }
    }

    # --- Decision indicators (Go, NoGo, Miss are mutually exclusive; Gray is the complement) ---
    ind_Go   <- (PrGo_mat   >= gamma_go) & (PrNoGo_mat <  gamma_nogo)
    ind_NoGo <- (PrGo_mat   <  gamma_go) & (PrNoGo_mat >= gamma_nogo)
    ind_Miss <- (PrGo_mat   >= gamma_go) & (PrNoGo_mat >= gamma_nogo)

    # ---------------------------------------------------------------------------
    # Section 3 (Exact): Stage 2 -- Compute operating characteristics per
    #            scenario by weighting Stage 1 decisions with multinomial
    #            probabilities
    # ---------------------------------------------------------------------------

    # Result matrix: columns = Go, NoGo, Miss
    result_mat <- matrix(0, nrow = n_scen, ncol = 3L)

    for (s in seq_len(n_scen)) {

      # Convert (pi, rho) to four-cell probability vector for each arm
      p_t <- getjointbin(pi1 = pi_t1[s], pi2 = pi_t2[s], rho = rho_t[s])

      if (design == 'uncontrolled') {
        # Control distribution is fixed (Dir(a2 + z)) and independent of x_c.
        # Stage 1 was computed with n_row_c = 1, so ind_Go/ind_NoGo/ind_Miss have
        # exactly one column.  Only treatment counts contribute multinomial
        # weights; the control column index is always 1.
        w_t <- apply(counts_t, 1L, function(x)
          dmultinom(x, size = n_t, prob = p_t))

        result_mat[s, 1L] <- sum(ind_Go[,   1L] * w_t)
        result_mat[s, 2L] <- sum(ind_NoGo[, 1L] * w_t)
        result_mat[s, 3L] <- sum(ind_Miss[, 1L] * w_t)

      } else {
        # Controlled or external: both arms contribute multinomial weights
        p_c <- getjointbin(pi1 = pi_c1[s], pi2 = pi_c2[s], rho = rho_c[s])

        w_t <- apply(counts_t, 1L, function(x)
          dmultinom(x, size = n_t, prob = p_t))
        w_c <- apply(counts_c, 1L, function(x)
          dmultinom(x, size = n_c, prob = p_c))

        # Outer product of weights: w[i,j] = w_t[i] * w_c[j]
        w_mat <- outer(w_t, w_c)

        result_mat[s, 1L] <- sum(ind_Go   * w_mat)
        result_mat[s, 2L] <- sum(ind_NoGo * w_mat)
        result_mat[s, 3L] <- sum(ind_Miss * w_mat)
      }
    }

  } # end if (CalcMethod == 'MC') ... else ...

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
  attr(results, 'prob')           <- prob
  attr(results, 'design')         <- design
  attr(results, 'GoRegions')      <- GoRegions
  attr(results, 'NoGoRegions')    <- NoGoRegions
  attr(results, 'gamma_go')       <- gamma_go
  attr(results, 'gamma_nogo')     <- gamma_nogo
  attr(results, 'n_t')            <- n_t
  attr(results, 'n_c')            <- n_c
  attr(results, 'a_t_00')         <- a_t_00
  attr(results, 'a_t_01')         <- a_t_01
  attr(results, 'a_t_10')         <- a_t_10
  attr(results, 'a_t_11')         <- a_t_11
  attr(results, 'a_c_00')         <- a_c_00
  attr(results, 'a_c_01')         <- a_c_01
  attr(results, 'a_c_10')         <- a_c_10
  attr(results, 'a_c_11')         <- a_c_11
  attr(results, 'm_t')            <- m_t
  attr(results, 'm_c')            <- m_c
  attr(results, 'theta_TV1')      <- theta_TV1
  attr(results, 'theta_MAV1')     <- theta_MAV1
  attr(results, 'theta_TV2')      <- theta_TV2
  attr(results, 'theta_MAV2')     <- theta_MAV2
  attr(results, 'theta_NULL1')    <- theta_NULL1
  attr(results, 'theta_NULL2')    <- theta_NULL2
  attr(results, 'z00')            <- z00
  attr(results, 'z01')            <- z01
  attr(results, 'z10')            <- z10
  attr(results, 'z11')            <- z11
  attr(results, 'xe_t_00')        <- xe_t_00
  attr(results, 'xe_t_01')        <- xe_t_01
  attr(results, 'xe_t_10')        <- xe_t_10
  attr(results, 'xe_t_11')        <- xe_t_11
  attr(results, 'xe_c_00')        <- xe_c_00
  attr(results, 'xe_c_01')        <- xe_c_01
  attr(results, 'xe_c_10')        <- xe_c_10
  attr(results, 'xe_c_11')        <- xe_c_11
  attr(results, 'alpha0e_t')      <- alpha0e_t
  attr(results, 'alpha0e_c')      <- alpha0e_c
  attr(results, 'nMC')            <- nMC
  attr(results, 'nsim')           <- nsim
  attr(results, 'CalcMethod')     <- CalcMethod
  attr(results, 'error_if_Miss')  <- error_if_Miss
  attr(results, 'Gray_inc_Miss')  <- Gray_inc_Miss

  class(results) <- c('pbayesdecisionprob2bin', 'data.frame')

  return(results)
}

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
  fmt <- function(v) if (is.null(v)) "NULL" else as.character(v)

  # Extract metadata
  prob          <- attr(x, "prob")
  design        <- attr(x, "design")
  GoRegions     <- attr(x, "GoRegions")
  NoGoRegions   <- attr(x, "NoGoRegions")
  gamma_go      <- attr(x, "gamma_go")
  gamma_nogo    <- attr(x, "gamma_nogo")
  n_t           <- attr(x, "n_t")
  n_c           <- attr(x, "n_c")
  a_t_00        <- attr(x, "a_t_00")
  a_t_01        <- attr(x, "a_t_01")
  a_t_10        <- attr(x, "a_t_10")
  a_t_11        <- attr(x, "a_t_11")
  a_c_00        <- attr(x, "a_c_00")
  a_c_01        <- attr(x, "a_c_01")
  a_c_10        <- attr(x, "a_c_10")
  a_c_11        <- attr(x, "a_c_11")
  m_t           <- attr(x, "m_t")
  m_c           <- attr(x, "m_c")
  theta_TV1     <- attr(x, "theta_TV1")
  theta_MAV1    <- attr(x, "theta_MAV1")
  theta_TV2     <- attr(x, "theta_TV2")
  theta_MAV2    <- attr(x, "theta_MAV2")
  theta_NULL1   <- attr(x, "theta_NULL1")
  theta_NULL2   <- attr(x, "theta_NULL2")
  z00           <- attr(x, "z00")
  z01           <- attr(x, "z01")
  z10           <- attr(x, "z10")
  z11           <- attr(x, "z11")
  xe_t_00       <- attr(x, "xe_t_00")
  xe_t_01       <- attr(x, "xe_t_01")
  xe_t_10       <- attr(x, "xe_t_10")
  xe_t_11       <- attr(x, "xe_t_11")
  xe_c_00       <- attr(x, "xe_c_00")
  xe_c_01       <- attr(x, "xe_c_01")
  xe_c_10       <- attr(x, "xe_c_10")
  xe_c_11       <- attr(x, "xe_c_11")
  alpha0e_t     <- attr(x, "alpha0e_t")
  alpha0e_c     <- attr(x, "alpha0e_c")
  nMC           <- attr(x, "nMC")
  nsim          <- attr(x, "nsim")
  CalcMethod    <- attr(x, "CalcMethod")
  error_if_Miss <- attr(x, "error_if_Miss")
  Gray_inc_Miss <- attr(x, "Gray_inc_Miss")

  # Build info lines with fixed label width (lw) for consistent alignment
  lw  <- 17L   # label field width
  pad <- "  "  # left margin

  lines <- character(0)
  lines <- c(lines, sprintf("%s%-*s: %s",           pad, lw, "Probability type", prob))
  lines <- c(lines, sprintf("%s%-*s: %s",           pad, lw, "Design",           design))

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

  # Dirichlet prior: treatment and control on separate lines
  lines <- c(lines, sprintf(
    "%s%-*s: a_t = (%s, %s, %s, %s)  [a_00, a_01, a_10, a_11]",
    pad, lw, "Prior (treatment)",
    fmt(a_t_00), fmt(a_t_01), fmt(a_t_10), fmt(a_t_11)))
  lines <- c(lines, sprintf(
    "%s%-*s: a_c = (%s, %s, %s, %s)  [a_00, a_01, a_10, a_11]",
    pad, lw, "Prior (control)  ",
    fmt(a_c_00), fmt(a_c_01), fmt(a_c_10), fmt(a_c_11)))

  if (design == "uncontrolled") {
    lines <- c(lines, sprintf(
      "%s%-*s: z = (%s, %s, %s, %s)  [z00, z01, z10, z11]",
      pad, lw, "Hyp. control",
      fmt(z00), fmt(z01), fmt(z10), fmt(z11)))
  }
  if (prob == "predictive") {
    lines <- c(lines, sprintf("%s%-*s: m_t = %s, m_c = %s",
                              pad, lw, "Future trial", fmt(m_t), fmt(m_c)))
  }
  if (design == "external") {
    # External data: treatment and control counts on separate lines,
    # followed by power prior weights
    lines <- c(lines, sprintf(
      "%s%-*s: xe_t = (%s, %s, %s, %s)  [xe_00, xe_01, xe_10, xe_11]",
      pad, lw, "External (treat.)",
      fmt(xe_t_00), fmt(xe_t_01), fmt(xe_t_10), fmt(xe_t_11)))
    lines <- c(lines, sprintf(
      "%s%-*s: xe_c = (%s, %s, %s, %s)  [xe_00, xe_01, xe_10, xe_11]",
      pad, lw, "External (cont.) ",
      fmt(xe_c_00), fmt(xe_c_01), fmt(xe_c_10), fmt(xe_c_11)))
    lines <- c(lines, sprintf("%s%-*s: alpha0e_t = %s, alpha0e_c = %s",
                              pad, lw, "Power prior",
                              fmt(alpha0e_t), fmt(alpha0e_c)))
  }
  lines <- c(lines, sprintf("%s%-*s: %s",           pad, lw, "Method",           fmt(CalcMethod)))
  lines <- c(lines, sprintf("%s%-*s: nMC = %s",     pad, lw, "MC draws",         fmt(nMC)))
  if (CalcMethod == "MC") {
    lines <- c(lines, sprintf("%s%-*s: nsim = %s",  pad, lw, "Sim size",         fmt(nsim)))
  }
  lines <- c(lines, sprintf("%s%-*s: error_if_Miss = %s, Gray_inc_Miss = %s",
                            pad, lw, "Miss handling",
                            fmt(error_if_Miss), fmt(Gray_inc_Miss)))

  # Determine separator width dynamically from the longest line
  title     <- "Go/NoGo/Gray Decision Probabilities (Two Binary Endpoints)"
  sep_width <- max(nchar(title), max(nchar(lines)))
  sep       <- strrep("-", sep_width)

  # Print header block
  cat(title, "\n")
  cat(sep, "\n")
  for (ln in lines) cat(ln, "\n")
  cat(sep, "\n")

  # Format probability columns only (not scenario columns)
  scenario_cols <- c("pi_t1", "pi_t2", "rho_t", "pi_c1", "pi_c2", "rho_c")
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
