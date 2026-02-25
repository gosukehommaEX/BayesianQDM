#' Find Optimal Go/NoGo Thresholds for Two Binary Endpoints
#'
#' Computes the optimal Go threshold \eqn{\gamma_1} and NoGo threshold
#' \eqn{\gamma_2} for two binary endpoints by searching over a
#' two-dimensional grid of candidate value pairs.  The search follows a
#' two-stage precompute-then-sweep strategy: region probabilities
#' \eqn{\hat{g}_{Go,ij}} and \eqn{\hat{g}_{NoGo,ij}} are precomputed once
#' for every possible multinomial outcome combination \eqn{(x_t, x_c)}, and
#' operating characteristics under each \eqn{(\gamma_1, \gamma_2)} pair are
#' obtained by a fast weighted summation without additional Dirichlet draws.
#'
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or
#'        \code{'external'}.
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
#' @param pi_t1 A numeric scalar in \code{(0, 1)} giving the true treatment
#'        response probability for Endpoint 1 for the scenario used to
#'        search for the optimal thresholds.
#' @param pi_t2 A numeric scalar in \code{(0, 1)} giving the true treatment
#'        response probability for Endpoint 2.
#' @param rho_t A numeric scalar giving the within-arm correlation in the
#'        treatment arm.
#' @param pi_c1 A numeric scalar in \code{(0, 1)} giving the true control
#'        response probability for Endpoint 1.  Required for
#'        \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param pi_c2 A numeric scalar in \code{(0, 1)} giving the true control
#'        response probability for Endpoint 2.  Required for
#'        \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param rho_c A numeric scalar giving the within-arm correlation in the
#'        control arm.  Required for \code{design = 'controlled'} or
#'        \code{'external'}; set to \code{NULL} for
#'        \code{design = 'uncontrolled'}.
#' @param target_go A numeric scalar in \code{(0, 1)} giving the target
#'        value for \eqn{\Pr(\mathrm{Go})} used to select the optimal
#'        \eqn{\gamma_1}.  For each candidate \eqn{\gamma_1}, the
#'        worst-case \eqn{\Pr(\mathrm{Go})} over all \eqn{\gamma_2} in
#'        \code{gamma2_grid} is compared against \code{target_go} using
#'        the operator specified by \code{crit_go}.
#' @param target_nogo A numeric scalar in \code{(0, 1)} giving the target
#'        value for \eqn{\Pr(\mathrm{NoGo})} used to select the optimal
#'        \eqn{\gamma_2}.  For each candidate \eqn{\gamma_2}, the
#'        worst-case \eqn{\Pr(\mathrm{NoGo})} over all \eqn{\gamma_1} in
#'        \code{gamma1_grid} is compared against \code{target_nogo} using
#'        the operator specified by \code{crit_nogo}.
#' @param crit_go A character string specifying the comparison operator
#'        applied to \eqn{\Pr(\mathrm{Go})} when searching for
#'        \eqn{\gamma_1}.  Must be one of \code{"<"}, \code{"<="},
#'        \code{">"}, or \code{">="}.  Default is \code{"<"}.
#' @param crit_nogo A character string specifying the comparison operator
#'        applied to \eqn{\Pr(\mathrm{NoGo})} when searching for
#'        \eqn{\gamma_2}.  Must be one of \code{"<"}, \code{"<="},
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
#'        (control arm) in the PoC trial.
#' @param a1_00 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (0,0) response pattern in group 1.
#' @param a1_01 A positive numeric scalar; see \code{a1_00}.
#' @param a1_10 A positive numeric scalar; see \code{a1_00}.
#' @param a1_11 A positive numeric scalar; see \code{a1_00}.
#' @param a2_00 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (0,0) response pattern in group 2.
#' @param a2_01 A positive numeric scalar; see \code{a2_00}.
#' @param a2_10 A positive numeric scalar; see \code{a2_00}.
#' @param a2_11 A positive numeric scalar; see \code{a2_00}.
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
#' @param z00 A non-negative integer giving the hypothetical control
#'        count for pattern (0,0).  Required when
#'        \code{design = 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param z01 A non-negative integer; see \code{z00}.
#' @param z10 A non-negative integer; see \code{z00}.
#' @param z11 A non-negative integer; see \code{z00}.
#' @param xe1_00 A non-negative integer giving the external group 1
#'        count for pattern (0,0).  Required when
#'        \code{design = 'external'} and external treatment data are
#'        used; otherwise \code{NULL}.
#' @param xe1_01 A non-negative integer; see \code{xe1_00}.
#' @param xe1_10 A non-negative integer; see \code{xe1_00}.
#' @param xe1_11 A non-negative integer; see \code{xe1_00}.
#' @param xe2_00 A non-negative integer giving the external group 2
#'        count for pattern (0,0).  Required when
#'        \code{design = 'external'} and external control data are
#'        used; otherwise \code{NULL}.
#' @param xe2_01 A non-negative integer; see \code{xe2_00}.
#' @param xe2_10 A non-negative integer; see \code{xe2_00}.
#' @param xe2_11 A non-negative integer; see \code{xe2_00}.
#' @param ae1 A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for external group 1 data.  Required when external
#'        treatment data are used; otherwise \code{NULL}.
#' @param ae2 A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for external group 2 data.  Required when external
#'        control data are used; otherwise \code{NULL}.
#' @param nMC A positive integer giving the number of Dirichlet draws
#'        used to evaluate region probabilities for each count
#'        combination in Stage 1.  Default is \code{1000L}.
#' @param method A character string specifying the computation method
#'        passed to \code{\link{pbayespostpred2bin}} internally.
#'        Must be \code{'Exact'} (default) or \code{'MC'}.
#' @param nsim A positive integer giving the number of outer Monte Carlo
#'        samples used when \code{method = 'MC'}.  Default is
#'        \code{10000L}.
#' @param gamma1_grid A numeric vector of candidate Go threshold values
#'        in \code{(0, 1)} to search over.  Defaults to
#'        \code{seq(0.01, 0.99, by = 0.01)}.
#' @param gamma2_grid A numeric vector of candidate NoGo threshold
#'        values in \code{(0, 1)} to search over.  Defaults to
#'        \code{seq(0.01, 0.99, by = 0.01)}.
#'
#' @return A list of class \code{getgamma2bin} with the following
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
#' The function uses the same two-stage precompute-then-sweep strategy
#' as \code{\link{pbayesdecisionprob2bin}}.
#'
#' \strong{Stage 1 (precomputation)}: \code{\link{pbayespostpred2bin}}
#' is called for every possible multinomial outcome combination
#' \eqn{(x_t, x_c)} enumerated by \code{\link{allmultinom}}.
#' The resulting 9-element region probability vector is summed over
#' \code{GoRegions} and \code{NoGoRegions} to obtain
#' \eqn{\hat{g}_{Go,ij}} and \eqn{\hat{g}_{NoGo,ij}}, the
#' Go/NoGo probabilities for each combination, independent of the
#' decision thresholds \eqn{(\gamma_1, \gamma_2)}.
#'
#' \strong{Stage 2 (gamma sweep)}: For each pair
#' \eqn{(\gamma_1, \gamma_2)} in the two-dimensional grid, the
#' operating characteristics are computed as:
#' \deqn{\Pr(\mathrm{Go}) = \sum_{i,j} w_{ij}
#'   \mathbf{1}\!\left[\hat{g}_{Go,ij} \ge \gamma_1,\;
#'   \hat{g}_{NoGo,ij} < \gamma_2\right]}
#' \deqn{\Pr(\mathrm{NoGo}) = \sum_{i,j} w_{ij}
#'   \mathbf{1}\!\left[\hat{g}_{NoGo,ij} \ge \gamma_2,\;
#'   \hat{g}_{Go,ij} < \gamma_1\right]}
#' where \eqn{w_{ij} = P(X_t = x_{t,i}) \times P(X_c = x_{c,j})} are
#' the multinomial probability weights.  No additional Dirichlet draws
#' are required for this step.
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
#' # Example 1: Controlled design, posterior probability
#' # gamma1: smallest gamma1 s.t. max_{gamma2} Pr(Go)   < 0.05
#' # gamma2: largest  gamma2 s.t. max_{gamma1} Pr(NoGo) < 0.20
#' \dontrun{
#' getgamma2bin(
#'   prob = 'posterior', design = 'controlled',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   pi_t1 = 0.15, pi_t2 = 0.20, rho_t = 0.0,
#'   pi_c1 = 0.15, pi_c2 = 0.20, rho_c = 0.0,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n1 = 10L, n2 = 10L,
#'   a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
#'   a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
#'   theta.TV1   = 0.15, theta.MAV1 = 0.10,
#'   theta.TV2   = 0.15, theta.MAV2 = 0.10,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   m1 = NULL, m2 = NULL,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 100L, method = 'Exact',
#'   gamma1_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma2_grid = seq(0.01, 0.99, by = 0.01)
#' )
#' }
#'
#' # Example 2: Uncontrolled design, posterior probability
#' \dontrun{
#' getgamma2bin(
#'   prob = 'posterior', design = 'uncontrolled',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   pi_t1 = 0.15, pi_t2 = 0.20, rho_t = 0.0,
#'   pi_c1 = NULL, pi_c2 = NULL, rho_c = NULL,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n1 = 10L, n2 = 10L,
#'   a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
#'   a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
#'   theta.TV1   = 0.15, theta.MAV1 = 0.10,
#'   theta.TV2   = 0.15, theta.MAV2 = 0.10,
#'   theta.NULL1 = NULL, theta.NULL2 = NULL,
#'   m1 = NULL, m2 = NULL,
#'   z00 = 3L, z01 = 2L, z10 = 3L, z11 = 2L,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 100L, method = 'Exact',
#'   gamma1_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma2_grid = seq(0.01, 0.99, by = 0.01)
#' )
#' }
#'
#' # Example 3: Controlled design, predictive probability
#' \dontrun{
#' getgamma2bin(
#'   prob = 'predictive', design = 'controlled',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   pi_t1 = 0.15, pi_t2 = 0.20, rho_t = 0.0,
#'   pi_c1 = 0.15, pi_c2 = 0.20, rho_c = 0.0,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n1 = 10L, n2 = 10L,
#'   a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
#'   a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
#'   theta.TV1   = NULL, theta.MAV1 = NULL,
#'   theta.TV2   = NULL, theta.MAV2 = NULL,
#'   theta.NULL1 = 0.10, theta.NULL2 = 0.10,
#'   m1 = 30L, m2 = 30L,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 100L, method = 'Exact',
#'   gamma1_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma2_grid = seq(0.01, 0.99, by = 0.01)
#' )
#' }
#'
#' # Example 4: External design, predictive probability
#' \dontrun{
#' getgamma2bin(
#'   prob = 'predictive', design = 'external',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   pi_t1 = 0.15, pi_t2 = 0.20, rho_t = 0.0,
#'   pi_c1 = 0.15, pi_c2 = 0.20, rho_c = 0.0,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n1 = 10L, n2 = 10L,
#'   a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
#'   a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
#'   theta.TV1   = NULL, theta.MAV1 = NULL,
#'   theta.TV2   = NULL, theta.MAV2 = NULL,
#'   theta.NULL1 = 0.10, theta.NULL2 = 0.10,
#'   m1 = 30L, m2 = 30L,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe1_00 = 3L, xe1_01 = 2L, xe1_10 = 3L, xe1_11 = 2L,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = 0.5, ae2 = NULL,
#'   nMC = 100L, method = 'Exact',
#'   gamma1_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma2_grid = seq(0.01, 0.99, by = 0.01)
#' )
#' }
#'
#' @importFrom stats dmultinom
#' @export
getgamma2bin <- function(prob = 'posterior', design = 'controlled',
                         GoRegions, NoGoRegions,
                         pi_t1, pi_t2, rho_t,
                         pi_c1 = NULL, pi_c2 = NULL, rho_c = NULL,
                         target_go, target_nogo,
                         crit_go  = '<',        crit_nogo  = '<',
                         sel_go   = 'smallest', sel_nogo   = 'largest',
                         n1, n2,
                         a1_00, a1_01, a1_10, a1_11,
                         a2_00, a2_01, a2_10, a2_11,
                         theta.TV1   = NULL, theta.MAV1  = NULL,
                         theta.TV2   = NULL, theta.MAV2  = NULL,
                         theta.NULL1 = NULL, theta.NULL2 = NULL,
                         m1 = NULL, m2 = NULL,
                         z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
                         xe1_00 = NULL, xe1_01 = NULL,
                         xe1_10 = NULL, xe1_11 = NULL,
                         xe2_00 = NULL, xe2_01 = NULL,
                         xe2_10 = NULL, xe2_11 = NULL,
                         ae1 = NULL, ae2 = NULL,
                         nMC    = 1000L,
                         method = 'Exact',
                         nsim   = 10000L,
                         gamma1_grid = seq(0.01, 0.99, by = 0.01),
                         gamma2_grid = seq(0.01, 0.99, by = 0.01)) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  if (!is.character(prob) || length(prob) != 1L ||
      !prob %in% c('posterior', 'predictive')) {
    stop("'prob' must be either 'posterior' or 'predictive'")
  }

  if (!is.character(design) || length(design) != 1L ||
      !design %in% c('controlled', 'uncontrolled', 'external')) {
    stop("'design' must be 'controlled', 'uncontrolled', or 'external'")
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

  if (!is.numeric(pi_t1) || length(pi_t1) != 1L || is.na(pi_t1) ||
      pi_t1 <= 0 || pi_t1 >= 1) {
    stop("'pi_t1' must be a single numeric value in (0, 1)")
  }
  if (!is.numeric(pi_t2) || length(pi_t2) != 1L || is.na(pi_t2) ||
      pi_t2 <= 0 || pi_t2 >= 1) {
    stop("'pi_t2' must be a single numeric value in (0, 1)")
  }
  if (!is.numeric(rho_t) || length(rho_t) != 1L || is.na(rho_t)) {
    stop("'rho_t' must be a single numeric value")
  }

  if (design != 'uncontrolled') {
    if (is.null(pi_c1) || !is.numeric(pi_c1) || length(pi_c1) != 1L ||
        is.na(pi_c1) || pi_c1 <= 0 || pi_c1 >= 1) {
      stop("'pi_c1' must be a single numeric value in (0, 1) for controlled or external design")
    }
    if (is.null(pi_c2) || !is.numeric(pi_c2) || length(pi_c2) != 1L ||
        is.na(pi_c2) || pi_c2 <= 0 || pi_c2 >= 1) {
      stop("'pi_c2' must be a single numeric value in (0, 1) for controlled or external design")
    }
    if (is.null(rho_c) || !is.numeric(rho_c) || length(rho_c) != 1L ||
        is.na(rho_c)) {
      stop("'rho_c' must be a single numeric value for controlled or external design")
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

  for (nm in c("n1", "n2")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val != floor(val) || val < 1L) {
      stop(paste0("'", nm, "' must be a single positive integer"))
    }
  }

  for (nm in c("a1_00", "a1_01", "a1_10", "a1_11",
               "a2_00", "a2_01", "a2_10", "a2_11")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0) {
      stop(paste0("'", nm, "' must be a single positive numeric value"))
    }
  }

  if (prob == 'posterior') {
    for (nm in c("theta.TV1", "theta.MAV1", "theta.TV2", "theta.MAV2")) {
      val <- get(nm)
      if (is.null(val)) {
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

  if (design == 'uncontrolled') {
    for (nm in c("z00", "z01", "z10", "z11")) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L ||
          is.na(val) || val != floor(val) || val < 0L) {
        stop(paste0("'", nm, "' must be a non-negative integer when design = 'uncontrolled'"))
      }
    }
  }

  if (!is.character(method) || length(method) != 1L ||
      !method %in% c('Exact', 'MC')) {
    stop("'method' must be either 'Exact' or 'MC'")
  }
  if (!is.numeric(nMC) || length(nMC) != 1L || is.na(nMC) ||
      nMC != floor(nMC) || nMC < 1L) {
    stop("'nMC' must be a single positive integer")
  }
  if (method == 'MC') {
    if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) ||
        nsim != floor(nsim) || nsim < 1L) {
      stop("'nsim' must be a single positive integer when method = 'MC'")
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
  # Stage 1: Enumerate all count combinations, compute PrGo_hat and
  #          PrNoGo_hat per combination via pbayespostpred2bin
  # ---------------------------------------------------------------------------
  counts_t <- allmultinom(n1)   # (n_t x 4) integer matrix
  n_t      <- nrow(counts_t)

  if (design == 'uncontrolled') {
    # Control distribution is fixed (determined by z00..z11); only
    # treatment outcomes vary.
    counts_c <- matrix(0L, nrow = 1L, ncol = 4L)
    n_c      <- 1L
  } else {
    counts_c <- allmultinom(n2)
    n_c      <- nrow(counts_c)
  }

  # Multinomial probability weights for treatment outcomes
  p_t     <- getjointbin(pi1 = pi_t1, pi2 = pi_t2, rho = rho_t)
  log_w_t <- apply(counts_t, 1L, function(x) {
    dmultinom(x, size = n1, prob = p_t, log = TRUE)
  })

  if (design != 'uncontrolled') {
    p_c     <- getjointbin(pi1 = pi_c1, pi2 = pi_c2, rho = rho_c)
    log_w_c <- apply(counts_c, 1L, function(x) {
      dmultinom(x, size = n2, prob = p_c, log = TRUE)
    })
    # w_mat[i, j] = P(X_t = counts_t[i,]) * P(X_c = counts_c[j,])
    w_mat <- exp(outer(log_w_t, log_w_c, `+`))
  } else {
    w_mat <- matrix(exp(log_w_t), nrow = n_t, ncol = 1L)
  }

  # Preallocate region probability storage: PrGo_hat[i,j] and PrNoGo_hat[i,j]
  PrGo_hat   <- matrix(0, nrow = n_t, ncol = n_c)
  PrNoGo_hat <- matrix(0, nrow = n_t, ncol = n_c)

  for (i in seq_len(n_t)) {
    x1_i <- as.integer(counts_t[i, ])

    if (design == 'uncontrolled') {
      # Control counts not used: pass x2 as NULL
      Pr_R <- pbayespostpred2bin(
        prob    = prob, design  = design,
        theta.TV1   = theta.TV1,   theta.MAV1  = theta.MAV1,
        theta.TV2   = theta.TV2,   theta.MAV2  = theta.MAV2,
        theta.NULL1 = theta.NULL1, theta.NULL2 = theta.NULL2,
        x1_00 = x1_i[1L], x1_01 = x1_i[2L],
        x1_10 = x1_i[3L], x1_11 = x1_i[4L],
        x2_00 = NULL, x2_01 = NULL, x2_10 = NULL, x2_11 = NULL,
        a1_00 = a1_00, a1_01 = a1_01, a1_10 = a1_10, a1_11 = a1_11,
        a2_00 = a2_00, a2_01 = a2_01, a2_10 = a2_10, a2_11 = a2_11,
        m1 = m1, m2 = m2,
        z00 = z00, z01 = z01, z10 = z10, z11 = z11,
        xe1_00 = xe1_00, xe1_01 = xe1_01,
        xe1_10 = xe1_10, xe1_11 = xe1_11,
        xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
        ae1 = ae1, ae2 = NULL,
        nMC = nMC
      )
      PrGo_hat[i, 1L]   <- sum(Pr_R[GoRegions])
      PrNoGo_hat[i, 1L] <- sum(Pr_R[NoGoRegions])

    } else {

      for (j in seq_len(n_c)) {
        x2_j <- as.integer(counts_c[j, ])

        Pr_R <- pbayespostpred2bin(
          prob    = prob, design  = design,
          theta.TV1   = theta.TV1,   theta.MAV1  = theta.MAV1,
          theta.TV2   = theta.TV2,   theta.MAV2  = theta.MAV2,
          theta.NULL1 = theta.NULL1, theta.NULL2 = theta.NULL2,
          x1_00 = x1_i[1L], x1_01 = x1_i[2L],
          x1_10 = x1_i[3L], x1_11 = x1_i[4L],
          x2_00 = x2_j[1L], x2_01 = x2_j[2L],
          x2_10 = x2_j[3L], x2_11 = x2_j[4L],
          a1_00 = a1_00, a1_01 = a1_01, a1_10 = a1_10, a1_11 = a1_11,
          a2_00 = a2_00, a2_01 = a2_01, a2_10 = a2_10, a2_11 = a2_11,
          m1 = m1, m2 = m2,
          z00 = z00, z01 = z01, z10 = z10, z11 = z11,
          xe1_00 = xe1_00, xe1_01 = xe1_01,
          xe1_10 = xe1_10, xe1_11 = xe1_11,
          xe2_00 = xe2_00, xe2_01 = xe2_01,
          xe2_10 = xe2_10, xe2_11 = xe2_11,
          ae1 = ae1, ae2 = ae2,
          nMC = nMC
        )
        PrGo_hat[i, j]   <- sum(Pr_R[GoRegions])
        PrNoGo_hat[i, j] <- sum(Pr_R[NoGoRegions])
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Stage 2: Sweep (gamma1_grid x gamma2_grid)
  #
  # For each (g1, g2):
  #   Pr(Go)   = sum w[i,j] * I(PrGo_hat[i,j] >= g1  AND PrNoGo_hat[i,j] <  g2)
  #   Pr(NoGo) = sum w[i,j] * I(PrNoGo_hat[i,j] >= g2 AND PrGo_hat[i,j]  <  g1)
  # ---------------------------------------------------------------------------
  PrGo_grid   <- matrix(NA_real_, nrow = ng1, ncol = ng2)
  PrNoGo_grid <- matrix(NA_real_, nrow = ng1, ncol = ng2)

  for (k1 in seq_len(ng1)) {
    g1 <- gamma1_grid[k1]
    for (k2 in seq_len(ng2)) {
      g2 <- gamma2_grid[k2]

      go_mask   <- (PrGo_hat   >= g1) & (PrNoGo_hat <  g2)
      nogo_mask <- (PrNoGo_hat >= g2) & (PrGo_hat   <  g1)
      go_mask[is.na(go_mask)]     <- FALSE
      nogo_mask[is.na(nogo_mask)] <- FALSE

      PrGo_grid[k1, k2]   <- sum(w_mat[go_mask])
      PrNoGo_grid[k1, k2] <- sum(w_mat[nogo_mask])
    }
  }

  # ---------------------------------------------------------------------------
  # Stage 3: Select optimal (gamma1, gamma2)
  #
  # gamma1: for each gamma1 candidate, take max Pr(Go) over gamma2_grid;
  #         select the sel_go index satisfying crit_go vs target_go.
  # gamma2: for each gamma2 candidate, take max Pr(NoGo) over gamma1_grid;
  #         select the sel_nogo index satisfying crit_nogo vs target_nogo.
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

  # Worst-case Pr(Go) over gamma2_grid for each gamma1 candidate
  max_PrGo_per_g1   <- apply(PrGo_grid,   1L, max, na.rm = TRUE)
  opt1 <- .select_idx(.compare(max_PrGo_per_g1,   crit_go,   target_go),   sel_go)

  # Worst-case Pr(NoGo) over gamma1_grid for each gamma2 candidate
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

  class(result) <- 'getgamma2bin'
  return(result)
}
