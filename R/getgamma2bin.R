#' Find Optimal Go/NoGo Thresholds for Two Binary Endpoints
#'
#' Computes the optimal Go threshold \eqn{\gamma_go} and NoGo threshold
#' \eqn{\gamma_nogo} for two binary endpoints by searching over a
#' two-dimensional grid of candidate value pairs.  The search follows a
#' two-stage precompute-then-sweep strategy: region probabilities
#' \eqn{\hat{g}_{Go,ij}} and \eqn{\hat{g}_{NoGo,ij}} are precomputed once
#' for every possible multinomial outcome combination \eqn{(x_t, x_c)}, and
#' operating characteristics under each \eqn{(\gamma_go, \gamma_nogo)} pair are
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
#' @param rho_t A numeric scalar giving the within-group correlation in the
#'        treatment group.
#' @param pi_c1 A numeric scalar in \code{(0, 1)} giving the true control
#'        response probability for Endpoint 1.  Required for
#'        \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param pi_c2 A numeric scalar in \code{(0, 1)} giving the true control
#'        response probability for Endpoint 2.  Required for
#'        \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param rho_c A numeric scalar giving the within-group correlation in the
#'        control group.  Required for \code{design = 'controlled'} or
#'        \code{'external'}; set to \code{NULL} for
#'        \code{design = 'uncontrolled'}.
#' @param target_go A numeric scalar in \code{(0, 1)} giving the target
#'        value for \eqn{\Pr(\mathrm{Go})} used to select the optimal
#'        \eqn{\gamma_go}.  For each candidate \eqn{\gamma_go}, the
#'        worst-case \eqn{\Pr(\mathrm{Go})} over all \eqn{\gamma_nogo} in
#'        \code{gamma_nogo_grid} is compared against \code{target_go} using
#'        the operator specified by \code{crit_go}.
#' @param target_nogo A numeric scalar in \code{(0, 1)} giving the target
#'        value for \eqn{\Pr(\mathrm{NoGo})} used to select the optimal
#'        \eqn{\gamma_nogo}.  For each candidate \eqn{\gamma_nogo}, the
#'        worst-case \eqn{\Pr(\mathrm{NoGo})} over all \eqn{\gamma_go} in
#'        \code{gamma_go_grid} is compared against \code{target_nogo} using
#'        the operator specified by \code{crit_nogo}.
#' @param crit_go A character string specifying the comparison operator
#'        applied to \eqn{\Pr(\mathrm{Go})} when searching for
#'        \eqn{\gamma_go}.  Must be one of \code{"<"}, \code{"<="},
#'        \code{">"}, or \code{">="}.  Default is \code{"<"}.
#' @param crit_nogo A character string specifying the comparison operator
#'        applied to \eqn{\Pr(\mathrm{NoGo})} when searching for
#'        \eqn{\gamma_nogo}.  Must be one of \code{"<"}, \code{"<="},
#'        \code{">"}, or \code{">="}.  Default is \code{"<"}.
#' @param sel_go A character string specifying whether to select the
#'        \code{"smallest"} or \code{"largest"} value in
#'        \code{gamma_go_grid} satisfying the \code{crit_go} criterion.
#'        Default is \code{"smallest"}.
#' @param sel_nogo A character string specifying whether to select the
#'        \code{"smallest"} or \code{"largest"} value in
#'        \code{gamma_nogo_grid} satisfying the \code{crit_nogo} criterion.
#'        Default is \code{"largest"}.
#' @param n_t A positive integer giving the number of patients in the
#'        treatment group (treatment group) in the PoC trial.
#' @param n_c A positive integer giving the number of patients in the
#'        control group (control group) in the PoC trial.
#' @param a_t_00 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (0,0) response pattern in the treatment group.
#' @param a_t_01 A positive numeric scalar; see \code{a_t_00}.
#' @param a_t_10 A positive numeric scalar; see \code{a_t_00}.
#' @param a_t_11 A positive numeric scalar; see \code{a_t_00}.
#' @param a_c_00 A positive numeric scalar giving the Dirichlet prior
#'        parameter for the (0,0) response pattern in the control group.
#' @param a_c_01 A positive numeric scalar; see \code{a_c_00}.
#' @param a_c_10 A positive numeric scalar; see \code{a_c_00}.
#' @param a_c_11 A positive numeric scalar; see \code{a_c_00}.
#' @param theta_TV1 A numeric scalar giving the TV threshold for
#'        Endpoint 1.  Required when \code{prob = 'posterior'};
#'        otherwise set to \code{NULL}.
#' @param theta_MAV1 A numeric scalar giving the MAV threshold for
#'        Endpoint 1.  Required when \code{prob = 'posterior'};
#'        otherwise set to \code{NULL}.
#' @param theta_TV2 A numeric scalar giving the TV threshold for
#'        Endpoint 2.  Required when \code{prob = 'posterior'};
#'        otherwise set to \code{NULL}.
#' @param theta_MAV2 A numeric scalar giving the MAV threshold for
#'        Endpoint 2.  Required when \code{prob = 'posterior'};
#'        otherwise set to \code{NULL}.
#' @param theta_NULL1 A numeric scalar giving the null hypothesis
#'        threshold for Endpoint 1.  Required when
#'        \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param theta_NULL2 A numeric scalar giving the null hypothesis
#'        threshold for Endpoint 2.  Required when
#'        \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param m_t A positive integer giving the future sample size for
#'        the treatment group.  Required when \code{prob = 'predictive'};
#'        set to \code{NULL} otherwise.
#' @param m_c A positive integer giving the future sample size for
#'        the control group.  Required when \code{prob = 'predictive'};
#'        set to \code{NULL} otherwise.
#' @param z00 A non-negative integer giving the hypothetical control
#'        count for pattern (0,0).  Required when
#'        \code{design = 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param z01 A non-negative integer; see \code{z00}.
#' @param z10 A non-negative integer; see \code{z00}.
#' @param z11 A non-negative integer; see \code{z00}.
#' @param xe_t_00 A non-negative integer giving the external treatment group
#'        count for pattern (0,0).  Required when
#'        \code{design = 'external'} and external treatment data are
#'        used; otherwise \code{NULL}.
#' @param xe_t_01 A non-negative integer; see \code{xe_t_00}.
#' @param xe_t_10 A non-negative integer; see \code{xe_t_00}.
#' @param xe_t_11 A non-negative integer; see \code{xe_t_00}.
#' @param xe_c_00 A non-negative integer giving the external control group
#'        count for pattern (0,0).  Required when
#'        \code{design = 'external'} and external control data are
#'        used; otherwise \code{NULL}.
#' @param xe_c_01 A non-negative integer; see \code{xe_c_00}.
#' @param xe_c_10 A non-negative integer; see \code{xe_c_00}.
#' @param xe_c_11 A non-negative integer; see \code{xe_c_00}.
#' @param ae_t A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for external treatment group data.  Required when external
#'        treatment data are used; otherwise \code{NULL}.
#' @param ae_c A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for external control group data.  Required when external
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
#' @param gamma_go_grid A numeric vector of candidate Go threshold values
#'        in \code{(0, 1)} to search over.  Defaults to
#'        \code{seq(0.01, 0.99, by = 0.01)}.
#' @param gamma_nogo_grid A numeric vector of candidate NoGo threshold
#'        values in \code{(0, 1)} to search over.  Defaults to
#'        \code{seq(0.01, 0.99, by = 0.01)}.
#'
#' @return A list of class \code{getgamma2bin} with the following
#'   elements:
#' \describe{
#'   \item{gamma_go}{Optimal Go threshold selected from
#'         \code{gamma_go_grid} according to \code{crit_go} and
#'         \code{sel_go}.  \code{NA} if no value satisfies the
#'         criterion.}
#'   \item{gamma_nogo}{Optimal NoGo threshold selected from
#'         \code{gamma_nogo_grid} according to \code{crit_nogo} and
#'         \code{sel_nogo}.  \code{NA} if no value satisfies the
#'         criterion.}
#'   \item{PrGo_at_gamma}{Pr(Go) evaluated at
#'         \code{(gamma_go, gamma_nogo)}.  \code{NA} if either threshold
#'         is \code{NA}.}
#'   \item{PrNoGo_at_gamma}{Pr(NoGo) evaluated at
#'         \code{(gamma_go, gamma_nogo)}.  \code{NA} if either threshold
#'         is \code{NA}.}
#'   \item{gamma_go_grid}{The candidate Go threshold grid used.}
#'   \item{gamma_nogo_grid}{The candidate NoGo threshold grid used.}
#'   \item{PrGo_grid}{Numeric matrix of dimensions
#'         \code{length(gamma_go_grid)} x \code{length(gamma_nogo_grid)}
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
#' decision thresholds \eqn{(\gamma_go, \gamma_nogo)}.
#'
#' \strong{Stage 2 (gamma sweep)}: For each pair
#' \eqn{(\gamma_go, \gamma_nogo)} in the two-dimensional grid, the
#' operating characteristics are computed as:
#' \deqn{\Pr(\mathrm{Go}) = \sum_{i,j} w_{ij}
#'   \mathbf{1}\!\left[\hat{g}_{Go,ij} \ge \gamma_go,\;
#'   \hat{g}_{NoGo,ij} < \gamma_nogo\right]}
#' \deqn{\Pr(\mathrm{NoGo}) = \sum_{i,j} w_{ij}
#'   \mathbf{1}\!\left[\hat{g}_{NoGo,ij} \ge \gamma_nogo,\;
#'   \hat{g}_{Go,ij} < \gamma_go\right]}
#' where \eqn{w_{ij} = P(X_t = x_{t,i}) \times P(X_c = x_{c,j})} are
#' the multinomial probability weights.  No additional Dirichlet draws
#' are required for this step.
#'
#' \strong{Stage 3 (optimal threshold selection)}: For each candidate
#' \eqn{\gamma_go}, the maximum \eqn{\Pr(\mathrm{Go})} over all
#' \eqn{\gamma_nogo} in \code{gamma_nogo_grid} is computed; the optimal
#' \eqn{\gamma_go} is then chosen as the \code{sel_go} value satisfying
#' the \code{crit_go} comparison against \code{target_go}.  Analogously,
#' the optimal \eqn{\gamma_nogo} is the \code{sel_nogo} value for which
#' the maximum \eqn{\Pr(\mathrm{NoGo})} over all \eqn{\gamma_go}
#' satisfies \code{crit_nogo} against \code{target_nogo}.
#'
#' @examples
#' # Example 1: Controlled design, posterior probability
#' # gamma_go: smallest gamma_go s.t. max_{gamma_nogo} Pr(Go)   < 0.05
#' # gamma_nogo: largest  gamma_nogo s.t. max_{gamma_go} Pr(NoGo) < 0.20
#' \dontrun{
#' getgamma2bin(
#'   prob = 'posterior', design = 'controlled',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   pi_t1 = 0.15, pi_t2 = 0.20, rho_t = 0.0,
#'   pi_c1 = 0.15, pi_c2 = 0.20, rho_c = 0.0,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n_t = 10L, n_c = 10L,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   theta_TV1   = 0.15, theta_MAV1 = 0.10,
#'   theta_TV2   = 0.15, theta_MAV2 = 0.10,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   m_t = NULL, m_c = NULL,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   ae_t = NULL, ae_c = NULL,
#'   nMC = 100L, method = 'Exact',
#'   gamma_go_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma_nogo_grid = seq(0.01, 0.99, by = 0.01)
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
#'   n_t = 10L, n_c = 10L,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   theta_TV1   = 0.15, theta_MAV1 = 0.10,
#'   theta_TV2   = 0.15, theta_MAV2 = 0.10,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   m_t = NULL, m_c = NULL,
#'   z00 = 3L, z01 = 2L, z10 = 3L, z11 = 2L,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   ae_t = NULL, ae_c = NULL,
#'   nMC = 100L, method = 'Exact',
#'   gamma_go_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma_nogo_grid = seq(0.01, 0.99, by = 0.01)
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
#'   n_t = 10L, n_c = 10L,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   theta_TV1   = NULL, theta_MAV1 = NULL,
#'   theta_TV2   = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 0.10, theta_NULL2 = 0.10,
#'   m_t = 30L, m_c = 30L,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   ae_t = NULL, ae_c = NULL,
#'   nMC = 100L, method = 'Exact',
#'   gamma_go_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma_nogo_grid = seq(0.01, 0.99, by = 0.01)
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
#'   n_t = 10L, n_c = 10L,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   theta_TV1   = NULL, theta_MAV1 = NULL,
#'   theta_TV2   = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 0.10, theta_NULL2 = 0.10,
#'   m_t = 30L, m_c = 30L,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe_t_00 = 3L, xe_t_01 = 2L, xe_t_10 = 3L, xe_t_11 = 2L,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   ae_t = 0.5, ae_c = NULL,
#'   nMC = 100L, method = 'Exact',
#'   gamma_go_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma_nogo_grid = seq(0.01, 0.99, by = 0.01)
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
                         n_t, n_c,
                         a_t_00, a_t_01, a_t_10, a_t_11,
                         a_c_00, a_c_01, a_c_10, a_c_11,
                         theta_TV1   = NULL, theta_MAV1  = NULL,
                         theta_TV2   = NULL, theta_MAV2  = NULL,
                         theta_NULL1 = NULL, theta_NULL2 = NULL,
                         m_t = NULL, m_c = NULL,
                         z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
                         xe_t_00 = NULL, xe_t_01 = NULL,
                         xe_t_10 = NULL, xe_t_11 = NULL,
                         xe_c_00 = NULL, xe_c_01 = NULL,
                         xe_c_10 = NULL, xe_c_11 = NULL,
                         ae_t = NULL, ae_c = NULL,
                         nMC    = 1000L,
                         method = 'Exact',
                         nsim   = 10000L,
                         gamma_go_grid = seq(0.01, 0.99, by = 0.01),
                         gamma_nogo_grid = seq(0.01, 0.99, by = 0.01)) {

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

  for (nm in c("n_t", "n_c")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val != floor(val) || val < 1L) {
      stop(paste0("'", nm, "' must be a single positive integer"))
    }
  }

  for (nm in c("a_t_00", "a_t_01", "a_t_10", "a_t_11",
               "a_c_00", "a_c_01", "a_c_10", "a_c_11")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0) {
      stop(paste0("'", nm, "' must be a single positive numeric value"))
    }
  }

  if (prob == 'posterior') {
    for (nm in c("theta_TV1", "theta_MAV1", "theta_TV2", "theta_MAV2")) {
      val <- get(nm)
      if (is.null(val)) {
        stop(paste0("'", nm, "' must be non-NULL when prob = 'posterior'"))
      }
    }
    if (theta_TV1 <= theta_MAV1) {
      stop("'theta_TV1' must be strictly greater than 'theta_MAV1'")
    }
    if (theta_TV2 <= theta_MAV2) {
      stop("'theta_TV2' must be strictly greater than 'theta_MAV2'")
    }
  } else {
    for (nm in c("theta_NULL1", "theta_NULL2")) {
      if (is.null(get(nm))) {
        stop(paste0("'", nm, "' must be non-NULL when prob = 'predictive'"))
      }
    }
    if (is.null(m_t) || !is.numeric(m_t) || length(m_t) != 1L || is.na(m_t) ||
        m_t != floor(m_t) || m_t < 1L) {
      stop("'m_t' must be a single positive integer when prob = 'predictive'")
    }
    if (is.null(m_c) || !is.numeric(m_c) || length(m_c) != 1L || is.na(m_c) ||
        m_c != floor(m_c) || m_c < 1L) {
      stop("'m_c' must be a single positive integer when prob = 'predictive'")
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

  for (gname in c("gamma_go_grid", "gamma_nogo_grid")) {
    gval <- get(gname)
    if (!is.numeric(gval) || length(gval) < 1L ||
        any(is.na(gval)) || any(gval <= 0) || any(gval >= 1)) {
      stop(paste0("'", gname, "' must be a numeric vector with all values in (0, 1)"))
    }
  }
  gamma_go_grid <- sort(unique(gamma_go_grid))
  gamma_nogo_grid <- sort(unique(gamma_nogo_grid))
  ng_go         <- length(gamma_go_grid)
  ng_nogo         <- length(gamma_nogo_grid)

  # ---------------------------------------------------------------------------
  # Stage 1: Enumerate all count combinations, compute PrGo_hat and
  #          PrNoGo_hat per combination via pbayespostpred2bin
  # ---------------------------------------------------------------------------
  counts_t <- allmultinom(n_t)   # (n_row_t x 4) integer matrix
  n_row_t      <- nrow(counts_t)

  if (design == 'uncontrolled') {
    # Control distribution is fixed (determined by z00..z11); only
    # treatment outcomes vary.
    counts_c <- matrix(0L, nrow = 1L, ncol = 4L)
    n_row_c      <- 1L
  } else {
    counts_c <- allmultinom(n_c)
    n_row_c      <- nrow(counts_c)
  }

  # Multinomial probability weights for treatment outcomes
  p_t     <- getjointbin(pi1 = pi_t1, pi2 = pi_t2, rho = rho_t)
  log_w_t <- apply(counts_t, 1L, function(x) {
    dmultinom(x, size = n_t, prob = p_t, log = TRUE)
  })

  if (design != 'uncontrolled') {
    p_c     <- getjointbin(pi1 = pi_c1, pi2 = pi_c2, rho = rho_c)
    log_w_c <- apply(counts_c, 1L, function(x) {
      dmultinom(x, size = n_c, prob = p_c, log = TRUE)
    })
    # w_mat[i, j] = P(X_t = counts_t[i,]) * P(X_c = counts_c[j,])
    w_mat <- exp(outer(log_w_t, log_w_c, `+`))
  } else {
    w_mat <- matrix(exp(log_w_t), nrow = n_row_t, ncol = 1L)
  }

  # Preallocate region probability storage: PrGo_hat[i,j] and PrNoGo_hat[i,j]
  PrGo_hat   <- matrix(0, nrow = n_row_t, ncol = n_row_c)
  PrNoGo_hat <- matrix(0, nrow = n_row_t, ncol = n_row_c)

  for (i in seq_len(n_row_t)) {
    x_t_i <- as.integer(counts_t[i, ])

    if (design == 'uncontrolled') {
      # Control counts not used: pass x2 as NULL
      Pr_R <- pbayespostpred2bin(
        prob    = prob, design  = design,
        theta_TV1   = theta_TV1,   theta_MAV1  = theta_MAV1,
        theta_TV2   = theta_TV2,   theta_MAV2  = theta_MAV2,
        theta_NULL1 = theta_NULL1, theta_NULL2 = theta_NULL2,
        x_t_00 = x_t_i[1L], x_t_01 = x_t_i[2L],
        x_t_10 = x_t_i[3L], x_t_11 = x_t_i[4L],
        x_c_00 = NULL, x_c_01 = NULL, x_c_10 = NULL, x_c_11 = NULL,
        a_t_00 = a_t_00, a_t_01 = a_t_01, a_t_10 = a_t_10, a_t_11 = a_t_11,
        a_c_00 = a_c_00, a_c_01 = a_c_01, a_c_10 = a_c_10, a_c_11 = a_c_11,
        m_t = m_t, m_c = m_c,
        z00 = z00, z01 = z01, z10 = z10, z11 = z11,
        xe_t_00 = xe_t_00, xe_t_01 = xe_t_01,
        xe_t_10 = xe_t_10, xe_t_11 = xe_t_11,
        xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
        ae_t = ae_t, ae_c = NULL,
        nMC = nMC
      )
      PrGo_hat[i, 1L]   <- sum(Pr_R[GoRegions])
      PrNoGo_hat[i, 1L] <- sum(Pr_R[NoGoRegions])

    } else {

      for (j in seq_len(n_row_c)) {
        x_c_j <- as.integer(counts_c[j, ])

        Pr_R <- pbayespostpred2bin(
          prob    = prob, design  = design,
          theta_TV1   = theta_TV1,   theta_MAV1  = theta_MAV1,
          theta_TV2   = theta_TV2,   theta_MAV2  = theta_MAV2,
          theta_NULL1 = theta_NULL1, theta_NULL2 = theta_NULL2,
          x_t_00 = x_t_i[1L], x_t_01 = x_t_i[2L],
          x_t_10 = x_t_i[3L], x_t_11 = x_t_i[4L],
          x_c_00 = x_c_j[1L], x_c_01 = x_c_j[2L],
          x_c_10 = x_c_j[3L], x_c_11 = x_c_j[4L],
          a_t_00 = a_t_00, a_t_01 = a_t_01, a_t_10 = a_t_10, a_t_11 = a_t_11,
          a_c_00 = a_c_00, a_c_01 = a_c_01, a_c_10 = a_c_10, a_c_11 = a_c_11,
          m_t = m_t, m_c = m_c,
          z00 = z00, z01 = z01, z10 = z10, z11 = z11,
          xe_t_00 = xe_t_00, xe_t_01 = xe_t_01,
          xe_t_10 = xe_t_10, xe_t_11 = xe_t_11,
          xe_c_00 = xe_c_00, xe_c_01 = xe_c_01,
          xe_c_10 = xe_c_10, xe_c_11 = xe_c_11,
          ae_t = ae_t, ae_c = ae_c,
          nMC = nMC
        )
        PrGo_hat[i, j]   <- sum(Pr_R[GoRegions])
        PrNoGo_hat[i, j] <- sum(Pr_R[NoGoRegions])
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Stage 2: Sweep (gamma_go_grid x gamma_nogo_grid)
  #
  # For each (g1, g2):
  #   Pr(Go)   = sum w[i,j] * I(PrGo_hat[i,j] >= g1  AND PrNoGo_hat[i,j] <  g2)
  #   Pr(NoGo) = sum w[i,j] * I(PrNoGo_hat[i,j] >= g2 AND PrGo_hat[i,j]  <  g1)
  # ---------------------------------------------------------------------------
  PrGo_grid   <- matrix(NA_real_, nrow = ng_go, ncol = ng_nogo)
  PrNoGo_grid <- matrix(NA_real_, nrow = ng_go, ncol = ng_nogo)

  for (k1 in seq_len(ng_go)) {
    g1 <- gamma_go_grid[k1]
    for (k2 in seq_len(ng_nogo)) {
      g2 <- gamma_nogo_grid[k2]

      go_mask   <- (PrGo_hat   >= g1) & (PrNoGo_hat <  g2)
      nogo_mask <- (PrNoGo_hat >= g2) & (PrGo_hat   <  g1)
      go_mask[is.na(go_mask)]     <- FALSE
      nogo_mask[is.na(nogo_mask)] <- FALSE

      PrGo_grid[k1, k2]   <- sum(w_mat[go_mask])
      PrNoGo_grid[k1, k2] <- sum(w_mat[nogo_mask])
    }
  }

  # ---------------------------------------------------------------------------
  # Stage 3: Select optimal (gamma_go, gamma_nogo)
  #
  # gamma_go: for each gamma_go candidate, take max Pr(Go) over gamma_nogo_grid;
  #         select the sel_go index satisfying crit_go vs target_go.
  # gamma_nogo: for each gamma_nogo candidate, take max Pr(NoGo) over gamma_go_grid;
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

  # Worst-case Pr(Go) over gamma_nogo_grid for each gamma_go candidate
  max_PrGo_per_g1   <- apply(PrGo_grid,   1L, max, na.rm = TRUE)
  opt1 <- .select_idx(.compare(max_PrGo_per_g1,   crit_go,   target_go),   sel_go)

  # Worst-case Pr(NoGo) over gamma_go_grid for each gamma_nogo candidate
  max_PrNoGo_per_g2 <- apply(PrNoGo_grid, 2L, max, na.rm = TRUE)
  opt2 <- .select_idx(.compare(max_PrNoGo_per_g2, crit_nogo, target_nogo), sel_nogo)

  if (is.na(opt1) || is.na(opt2)) {
    gamma_go          <- if (!is.na(opt1)) gamma_go_grid[opt1] else NA_real_
    gamma_nogo          <- if (!is.na(opt2)) gamma_nogo_grid[opt2] else NA_real_
    PrGo_at_gamma   <- NA_real_
    PrNoGo_at_gamma <- NA_real_
  } else {
    gamma_go          <- gamma_go_grid[opt1]
    gamma_nogo          <- gamma_nogo_grid[opt2]
    PrGo_at_gamma   <- PrGo_grid[opt1, opt2]
    PrNoGo_at_gamma <- PrNoGo_grid[opt1, opt2]
  }

  # ---------------------------------------------------------------------------
  # Build and return result
  # ---------------------------------------------------------------------------
  result <- list(
    gamma_go          = gamma_go,
    gamma_nogo          = gamma_nogo,
    PrGo_at_gamma   = PrGo_at_gamma,
    PrNoGo_at_gamma = PrNoGo_at_gamma,
    gamma_go_grid     = gamma_go_grid,
    gamma_nogo_grid     = gamma_nogo_grid,
    PrGo_grid       = PrGo_grid,
    PrNoGo_grid     = PrNoGo_grid
  )

  class(result) <- 'getgamma2bin'
  return(result)
}
