#' Find Optimal Go/NoGo Thresholds for Two Binary Endpoints
#'
#' Computes the optimal Go threshold \eqn{\gamma_{\mathrm{go}}} and NoGo
#' threshold \eqn{\gamma_{\mathrm{nogo}}} for two binary endpoints by
#' searching over a two-dimensional grid of candidate value pairs.  The two
#' thresholds are calibrated independently under separate scenarios:
#' \itemize{
#'   \item \eqn{\gamma_{\mathrm{go}}} is the \strong{smallest} value in
#'         \code{gamma_go_grid} such that the worst-case marginal Go
#'         probability over all \eqn{\gamma_{\mathrm{nogo}}} in
#'         \code{gamma_nogo_grid} is strictly less than \code{target_go}
#'         under the Go-calibration scenario (\code{pi_t1_go},
#'         \code{pi_t2_go}, \code{rho_t_go}, \code{pi_c1_go},
#'         \code{pi_c2_go}, \code{rho_c_go}); typically the Null scenario.
#'   \item \eqn{\gamma_{\mathrm{nogo}}} is the \strong{smallest} value in
#'         \code{gamma_nogo_grid} such that the worst-case marginal NoGo
#'         probability over all \eqn{\gamma_{\mathrm{go}}} in
#'         \code{gamma_go_grid} is strictly less than \code{target_nogo}
#'         under the NoGo-calibration scenario (\code{pi_t1_nogo},
#'         \code{pi_t2_nogo}, \code{rho_t_nogo}, \code{pi_c1_nogo},
#'         \code{pi_c2_nogo}, \code{rho_c_nogo}); typically the Alternative
#'         scenario.
#' }
#'
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or
#'        \code{'external'}.
#' @param GoRegions An integer vector of region indices (subset of
#'        \code{1:9} for \code{prob = 'posterior'} or \code{1:4} for
#'        \code{prob = 'predictive'}) that constitute the Go region.
#' @param NoGoRegions An integer vector of region indices that constitute
#'        the NoGo region.  Must be disjoint from \code{GoRegions}.
#' @param pi_t1_go A numeric scalar in \code{(0, 1)} giving the true
#'        treatment response probability for Endpoint 1 under the
#'        Go-calibration scenario (typically Null).
#' @param pi_t2_go A numeric scalar in \code{(0, 1)} giving the true
#'        treatment response probability for Endpoint 2 under the
#'        Go-calibration scenario.
#' @param rho_t_go A numeric scalar giving the within-group correlation in
#'        the treatment group under the Go-calibration scenario.
#' @param pi_c1_go A numeric scalar in \code{(0, 1)} giving the true
#'        control response probability for Endpoint 1 under the
#'        Go-calibration scenario.  Required for
#'        \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param pi_c2_go A numeric scalar in \code{(0, 1)} giving the true
#'        control response probability for Endpoint 2 under the
#'        Go-calibration scenario.  Required for
#'        \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param rho_c_go A numeric scalar giving the within-group correlation in
#'        the control group under the Go-calibration scenario.  Required for
#'        \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param pi_t1_nogo A numeric scalar in \code{(0, 1)} giving the true
#'        treatment response probability for Endpoint 1 under the
#'        NoGo-calibration scenario (typically Alternative).
#' @param pi_t2_nogo A numeric scalar in \code{(0, 1)} giving the true
#'        treatment response probability for Endpoint 2 under the
#'        NoGo-calibration scenario.
#' @param rho_t_nogo A numeric scalar giving the within-group correlation in
#'        the treatment group under the NoGo-calibration scenario.
#' @param pi_c1_nogo A numeric scalar in \code{(0, 1)} giving the true
#'        control response probability for Endpoint 1 under the
#'        NoGo-calibration scenario.  Required for
#'        \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param pi_c2_nogo A numeric scalar in \code{(0, 1)} giving the true
#'        control response probability for Endpoint 2 under the
#'        NoGo-calibration scenario.  Required for
#'        \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param rho_c_nogo A numeric scalar giving the within-group correlation in
#'        the control group under the NoGo-calibration scenario.  Required for
#'        \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param target_go A numeric scalar in \code{(0, 1)} giving the upper bound
#'        on the worst-case marginal Go probability under the Go-calibration
#'        scenario.  The optimal \eqn{\gamma_{\mathrm{go}}} is the smallest
#'        grid value satisfying the constraint.
#' @param target_nogo A numeric scalar in \code{(0, 1)} giving the upper
#'        bound on the worst-case marginal NoGo probability under the
#'        NoGo-calibration scenario.  The optimal \eqn{\gamma_{\mathrm{nogo}}}
#'        is the smallest grid value satisfying the constraint.
#' @param n_t A positive integer giving the number of patients in the
#'        treatment group in the PoC trial.
#' @param n_c A positive integer giving the number of patients in the
#'        control group in the PoC trial.
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
#' @param theta_NULL1 A numeric scalar giving the null hypothesis threshold
#'        for Endpoint 1.  Required when \code{prob = 'predictive'};
#'        otherwise set to \code{NULL}.
#' @param theta_NULL2 A numeric scalar giving the null hypothesis threshold
#'        for Endpoint 2.  Required when \code{prob = 'predictive'};
#'        otherwise set to \code{NULL}.
#' @param m_t A positive integer giving the future sample size for the
#'        treatment group.  Required when \code{prob = 'predictive'};
#'        set to \code{NULL} otherwise.
#' @param m_c A positive integer giving the future sample size for the
#'        control group.  Required when \code{prob = 'predictive'};
#'        set to \code{NULL} otherwise.
#' @param z00 A non-negative integer giving the hypothetical control count
#'        for pattern (0,0).  Required when
#'        \code{design = 'uncontrolled'}; otherwise set to \code{NULL}.
#' @param z01 A non-negative integer; see \code{z00}.
#' @param z10 A non-negative integer; see \code{z00}.
#' @param z11 A non-negative integer; see \code{z00}.
#' @param xe_t_00 A non-negative integer giving the external treatment group
#'        count for pattern (0,0).  Required when
#'        \code{design = 'external'} and external treatment data are used;
#'        otherwise \code{NULL}.
#' @param xe_t_01 A non-negative integer; see \code{xe_t_00}.
#' @param xe_t_10 A non-negative integer; see \code{xe_t_00}.
#' @param xe_t_11 A non-negative integer; see \code{xe_t_00}.
#' @param xe_c_00 A non-negative integer giving the external control group
#'        count for pattern (0,0).  Required when
#'        \code{design = 'external'} and external control data are used;
#'        otherwise \code{NULL}.
#' @param xe_c_01 A non-negative integer; see \code{xe_c_00}.
#' @param xe_c_10 A non-negative integer; see \code{xe_c_00}.
#' @param xe_c_11 A non-negative integer; see \code{xe_c_00}.
#' @param alpha0e_t A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for external treatment group data.  Required when external
#'        treatment data are used; otherwise \code{NULL}.
#' @param alpha0e_c A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for external control group data.  Required when external
#'        control data are used; otherwise \code{NULL}.
#' @param nMC A positive integer giving the number of Dirichlet draws used
#'        to evaluate region probabilities for each count combination in
#'        Stage 1.  Default is \code{1000L}.
#' @param gamma_go_grid A numeric vector of candidate Go threshold values
#'        in \code{(0, 1)} to search over.  Defaults to
#'        \code{seq(0.01, 0.99, by = 0.01)}.
#' @param gamma_nogo_grid A numeric vector of candidate NoGo threshold
#'        values in \code{(0, 1)} to search over.  Defaults to
#'        \code{seq(0.01, 0.99, by = 0.01)}.
#'
#' @return A list of class \code{getgamma2bin} with the following elements:
#' \describe{
#'   \item{gamma_go}{Optimal Go threshold: the smallest value in
#'         \code{gamma_go_grid} for which the marginal
#'         \eqn{\Pr(\mathrm{Go}) < \code{target\_go}} under the
#'         Go-calibration scenario.  \code{NA} if no such value exists.}
#'   \item{gamma_nogo}{Optimal NoGo threshold: the smallest value in
#'         \code{gamma_nogo_grid} for which the marginal
#'         \eqn{\Pr(\mathrm{NoGo}) < \code{target\_nogo}} under the
#'         NoGo-calibration scenario.  \code{NA} if no such value exists.}
#'   \item{PrGo_opt}{Marginal \eqn{\Pr(\mathrm{Go})} at
#'         \code{gamma_go} under the Go-calibration scenario.
#'         \code{NA} if \code{gamma_go} is \code{NA}.}
#'   \item{PrNoGo_opt}{Marginal \eqn{\Pr(\mathrm{NoGo})} at
#'         \code{gamma_nogo} under the NoGo-calibration scenario.
#'         \code{NA} if \code{gamma_nogo} is \code{NA}.}
#'   \item{target_go}{The value of \code{target_go} supplied by the user.}
#'   \item{target_nogo}{The value of \code{target_nogo} supplied by the user.}
#'   \item{grid_results}{A data frame with columns \code{gamma_grid},
#'         \code{PrGo_grid} (marginal Go probability under the Go-calibration
#'         scenario), and \code{PrNoGo_grid} (marginal NoGo probability under
#'         the NoGo-calibration scenario).}
#' }
#'
#' @details
#' The function uses the same two-stage precompute-then-sweep strategy as
#' \code{\link{pbayesdecisionprob2bin}}.
#'
#' \strong{Stage 1 (precomputation)}: \code{\link{pbayespostpred2bin}} is
#' called for every possible multinomial outcome combination \eqn{(x_t, x_c)}
#' enumerated by \code{\link{allmultinom}}.  The resulting region probability
#' vector is summed over \code{GoRegions} and \code{NoGoRegions} to obtain
#' \eqn{\hat{g}_{Go,ij}} and \eqn{\hat{g}_{NoGo,ij}}.  These are independent
#' of the calibration scenario; only the multinomial weights differ.
#'
#' \strong{Stage 2 (gamma sweep)}: For each pair
#' \eqn{(\gamma_{\mathrm{go}}, \gamma_{\mathrm{nogo}})} in the
#' two-dimensional grid, operating characteristics are computed separately
#' under each calibration scenario using the respective multinomial weights:
#' \deqn{\Pr(\mathrm{Go}) = \sum_{i,j} w_{ij}^{(\mathrm{go})}
#'   \mathbf{1}\!\left[\hat{g}_{Go,ij} \ge \gamma_{\mathrm{go}},\;
#'   \hat{g}_{NoGo,ij} < \gamma_{\mathrm{nogo}}\right]}
#' \deqn{\Pr(\mathrm{NoGo}) = \sum_{i,j} w_{ij}^{(\mathrm{nogo})}
#'   \mathbf{1}\!\left[\hat{g}_{NoGo,ij} \ge \gamma_{\mathrm{nogo}},\;
#'   \hat{g}_{Go,ij} < \gamma_{\mathrm{go}}\right]}
#'
#' \strong{Stage 3 (optimal threshold selection)}: For each candidate
#' \eqn{\gamma_{\mathrm{go}}}, the worst-case \eqn{\Pr(\mathrm{Go})} over all
#' \eqn{\gamma_{\mathrm{nogo}}} in \code{gamma_nogo_grid} is computed; the
#' optimal \eqn{\gamma_{\mathrm{go}}} is the \emph{smallest} grid value for
#' which this worst-case probability is less than \code{target_go}.
#' Analogously, the optimal \eqn{\gamma_{\mathrm{nogo}}} is the
#' \emph{smallest} grid value for which the worst-case \eqn{\Pr(\mathrm{NoGo})}
#' is less than \code{target_nogo}.
#'
#' @examples
#' # Example 1: Controlled design, posterior probability
#' # gamma_go  : smallest gamma_go   s.t. max_{gamma_nogo} Pr(Go)   < 0.05 under Null
#' # gamma_nogo: smallest gamma_nogo s.t. max_{gamma_go}   Pr(NoGo) < 0.20 under Alt
#' \dontrun{
#' getgamma2bin(
#'   prob = 'posterior', design = 'controlled',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   pi_t1_go = 0.15, pi_t2_go = 0.20, rho_t_go = 0.0,
#'   pi_c1_go = 0.15, pi_c2_go = 0.20, rho_c_go = 0.0,
#'   pi_t1_nogo = 0.35, pi_t2_nogo = 0.40, rho_t_nogo = 0.0,
#'   pi_c1_nogo = 0.15, pi_c2_nogo = 0.20, rho_c_nogo = 0.0,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 7L, n_c = 7L,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   theta_TV1   = 0.15, theta_MAV1 = 0.10,
#'   theta_TV2   = 0.15, theta_MAV2 = 0.10,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   m_t = NULL, m_c = NULL,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL,
#'   nMC = 100L,
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
#'   pi_t1_go = 0.15, pi_t2_go = 0.20, rho_t_go = 0.0,
#'   pi_c1_go = NULL, pi_c2_go = NULL, rho_c_go = NULL,
#'   pi_t1_nogo = 0.35, pi_t2_nogo = 0.40, rho_t_nogo = 0.0,
#'   pi_c1_nogo = NULL, pi_c2_nogo = NULL, rho_c_nogo = NULL,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 7L, n_c = 7L,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   theta_TV1   = 0.15, theta_MAV1 = 0.10,
#'   theta_TV2   = 0.15, theta_MAV2 = 0.10,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   m_t = NULL, m_c = NULL,
#'   z00 = 3L, z01 = 2L, z10 = 3L, z11 = 2L,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL,
#'   nMC = 100L,
#'   gamma_go_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma_nogo_grid = seq(0.01, 0.99, by = 0.01)
#' )
#' }
#'
#' # Example 3: External design, posterior probability
#' \dontrun{
#' getgamma2bin(
#'   prob = 'posterior', design = 'external',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   pi_t1_go = 0.15, pi_t2_go = 0.20, rho_t_go = 0.0,
#'   pi_c1_go = 0.15, pi_c2_go = 0.20, rho_c_go = 0.0,
#'   pi_t1_nogo = 0.35, pi_t2_nogo = 0.40, rho_t_nogo = 0.0,
#'   pi_c1_nogo = 0.15, pi_c2_nogo = 0.20, rho_c_nogo = 0.0,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 7L, n_c = 7L,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   theta_TV1   = 0.15, theta_MAV1 = 0.10,
#'   theta_TV2   = 0.15, theta_MAV2 = 0.10,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   m_t = NULL, m_c = NULL,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = 3L, xe_c_01 = 2L, xe_c_10 = 3L, xe_c_11 = 2L,
#'   alpha0e_t = NULL, alpha0e_c = 0.5,
#'   nMC = 100L,
#'   gamma_go_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma_nogo_grid = seq(0.01, 0.99, by = 0.01)
#' )
#' }
#'
#' # Example 4: Controlled design, predictive probability
#' \dontrun{
#' getgamma2bin(
#'   prob = 'predictive', design = 'controlled',
#'   GoRegions = 1L, NoGoRegions = 4L,
#'   pi_t1_go = 0.15, pi_t2_go = 0.20, rho_t_go = 0.0,
#'   pi_c1_go = 0.15, pi_c2_go = 0.20, rho_c_go = 0.0,
#'   pi_t1_nogo = 0.35, pi_t2_nogo = 0.40, rho_t_nogo = 0.0,
#'   pi_c1_nogo = 0.15, pi_c2_nogo = 0.20, rho_c_nogo = 0.0,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 7L, n_c = 7L,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   theta_TV1   = NULL, theta_MAV1 = NULL,
#'   theta_TV2   = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 0.10, theta_NULL2 = 0.10,
#'   m_t = 5L, m_c = 5L,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL,
#'   nMC = 100L,
#'   gamma_go_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma_nogo_grid = seq(0.01, 0.99, by = 0.01)
#' )
#' }
#'
#' # Example 5: Uncontrolled design, predictive probability
#' \dontrun{
#' getgamma2bin(
#'   prob = 'predictive', design = 'uncontrolled',
#'   GoRegions = 1L, NoGoRegions = 4L,
#'   pi_t1_go = 0.15, pi_t2_go = 0.20, rho_t_go = 0.0,
#'   pi_c1_go = NULL, pi_c2_go = NULL, rho_c_go = NULL,
#'   pi_t1_nogo = 0.35, pi_t2_nogo = 0.40, rho_t_nogo = 0.0,
#'   pi_c1_nogo = NULL, pi_c2_nogo = NULL, rho_c_nogo = NULL,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 7L, n_c = 7L,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   theta_TV1   = NULL, theta_MAV1 = NULL,
#'   theta_TV2   = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 0.10, theta_NULL2 = 0.10,
#'   m_t = 5L, m_c = 5L,
#'   z00 = 3L, z01 = 2L, z10 = 3L, z11 = 2L,
#'   xe_t_00 = NULL, xe_t_01 = NULL, xe_t_10 = NULL, xe_t_11 = NULL,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL,
#'   nMC = 100L,
#'   gamma_go_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma_nogo_grid = seq(0.01, 0.99, by = 0.01)
#' )
#' }
#'
#' # Example 6: External design, predictive probability
#' \dontrun{
#' getgamma2bin(
#'   prob = 'predictive', design = 'external',
#'   GoRegions = 1L, NoGoRegions = 4L,
#'   pi_t1_go = 0.15, pi_t2_go = 0.20, rho_t_go = 0.0,
#'   pi_c1_go = 0.15, pi_c2_go = 0.20, rho_c_go = 0.0,
#'   pi_t1_nogo = 0.35, pi_t2_nogo = 0.40, rho_t_nogo = 0.0,
#'   pi_c1_nogo = 0.15, pi_c2_nogo = 0.20, rho_c_nogo = 0.0,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 7L, n_c = 7L,
#'   a_t_00 = 0.25, a_t_01 = 0.25, a_t_10 = 0.25, a_t_11 = 0.25,
#'   a_c_00 = 0.25, a_c_01 = 0.25, a_c_10 = 0.25, a_c_11 = 0.25,
#'   theta_TV1   = NULL, theta_MAV1 = NULL,
#'   theta_TV2   = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 0.10, theta_NULL2 = 0.10,
#'   m_t = 5L, m_c = 5L,
#'   z00 = NULL, z01 = NULL, z10 = NULL, z11 = NULL,
#'   xe_t_00 = 3L, xe_t_01 = 2L, xe_t_10 = 3L, xe_t_11 = 2L,
#'   xe_c_00 = NULL, xe_c_01 = NULL, xe_c_10 = NULL, xe_c_11 = NULL,
#'   alpha0e_t = 0.5, alpha0e_c = NULL,
#'   nMC = 100L,
#'   gamma_go_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma_nogo_grid = seq(0.01, 0.99, by = 0.01)
#' )
#' }
#'
#' @importFrom stats dmultinom
#' @export
getgamma2bin <- function(prob = 'posterior', design = 'controlled',
                         GoRegions, NoGoRegions,
                         pi_t1_go, pi_t2_go, rho_t_go,
                         pi_c1_go = NULL, pi_c2_go = NULL, rho_c_go = NULL,
                         pi_t1_nogo, pi_t2_nogo, rho_t_nogo,
                         pi_c1_nogo = NULL, pi_c2_nogo = NULL, rho_c_nogo = NULL,
                         target_go, target_nogo,
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
                         alpha0e_t = NULL, alpha0e_c = NULL,
                         nMC    = 1000L,
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

  max_region <- if (prob == 'posterior') 9L else 4L
  if (!is.integer(GoRegions) || length(GoRegions) < 1L ||
      any(is.na(GoRegions)) || any(GoRegions < 1L) || any(GoRegions > max_region)) {
    stop(sprintf("'GoRegions' must be an integer vector with all values in 1:%d", max_region))
  }
  if (!is.integer(NoGoRegions) || length(NoGoRegions) < 1L ||
      any(is.na(NoGoRegions)) || any(NoGoRegions < 1L) ||
      any(NoGoRegions > max_region)) {
    stop(sprintf("'NoGoRegions' must be an integer vector with all values in 1:%d", max_region))
  }
  if (length(intersect(GoRegions, NoGoRegions)) > 0L) {
    stop("'GoRegions' and 'NoGoRegions' must be disjoint")
  }

  # Validate Go-calibration scenario parameters
  if (!is.numeric(pi_t1_go) || length(pi_t1_go) != 1L || is.na(pi_t1_go) ||
      pi_t1_go <= 0 || pi_t1_go >= 1) {
    stop("'pi_t1_go' must be a single numeric value in (0, 1)")
  }
  if (!is.numeric(pi_t2_go) || length(pi_t2_go) != 1L || is.na(pi_t2_go) ||
      pi_t2_go <= 0 || pi_t2_go >= 1) {
    stop("'pi_t2_go' must be a single numeric value in (0, 1)")
  }
  if (!is.numeric(rho_t_go) || length(rho_t_go) != 1L || is.na(rho_t_go)) {
    stop("'rho_t_go' must be a single numeric value")
  }

  # Validate NoGo-calibration scenario parameters
  if (!is.numeric(pi_t1_nogo) || length(pi_t1_nogo) != 1L || is.na(pi_t1_nogo) ||
      pi_t1_nogo <= 0 || pi_t1_nogo >= 1) {
    stop("'pi_t1_nogo' must be a single numeric value in (0, 1)")
  }
  if (!is.numeric(pi_t2_nogo) || length(pi_t2_nogo) != 1L || is.na(pi_t2_nogo) ||
      pi_t2_nogo <= 0 || pi_t2_nogo >= 1) {
    stop("'pi_t2_nogo' must be a single numeric value in (0, 1)")
  }
  if (!is.numeric(rho_t_nogo) || length(rho_t_nogo) != 1L || is.na(rho_t_nogo)) {
    stop("'rho_t_nogo' must be a single numeric value")
  }

  if (design != 'uncontrolled') {
    if (is.null(pi_c1_go) || !is.numeric(pi_c1_go) || length(pi_c1_go) != 1L ||
        is.na(pi_c1_go) || pi_c1_go <= 0 || pi_c1_go >= 1) {
      stop("'pi_c1_go' must be a single numeric value in (0, 1) for controlled or external design")
    }
    if (is.null(pi_c2_go) || !is.numeric(pi_c2_go) || length(pi_c2_go) != 1L ||
        is.na(pi_c2_go) || pi_c2_go <= 0 || pi_c2_go >= 1) {
      stop("'pi_c2_go' must be a single numeric value in (0, 1) for controlled or external design")
    }
    if (is.null(rho_c_go) || !is.numeric(rho_c_go) || length(rho_c_go) != 1L ||
        is.na(rho_c_go)) {
      stop("'rho_c_go' must be a single numeric value for controlled or external design")
    }
    if (is.null(pi_c1_nogo) || !is.numeric(pi_c1_nogo) || length(pi_c1_nogo) != 1L ||
        is.na(pi_c1_nogo) || pi_c1_nogo <= 0 || pi_c1_nogo >= 1) {
      stop("'pi_c1_nogo' must be a single numeric value in (0, 1) for controlled or external design")
    }
    if (is.null(pi_c2_nogo) || !is.numeric(pi_c2_nogo) || length(pi_c2_nogo) != 1L ||
        is.na(pi_c2_nogo) || pi_c2_nogo <= 0 || pi_c2_nogo >= 1) {
      stop("'pi_c2_nogo' must be a single numeric value in (0, 1) for controlled or external design")
    }
    if (is.null(rho_c_nogo) || !is.numeric(rho_c_nogo) || length(rho_c_nogo) != 1L ||
        is.na(rho_c_nogo)) {
      stop("'rho_c_nogo' must be a single numeric value for controlled or external design")
    }
  }

  for (nm in c("target_go", "target_nogo")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val <= 0 || val >= 1) {
      stop(paste0("'", nm, "' must be a single numeric value in (0, 1)"))
    }
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

  if (!is.numeric(nMC) || length(nMC) != 1L || is.na(nMC) ||
      nMC != floor(nMC) || nMC < 1L) {
    stop("'nMC' must be a single positive integer")
  }

  for (gname in c("gamma_go_grid", "gamma_nogo_grid")) {
    gval <- get(gname)
    if (!is.numeric(gval) || length(gval) < 1L ||
        any(is.na(gval)) || any(gval <= 0) || any(gval >= 1)) {
      stop(paste0("'", gname, "' must be a numeric vector with all values in (0, 1)"))
    }
  }
  gamma_go_grid   <- sort(unique(gamma_go_grid))
  gamma_nogo_grid <- sort(unique(gamma_nogo_grid))
  if (!identical(gamma_go_grid, gamma_nogo_grid)) {
    stop("'gamma_go_grid' and 'gamma_nogo_grid' must be identical vectors")
  }
  gamma_grid <- gamma_go_grid
  ng         <- length(gamma_grid)

  # ---------------------------------------------------------------------------
  # Stage 1: Enumerate all count combinations, compute PrGo_hat and
  #          PrNoGo_hat per combination via pbayespostpred2bin
  # ---------------------------------------------------------------------------
  counts_t <- allmultinom(n_t)   # (n_row_t x 4) integer matrix
  n_row_t  <- nrow(counts_t)

  if (design == 'uncontrolled') {
    # Control distribution is fixed (determined by z00..z11); only
    # treatment outcomes vary.
    counts_c <- matrix(0L, nrow = 1L, ncol = 4L)
    n_row_c  <- 1L
  } else {
    counts_c <- allmultinom(n_c)
    n_row_c  <- nrow(counts_c)
  }

  # Multinomial probability weights for Go-calibration scenario
  p_t_go     <- getjointbin(pi1 = pi_t1_go, pi2 = pi_t2_go, rho = rho_t_go)
  log_w_t_go <- apply(counts_t, 1L, function(x) {
    dmultinom(x, size = n_t, prob = p_t_go, log = TRUE)
  })

  # Multinomial probability weights for NoGo-calibration scenario
  p_t_nogo     <- getjointbin(pi1 = pi_t1_nogo, pi2 = pi_t2_nogo, rho = rho_t_nogo)
  log_w_t_nogo <- apply(counts_t, 1L, function(x) {
    dmultinom(x, size = n_t, prob = p_t_nogo, log = TRUE)
  })

  if (design != 'uncontrolled') {
    p_c_go     <- getjointbin(pi1 = pi_c1_go, pi2 = pi_c2_go, rho = rho_c_go)
    log_w_c_go <- apply(counts_c, 1L, function(x) {
      dmultinom(x, size = n_c, prob = p_c_go, log = TRUE)
    })
    p_c_nogo     <- getjointbin(pi1 = pi_c1_nogo, pi2 = pi_c2_nogo, rho = rho_c_nogo)
    log_w_c_nogo <- apply(counts_c, 1L, function(x) {
      dmultinom(x, size = n_c, prob = p_c_nogo, log = TRUE)
    })
    # w_mat_go[i,j]   = P_go(X_t = counts_t[i,])   * P_go(X_c = counts_c[j,])
    # w_mat_nogo[i,j] = P_nogo(X_t = counts_t[i,]) * P_nogo(X_c = counts_c[j,])
    w_mat_go   <- exp(outer(log_w_t_go,   log_w_c_go,   `+`))
    w_mat_nogo <- exp(outer(log_w_t_nogo, log_w_c_nogo, `+`))
  } else {
    w_mat_go   <- matrix(exp(log_w_t_go),   nrow = n_row_t, ncol = 1L)
    w_mat_nogo <- matrix(exp(log_w_t_nogo), nrow = n_row_t, ncol = 1L)
  }

  # Preallocate region probability storage: PrGo_hat[i,j] and PrNoGo_hat[i,j]
  # These are functions of data only (independent of calibration scenario weights)
  PrGo_hat   <- matrix(0, nrow = n_row_t, ncol = n_row_c)
  PrNoGo_hat <- matrix(0, nrow = n_row_t, ncol = n_row_c)

  for (i in seq_len(n_row_t)) {
    x_t_i <- as.integer(counts_t[i, ])

    if (design == 'uncontrolled') {
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
        alpha0e_t = alpha0e_t, alpha0e_c = NULL,
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
          alpha0e_t = alpha0e_t, alpha0e_c = alpha0e_c,
          nMC = nMC
        )
        PrGo_hat[i, j]   <- sum(Pr_R[GoRegions])
        PrNoGo_hat[i, j] <- sum(Pr_R[NoGoRegions])
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Stage 2: Marginal sweep over gamma_go_grid and gamma_nogo_grid
  #
  # Go and NoGo probabilities are computed marginally (independently):
  #   Pr(Go)   = sum w_go[i,j]   * I(PrGo_hat[i,j]   >= g1)
  #   Pr(NoGo) = sum w_nogo[i,j] * I(PrNoGo_hat[i,j] >= g2)
  # ---------------------------------------------------------------------------
  PrGo_grid   <- numeric(ng)
  PrNoGo_grid <- numeric(ng)

  for (k in seq_len(ng)) {
    g <- gamma_grid[k]
    go_mask <- (PrGo_hat >= g)
    go_mask[is.na(go_mask)] <- FALSE
    PrGo_grid[k] <- sum(w_mat_go[go_mask])
  }

  for (k in seq_len(ng)) {
    g <- gamma_grid[k]
    nogo_mask <- (PrNoGo_hat >= g)
    nogo_mask[is.na(nogo_mask)] <- FALSE
    PrNoGo_grid[k] <- sum(w_mat_nogo[nogo_mask])
  }

  # ---------------------------------------------------------------------------
  # Stage 3: Select optimal (gamma_go, gamma_nogo)
  #
  # gamma_go  : smallest value in gamma_go_grid   s.t. Pr(Go)   < target_go
  # gamma_nogo: smallest value in gamma_nogo_grid s.t. Pr(NoGo) < target_nogo
  # Both curves are non-increasing in their respective gamma.
  # ---------------------------------------------------------------------------
  idx_go <- which(PrGo_grid < target_go)
  if (length(idx_go) == 0L) {
    gamma_go <- NA_real_
    PrGo_opt <- NA_real_
  } else {
    opt1     <- min(idx_go)
    gamma_go <- gamma_grid[opt1]
    PrGo_opt <- PrGo_grid[opt1]
  }

  idx_nogo <- which(PrNoGo_grid < target_nogo)
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
    target_go    = target_go,
    target_nogo  = target_nogo,
    grid_results = data.frame(
      gamma_grid   = gamma_grid,
      PrGo_grid    = PrGo_grid,
      PrNoGo_grid  = PrNoGo_grid
    )
  )

  class(result) <- 'getgamma2bin'
  return(result)
}
