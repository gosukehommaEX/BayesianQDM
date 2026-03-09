#' Find Optimal Go/NoGo Thresholds for Two Continuous Endpoints
#'
#' Computes the optimal Go threshold \eqn{\gamma_{\mathrm{go}}} and NoGo
#' threshold \eqn{\gamma_{\mathrm{nogo}}} for two continuous endpoints by
#' searching over a two-dimensional grid of candidate value pairs.  The two
#' thresholds are calibrated independently under separate scenarios:
#' \itemize{
#'   \item \eqn{\gamma_{\mathrm{go}}} is the \strong{smallest} value in
#'         \code{gamma_go_grid} such that the worst-case marginal Go
#'         probability over all \eqn{\gamma_{\mathrm{nogo}}} in
#'         \code{gamma_nogo_grid} is strictly less than \code{target_go}
#'         under the Go-calibration scenario (\code{mu_t_go},
#'         \code{Sigma_t_go}, \code{mu_c_go}, \code{Sigma_c_go});
#'         typically the Null scenario.
#'   \item \eqn{\gamma_{\mathrm{nogo}}} is the \strong{smallest} value in
#'         \code{gamma_nogo_grid} such that the worst-case marginal NoGo
#'         probability over all \eqn{\gamma_{\mathrm{go}}} in
#'         \code{gamma_go_grid} is strictly less than \code{target_nogo}
#'         under the NoGo-calibration scenario (\code{mu_t_nogo},
#'         \code{Sigma_t_nogo}, \code{mu_c_nogo}, \code{Sigma_c_nogo});
#'         typically the Alternative scenario.
#' }
#'
#' @param nsim A positive integer giving the number of Monte Carlo
#'        datasets to simulate per calibration scenario.  Default is
#'        \code{10000L}.
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
#' @param mu_t_go A length-2 numeric vector giving the true bivariate mean
#'        for the treatment group under the Go-calibration scenario
#'        (typically Null).
#' @param Sigma_t_go A 2x2 positive-definite numeric matrix giving the
#'        true within-group covariance in the treatment group under the
#'        Go-calibration scenario.
#' @param mu_c_go A length-2 numeric vector giving the true bivariate mean
#'        for the control group under the Go-calibration scenario.  Required
#'        for \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param Sigma_c_go A 2x2 positive-definite numeric matrix giving the
#'        true within-group covariance in the control group under the
#'        Go-calibration scenario.  Required for
#'        \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param mu_t_nogo A length-2 numeric vector giving the true bivariate mean
#'        for the treatment group under the NoGo-calibration scenario
#'        (typically Alternative).
#' @param Sigma_t_nogo A 2x2 positive-definite numeric matrix giving the
#'        true within-group covariance in the treatment group under the
#'        NoGo-calibration scenario.
#' @param mu_c_nogo A length-2 numeric vector giving the true bivariate mean
#'        for the control group under the NoGo-calibration scenario.  Required
#'        for \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param Sigma_c_nogo A 2x2 positive-definite numeric matrix giving the
#'        true within-group covariance in the control group under the
#'        NoGo-calibration scenario.  Required for
#'        \code{design = 'controlled'} or \code{'external'};
#'        set to \code{NULL} for \code{design = 'uncontrolled'}.
#' @param target_go A numeric scalar in \code{(0, 1)} giving the upper bound
#'        on the worst-case marginal Go probability under the Go-calibration
#'        scenario.  The optimal \eqn{\gamma_{\mathrm{go}}} is the smallest
#'        grid value satisfying the constraint.
#' @param target_nogo A numeric scalar in \code{(0, 1)} giving the upper
#'        bound on the worst-case marginal NoGo probability under the
#'        NoGo-calibration scenario.  The optimal
#'        \eqn{\gamma_{\mathrm{nogo}}} is the smallest grid value satisfying
#'        the constraint.
#' @param n_t A positive integer giving the number of patients in the
#'        treatment group in the PoC trial.
#' @param n_c A positive integer giving the number of patients in the
#'        control group in the PoC trial.  Set to \code{NULL} for
#'        \code{design = 'uncontrolled'}.
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
#' @param kappa0_t Positive numeric scalar.  NIW prior hyperparameter
#'        \eqn{\kappa_{01}} for the treatment group.  Required when
#'        \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param nu0_t Positive numeric scalar.  NIW prior degrees of freedom
#'        \eqn{\nu_{01}} for the treatment group.  Required when
#'        \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param mu0_t Length-2 numeric vector.  NIW prior mean \eqn{\mu_{01}}
#'        for the treatment group.  Required when
#'        \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param Lambda0_t A 2x2 positive-definite numeric matrix.  NIW prior
#'        scale matrix \eqn{\Lambda_{01}} for the treatment group.  Required
#'        when \code{prior = 'N-Inv-Wishart'}; otherwise \code{NULL}.
#' @param kappa0_c Positive numeric scalar; see \code{kappa0_t}.
#'        For the control group.
#' @param nu0_c Positive numeric scalar; see \code{nu0_t}.
#'        For the control group.
#' @param mu0_c Length-2 numeric vector; see \code{mu0_t}.  For the control
#'        group.  May be required for the vague prior uncontrolled design;
#'        see \code{\link{pbayesdecisionprob2cont}}.
#' @param Lambda0_c A 2x2 matrix; see \code{Lambda0_t}.
#'        For the control group.
#' @param r A positive numeric scalar giving the power prior weight for
#'        the control group when \code{design = 'uncontrolled'} and
#'        \code{prior = 'vague'}.  Otherwise \code{NULL}.
#' @param ne_t A positive integer giving the external treatment sample size.
#'        Required when \code{design = 'external'} and external treatment
#'        data are used; otherwise \code{NULL}.
#' @param ne_c A positive integer giving the external control sample size.
#'        Required when \code{design = 'external'} and external control
#'        data are used; otherwise \code{NULL}.
#' @param alpha0e_t A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for external treatment data.  Required when external
#'        treatment data are used; otherwise \code{NULL}.
#' @param alpha0e_c A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for external control data.  Required when external control
#'        data are used; otherwise \code{NULL}.
#' @param bar_ye_t A length-2 numeric vector.  External treatment sample
#'        mean.  Required when external treatment data are used;
#'        otherwise \code{NULL}.
#' @param bar_ye_c A length-2 numeric vector.  External control sample mean.
#'        Required when external control data are used; otherwise \code{NULL}.
#' @param se_t A 2x2 numeric matrix.  External treatment sum-of-squares
#'        matrix.  Required when external treatment data are used;
#'        otherwise \code{NULL}.
#' @param se_c A 2x2 numeric matrix.  External control sum-of-squares
#'        matrix.  Required when external control data are used;
#'        otherwise \code{NULL}.
#' @param nMC A positive integer giving the number of Monte Carlo draws
#'        passed to \code{\link{pbayespostpred2cont}}.  Required when
#'        \code{CalcMethod = 'MC'}.  May be set to \code{NULL} when
#'        \code{CalcMethod = 'MM'} and \eqn{\nu_k > 4}; if
#'        \code{CalcMethod = 'MM'} but \eqn{\nu_k \le 4} causes a fallback
#'        to MC, \code{nMC} must be a positive integer.
#'        Default is \code{NULL}.
#' @param CalcMethod A character string specifying the computation method
#'        passed to \code{\link{pbayespostpred2cont}}.  Must be
#'        \code{'MC'} (default) or \code{'MM'}.
#' @param gamma_go_grid A numeric vector of candidate Go threshold values
#'        in \code{(0, 1)} to search over.  Defaults to
#'        \code{seq(0.01, 0.99, by = 0.01)}.
#' @param gamma_nogo_grid A numeric vector of candidate NoGo threshold
#'        values in \code{(0, 1)} to search over.  Defaults to
#'        \code{seq(0.01, 0.99, by = 0.01)}.
#' @param seed A numeric scalar for reproducible random number generation.
#'        The Go-calibration simulation uses \code{seed} and the
#'        NoGo-calibration simulation uses \code{seed + 1} to ensure
#'        independence between the two scenarios.
#'
#' @return A list of class \code{getgamma2cont} with the following elements:
#' \describe{
#'   \item{gamma_go}{Optimal Go threshold: the smallest value in
#'         \code{gamma_go_grid} for which the worst-case
#'         \eqn{\Pr(\mathrm{Go}) < \code{target\_go}} under the
#'         Go-calibration scenario.  \code{NA} if no such value exists.}
#'   \item{gamma_nogo}{Optimal NoGo threshold: the smallest value in
#'         \code{gamma_nogo_grid} for which the worst-case
#'         \eqn{\Pr(\mathrm{NoGo}) < \code{target\_nogo}} under the
#'         NoGo-calibration scenario.  \code{NA} if no such value exists.}
#'   \item{PrGo_opt}{Worst-case \eqn{\Pr(\mathrm{Go})} at
#'         \code{gamma_go} under the Go-calibration scenario.
#'         \code{NA} if \code{gamma_go} is \code{NA}.}
#'   \item{PrNoGo_opt}{Worst-case \eqn{\Pr(\mathrm{NoGo})} at
#'         \code{gamma_nogo} under the NoGo-calibration scenario.
#'         \code{NA} if \code{gamma_nogo} is \code{NA}.}
#'   \item{grid_results}{A list with elements \code{gamma_go_grid},
#'         \code{gamma_nogo_grid}, \code{PrGo_grid} (matrix of dimensions
#'         \code{length(gamma_go_grid)} x \code{length(gamma_nogo_grid)}
#'         under the Go-calibration scenario), and \code{PrNoGo_grid}
#'         (matrix of the same dimensions under the NoGo-calibration
#'         scenario).}
#' }
#'
#' @details
#' The function uses a two-stage simulate-then-sweep strategy:
#'
#' \strong{Stage 1 (simulation and precomputation)}: \code{nsim} bivariate
#' datasets are generated independently for each calibration scenario.
#' For the Go-calibration scenario, datasets are drawn from
#' \eqn{N_2(\mu_{t,\mathrm{go}}, \Sigma_{t,\mathrm{go}})} (and
#' \eqn{N_2(\mu_{c,\mathrm{go}}, \Sigma_{c,\mathrm{go}})} for
#' controlled/external designs); for the NoGo-calibration scenario,
#' the corresponding \code{_nogo} parameters are used.
#' \code{\link{pbayespostpred2cont}} is called once per scenario in
#' vectorised mode to return an \eqn{nsim \times 9} matrix of region
#' probabilities.  The probabilities are summed over \code{GoRegions}
#' (for the Go scenario) and \code{NoGoRegions} (for the NoGo scenario)
#' to obtain \eqn{\hat{g}_{Go,i}} and \eqn{\hat{g}_{NoGo,i}},
#' independent of the decision thresholds.
#'
#' \strong{Stage 2 (gamma sweep)}: For each pair
#' \eqn{(\gamma_{\mathrm{go}}, \gamma_{\mathrm{nogo}})} in the
#' two-dimensional grid, operating characteristics are computed separately
#' under each calibration scenario:
#' \deqn{\Pr(\mathrm{Go}) = \frac{1}{n_{\mathrm{sim}}} \sum_{i=1}^{n_{\mathrm{sim}}}
#'   \mathbf{1}\!\left[\hat{g}_{Go,i} \ge \gamma_{\mathrm{go}},\;
#'   \hat{g}_{NoGo,i} < \gamma_{\mathrm{nogo}}\right]}
#' \deqn{\Pr(\mathrm{NoGo}) = \frac{1}{n_{\mathrm{sim}}} \sum_{i=1}^{n_{\mathrm{sim}}}
#'   \mathbf{1}\!\left[\hat{g}_{NoGo,i} \ge \gamma_{\mathrm{nogo}},\;
#'   \hat{g}_{Go,i} < \gamma_{\mathrm{go}}\right]}
#'
#' \strong{Stage 3 (optimal threshold selection)}: For each candidate
#' \eqn{\gamma_{\mathrm{go}}}, the worst-case \eqn{\Pr(\mathrm{Go})} over
#' all \eqn{\gamma_{\mathrm{nogo}}} in \code{gamma_nogo_grid} is computed;
#' the optimal \eqn{\gamma_{\mathrm{go}}} is the \emph{smallest} grid value
#' for which this worst-case probability is less than \code{target_go}.
#' Analogously, the optimal \eqn{\gamma_{\mathrm{nogo}}} is the
#' \emph{smallest} grid value for which the worst-case
#' \eqn{\Pr(\mathrm{NoGo})} is less than \code{target_nogo}.
#'
#' @examples
#' # Example 1: Controlled design, posterior probability, vague prior
#' # gamma_go  : smallest gamma_go   s.t. max_{gamma_nogo} Pr(Go)   < 0.05 under Null
#' # gamma_nogo: smallest gamma_nogo s.t. max_{gamma_go}   Pr(NoGo) < 0.20 under Alt
#' \dontrun{
#' Sigma_null <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' Sigma_alt  <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' getgamma2cont(
#'   nsim = 1000L, prob = 'posterior', design = 'controlled',
#'   prior = 'vague',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   mu_t_go = c(-5.0, 0.0), Sigma_t_go = Sigma_null,
#'   mu_c_go = c(-10.0, -1.0), Sigma_c_go = Sigma_null,
#'   mu_t_nogo = c(5.0, 1.0), Sigma_t_nogo = Sigma_alt,
#'   mu_c_nogo = c(-10.0, -1.0), Sigma_c_nogo = Sigma_alt,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 30L, n_c = 30L,
#'   theta_TV1 = 10.0, theta_MAV1 = 5.0,
#'   theta_TV2 = 2.0,  theta_MAV2 = 1.0,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   m_t = NULL, m_c = NULL,
#'   kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
#'   kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
#'   r = NULL,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   nMC = 500L, CalcMethod = 'MC',
#'   gamma_go_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma_nogo_grid = seq(0.01, 0.99, by = 0.01),
#'   seed = 1L
#' )
#' }
#'
#' # Example 2: Uncontrolled design, posterior probability, vague prior
#' \dontrun{
#' Sigma_null <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' Sigma_alt  <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' getgamma2cont(
#'   nsim = 1000L, prob = 'posterior', design = 'uncontrolled',
#'   prior = 'vague',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   mu_t_go = c(-5.0, 0.0), Sigma_t_go = Sigma_null,
#'   mu_c_go = NULL, Sigma_c_go = NULL,
#'   mu_t_nogo = c(5.0, 1.0), Sigma_t_nogo = Sigma_alt,
#'   mu_c_nogo = NULL, Sigma_c_nogo = NULL,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 30L, n_c = NULL,
#'   theta_TV1 = 10.0, theta_MAV1 = 5.0,
#'   theta_TV2 = 2.0,  theta_MAV2 = 1.0,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   m_t = NULL, m_c = NULL,
#'   kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
#'   kappa0_c = NULL, nu0_c = NULL, mu0_c = c(-10.0, -1.0), Lambda0_c = NULL,
#'   r = 1.0,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   nMC = NULL, CalcMethod = 'MM',
#'   gamma_go_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma_nogo_grid = seq(0.01, 0.99, by = 0.01),
#'   seed = 1L
#' )
#' }
#'
#' # Example 3: External design (control only), posterior probability, NIW prior
#' \dontrun{
#' Sigma  <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' Lambda <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' se_c   <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' getgamma2cont(
#'   nsim = 1000L, prob = 'posterior', design = 'external',
#'   prior = 'N-Inv-Wishart',
#'   GoRegions = 1L, NoGoRegions = 9L,
#'   mu_t_go = c(-5.0, 0.0), Sigma_t_go = Sigma,
#'   mu_c_go = c(-10.0, -1.0), Sigma_c_go = Sigma,
#'   mu_t_nogo = c(5.0, 1.0), Sigma_t_nogo = Sigma,
#'   mu_c_nogo = c(-10.0, -1.0), Sigma_c_nogo = Sigma,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 30L, n_c = 30L,
#'   theta_TV1 = 10.0, theta_MAV1 = 5.0,
#'   theta_TV2 = 2.0,  theta_MAV2 = 1.0,
#'   theta_NULL1 = NULL, theta_NULL2 = NULL,
#'   m_t = NULL, m_c = NULL,
#'   kappa0_t = 0.1, nu0_t = 4.0, mu0_t = c(0.0, 1.0),  Lambda0_t = Lambda,
#'   kappa0_c = 0.1, nu0_c = 4.0, mu0_c = c(-10.0, -1.0), Lambda0_c = Lambda,
#'   r = NULL,
#'   ne_t = NULL, ne_c = 10L, alpha0e_t = NULL, alpha0e_c = 0.5,
#'   bar_ye_t = NULL, bar_ye_c = c(-10.0, -1.0), se_t = NULL, se_c = se_c,
#'   nMC = 500L, CalcMethod = 'MC',
#'   gamma_go_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma_nogo_grid = seq(0.01, 0.99, by = 0.01),
#'   seed = 1L
#' )
#' }
#'
#' # Example 4: Controlled design, predictive probability, vague prior
#' \dontrun{
#' Sigma_null <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' Sigma_alt  <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' getgamma2cont(
#'   nsim = 1000L, prob = 'predictive', design = 'controlled',
#'   prior = 'vague',
#'   GoRegions = 1L, NoGoRegions = 4L,
#'   mu_t_go = c(-5.0, 0.0), Sigma_t_go = Sigma_null,
#'   mu_c_go = c(-10.0, -1.0), Sigma_c_go = Sigma_null,
#'   mu_t_nogo = c(5.0, 1.0), Sigma_t_nogo = Sigma_alt,
#'   mu_c_nogo = c(-10.0, -1.0), Sigma_c_nogo = Sigma_alt,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 30L, n_c = 30L,
#'   theta_TV1 = NULL, theta_MAV1 = NULL,
#'   theta_TV2 = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 5.0, theta_NULL2 = 1.0,
#'   m_t = 100L, m_c = 100L,
#'   kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
#'   kappa0_c = NULL, nu0_c = NULL, mu0_c = NULL, Lambda0_c = NULL,
#'   r = NULL,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   nMC = 500L, CalcMethod = 'MC',
#'   gamma_go_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma_nogo_grid = seq(0.01, 0.99, by = 0.01),
#'   seed = 1L
#' )
#' }
#'
#' # Example 5: Uncontrolled design, predictive probability, vague prior
#' \dontrun{
#' Sigma_null <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' Sigma_alt  <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' getgamma2cont(
#'   nsim = 1000L, prob = 'predictive', design = 'uncontrolled',
#'   prior = 'vague',
#'   GoRegions = 1L, NoGoRegions = 4L,
#'   mu_t_go = c(-5.0, 0.0), Sigma_t_go = Sigma_null,
#'   mu_c_go = NULL, Sigma_c_go = NULL,
#'   mu_t_nogo = c(5.0, 1.0), Sigma_t_nogo = Sigma_alt,
#'   mu_c_nogo = NULL, Sigma_c_nogo = NULL,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 30L, n_c = NULL,
#'   theta_TV1 = NULL, theta_MAV1 = NULL,
#'   theta_TV2 = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 5.0, theta_NULL2 = 1.0,
#'   m_t = 100L, m_c = 100L,
#'   kappa0_t = NULL, nu0_t = NULL, mu0_t = NULL, Lambda0_t = NULL,
#'   kappa0_c = NULL, nu0_c = NULL, mu0_c = c(-10.0, -1.0), Lambda0_c = NULL,
#'   r = 1.0,
#'   ne_t = NULL, ne_c = NULL, alpha0e_t = NULL, alpha0e_c = NULL,
#'   bar_ye_t = NULL, bar_ye_c = NULL, se_t = NULL, se_c = NULL,
#'   nMC = 500L, CalcMethod = 'MC',
#'   gamma_go_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma_nogo_grid = seq(0.01, 0.99, by = 0.01),
#'   seed = 1L
#' )
#' }
#'
#' # Example 6: External design (control only), predictive probability, NIW prior
#' \dontrun{
#' Sigma  <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' Lambda <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' se_c   <- matrix(c(6400.0, 15.0, 15.0, 36.0), 2, 2)
#' getgamma2cont(
#'   nsim = 1000L, prob = 'predictive', design = 'external',
#'   prior = 'N-Inv-Wishart',
#'   GoRegions = 1L, NoGoRegions = 4L,
#'   mu_t_go = c(-5.0, 0.0), Sigma_t_go = Sigma,
#'   mu_c_go = c(-10.0, -1.0), Sigma_c_go = Sigma,
#'   mu_t_nogo = c(5.0, 1.0), Sigma_t_nogo = Sigma,
#'   mu_c_nogo = c(-10.0, -1.0), Sigma_c_nogo = Sigma,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 30L, n_c = 30L,
#'   theta_TV1 = NULL, theta_MAV1 = NULL,
#'   theta_TV2 = NULL, theta_MAV2 = NULL,
#'   theta_NULL1 = 5.0, theta_NULL2 = 1.0,
#'   m_t = 100L, m_c = 100L,
#'   kappa0_t = 0.1, nu0_t = 4.0, mu0_t = c(0.0, 1.0),  Lambda0_t = Lambda,
#'   kappa0_c = 0.1, nu0_c = 4.0, mu0_c = c(-10.0, -1.0), Lambda0_c = Lambda,
#'   r = NULL,
#'   ne_t = NULL, ne_c = 10L, alpha0e_t = NULL, alpha0e_c = 0.5,
#'   bar_ye_t = NULL, bar_ye_c = c(-10.0, -1.0), se_t = NULL, se_c = se_c,
#'   nMC = 500L, CalcMethod = 'MC',
#'   gamma_go_grid = seq(0.01, 0.99, by = 0.01),
#'   gamma_nogo_grid = seq(0.01, 0.99, by = 0.01),
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
                          mu_t_go, Sigma_t_go,
                          mu_c_go    = NULL, Sigma_c_go  = NULL,
                          mu_t_nogo, Sigma_t_nogo,
                          mu_c_nogo  = NULL, Sigma_c_nogo = NULL,
                          target_go, target_nogo,
                          n_t, n_c = NULL,
                          theta_TV1   = NULL, theta_MAV1  = NULL,
                          theta_TV2   = NULL, theta_MAV2  = NULL,
                          theta_NULL1 = NULL, theta_NULL2 = NULL,
                          m_t = NULL, m_c = NULL,
                          kappa0_t  = NULL, nu0_t     = NULL,
                          mu0_t     = NULL, Lambda0_t = NULL,
                          kappa0_c  = NULL, nu0_c     = NULL,
                          mu0_c     = NULL, Lambda0_c = NULL,
                          r        = NULL,
                          ne_t = NULL, ne_c = NULL,
                          alpha0e_t = NULL, alpha0e_c = NULL,
                          bar_ye_t  = NULL, bar_ye_c  = NULL,
                          se_t      = NULL, se_c      = NULL,
                          nMC    = NULL,
                          CalcMethod = 'MC',
                          gamma_go_grid = seq(0.01, 0.99, by = 0.01),
                          gamma_nogo_grid = seq(0.01, 0.99, by = 0.01),
                          seed) {

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

  # Validate Go-calibration scenario parameters
  if (!is.numeric(mu_t_go) || length(mu_t_go) != 2L || any(is.na(mu_t_go))) {
    stop("'mu_t_go' must be a length-2 numeric vector")
  }
  if (!is.matrix(Sigma_t_go) || !all(dim(Sigma_t_go) == c(2L, 2L)) ||
      !is.numeric(Sigma_t_go) || any(is.na(Sigma_t_go))) {
    stop("'Sigma_t_go' must be a 2x2 numeric matrix")
  }

  # Validate NoGo-calibration scenario parameters
  if (!is.numeric(mu_t_nogo) || length(mu_t_nogo) != 2L || any(is.na(mu_t_nogo))) {
    stop("'mu_t_nogo' must be a length-2 numeric vector")
  }
  if (!is.matrix(Sigma_t_nogo) || !all(dim(Sigma_t_nogo) == c(2L, 2L)) ||
      !is.numeric(Sigma_t_nogo) || any(is.na(Sigma_t_nogo))) {
    stop("'Sigma_t_nogo' must be a 2x2 numeric matrix")
  }

  if (design != 'uncontrolled') {
    if (is.null(mu_c_go) || !is.numeric(mu_c_go) || length(mu_c_go) != 2L ||
        any(is.na(mu_c_go))) {
      stop("'mu_c_go' must be a length-2 numeric vector for controlled or external design")
    }
    if (is.null(Sigma_c_go) || !is.matrix(Sigma_c_go) ||
        !all(dim(Sigma_c_go) == c(2L, 2L)) || !is.numeric(Sigma_c_go) ||
        any(is.na(Sigma_c_go))) {
      stop("'Sigma_c_go' must be a 2x2 numeric matrix for controlled or external design")
    }
    if (is.null(mu_c_nogo) || !is.numeric(mu_c_nogo) || length(mu_c_nogo) != 2L ||
        any(is.na(mu_c_nogo))) {
      stop("'mu_c_nogo' must be a length-2 numeric vector for controlled or external design")
    }
    if (is.null(Sigma_c_nogo) || !is.matrix(Sigma_c_nogo) ||
        !all(dim(Sigma_c_nogo) == c(2L, 2L)) || !is.numeric(Sigma_c_nogo) ||
        any(is.na(Sigma_c_nogo))) {
      stop("'Sigma_c_nogo' must be a 2x2 numeric matrix for controlled or external design")
    }
  }

  for (nm in c("target_go", "target_nogo")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val <= 0 || val >= 1) {
      stop(paste0("'", nm, "' must be a single numeric value in (0, 1)"))
    }
  }

  if (!is.numeric(n_t) || length(n_t) != 1L || is.na(n_t) ||
      n_t != floor(n_t) || n_t < 1L) {
    stop("'n_t' must be a single positive integer")
  }
  if (design != 'uncontrolled') {
    if (is.null(n_c) || !is.numeric(n_c) || length(n_c) != 1L || is.na(n_c) ||
        n_c != floor(n_c) || n_c < 1L) {
      stop("'n_c' must be a single positive integer for controlled or external design")
    }
  }

  if (prob == 'posterior') {
    for (nm in c("theta_TV1", "theta_MAV1", "theta_TV2", "theta_MAV2")) {
      if (is.null(get(nm))) {
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

  if (!is.character(CalcMethod) || length(CalcMethod) != 1L ||
      !CalcMethod %in% c('MC', 'MM')) {
    stop("'CalcMethod' must be either 'MC' or 'MM'")
  }
  if (CalcMethod == 'MC') {
    if (is.null(nMC) || !is.numeric(nMC) || length(nMC) != 1L ||
        is.na(nMC) || nMC != floor(nMC) || nMC < 1L) {
      stop("'nMC' must be a single positive integer when CalcMethod = 'MC'")
    }
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
  ng_go           <- length(gamma_go_grid)
  ng_nogo         <- length(gamma_nogo_grid)

  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("'seed' must be a single numeric value")
  }

  # ---------------------------------------------------------------------------
  # Internal helper: simulate nsim datasets and return PrGo_vec, PrNoGo_vec
  # ---------------------------------------------------------------------------
  .simulate_g <- function(mu_t_s, Sigma_t_s, mu_c_s, Sigma_c_s, seed_val) {
    set.seed(seed_val)

    n_t_int <- as.integer(n_t)
    R_Sigma_t   <- chol(Sigma_t_s)
    Z_t_raw     <- matrix(rnorm(nsim * n_t_int * 2L),
                          nrow = nsim * n_t_int, ncol = 2L) %*% R_Sigma_t
    block_t     <- rep(seq_len(nsim), each = n_t_int)
    Z_t_colsums <- apply(Z_t_raw, 2L, function(col) tapply(col, block_t, sum))
    Z_t_colmeans_rep <- Z_t_colsums[block_t, ] / n_t_int
    Z_t_centered     <- Z_t_raw - Z_t_colmeans_rep
    S_t_11 <- tapply(Z_t_centered[, 1L] ^ 2,                  block_t, sum)
    S_t_12 <- tapply(Z_t_centered[, 1L] * Z_t_centered[, 2L], block_t, sum)
    S_t_22 <- tapply(Z_t_centered[, 2L] ^ 2,                  block_t, sum)
    ybar_t <- sweep(Z_t_colsums / n_t_int, 2L, mu_t_s, '+')
    S_t_list <- vector('list', nsim)
    for (i in seq_len(nsim)) {
      S_t_list[[i]] <- matrix(c(S_t_11[i], S_t_12[i], S_t_12[i], S_t_22[i]),
                              nrow = 2L, ncol = 2L)
    }

    if (design %in% c('controlled', 'external')) {
      n_c_int <- as.integer(n_c)
      R_Sigma_c   <- chol(Sigma_c_s)
      Z_c_raw     <- matrix(rnorm(nsim * n_c_int * 2L),
                            nrow = nsim * n_c_int, ncol = 2L) %*% R_Sigma_c
      block_c     <- rep(seq_len(nsim), each = n_c_int)
      Z_c_colsums <- apply(Z_c_raw, 2L, function(col) tapply(col, block_c, sum))
      Z_c_colmeans_rep <- Z_c_colsums[block_c, ] / n_c_int
      Z_c_centered     <- Z_c_raw - Z_c_colmeans_rep
      S_c_11 <- tapply(Z_c_centered[, 1L] ^ 2,                  block_c, sum)
      S_c_12 <- tapply(Z_c_centered[, 1L] * Z_c_centered[, 2L], block_c, sum)
      S_c_22 <- tapply(Z_c_centered[, 2L] ^ 2,                  block_c, sum)
      ybar_c <- sweep(Z_c_colsums / n_c_int, 2L, mu_c_s, '+')
      S_c_list <- vector('list', nsim)
      for (i in seq_len(nsim)) {
        S_c_list[[i]] <- matrix(c(S_c_11[i], S_c_12[i], S_c_12[i], S_c_22[i]),
                                nrow = 2L, ncol = 2L)
      }
    } else {
      ybar_c   <- NULL
      S_c_list <- NULL
      n_c_int  <- NULL
    }

    # Vectorised call to pbayespostpred2cont: returns nsim x n_regions matrix
    Pr_R_mat <- pbayespostpred2cont(
      prob        = prob,
      design      = design,
      prior       = prior,
      theta_TV1   = theta_TV1,   theta_MAV1  = theta_MAV1,
      theta_TV2   = theta_TV2,   theta_MAV2  = theta_MAV2,
      theta_NULL1 = theta_NULL1, theta_NULL2 = theta_NULL2,
      n_t = n_t_int, n_c = n_c_int,
      ybar_t = ybar_t, S_t = S_t_list,
      ybar_c = ybar_c, S_c = S_c_list,
      m_t = m_t, m_c = m_c,
      kappa0_t  = kappa0_t,  nu0_t    = nu0_t,
      mu0_t     = mu0_t,     Lambda0_t = Lambda0_t,
      kappa0_c  = kappa0_c,  nu0_c    = nu0_c,
      mu0_c     = mu0_c,     Lambda0_c = Lambda0_c,
      r         = r,
      ne_t = ne_t, ne_c = ne_c,
      alpha0e_t = alpha0e_t, alpha0e_c = alpha0e_c,
      bar_ye_t  = bar_ye_t,  bar_ye_c  = bar_ye_c,
      se_t      = se_t,      se_c      = se_c,
      nMC    = nMC,
      CalcMethod = CalcMethod
    )

    # Sum region probabilities over GoRegions and NoGoRegions
    PrGo_vec   <- rowSums(Pr_R_mat[, GoRegions,   drop = FALSE])
    PrNoGo_vec <- rowSums(Pr_R_mat[, NoGoRegions, drop = FALSE])

    list(PrGo_vec = PrGo_vec, PrNoGo_vec = PrNoGo_vec)
  }

  # ---------------------------------------------------------------------------
  # Stage 1: Simulate nsim datasets for each calibration scenario independently
  # ---------------------------------------------------------------------------
  # Go-calibration scenario (typically Null): use seed
  sim_go   <- .simulate_g(mu_t_go,   Sigma_t_go,   mu_c_go,   Sigma_c_go,
                          seed_val = seed)
  # NoGo-calibration scenario (typically Alternative): use seed + 1
  sim_nogo <- .simulate_g(mu_t_nogo, Sigma_t_nogo, mu_c_nogo, Sigma_c_nogo,
                          seed_val = seed + 1L)

  PrGo_vec   <- sim_go$PrGo_vec
  PrNoGo_vec <- sim_nogo$PrNoGo_vec

  # ---------------------------------------------------------------------------
  # Stage 2: Sweep (gamma_go_grid x gamma_nogo_grid)
  #
  # For each (g_go, g_nogo), compute under respective calibration scenarios:
  #   Pr(Go)   = mean( I(PrGo_vec >= g_go   AND PrNoGo_vec_go   <  g_nogo) )
  #   Pr(NoGo) = mean( I(PrNoGo_vec >= g_nogo AND PrGo_vec_nogo  <  g_go) )
  # ---------------------------------------------------------------------------
  PrGo_vec_go     <- sim_go$PrGo_vec
  PrNoGo_vec_go   <- sim_go$PrNoGo_vec
  PrGo_vec_nogo   <- sim_nogo$PrGo_vec
  PrNoGo_vec_nogo <- sim_nogo$PrNoGo_vec

  PrGo_grid   <- matrix(NA_real_, nrow = ng_go, ncol = ng_nogo)
  PrNoGo_grid <- matrix(NA_real_, nrow = ng_go, ncol = ng_nogo)

  for (k_go in seq_len(ng_go)) {
    g_go <- gamma_go_grid[k_go]
    for (k_nogo in seq_len(ng_nogo)) {
      g_nogo <- gamma_nogo_grid[k_nogo]

      # Pr(Go) evaluated under Go-calibration scenario
      go_mask   <- (PrGo_vec_go   >= g_go) & (PrNoGo_vec_go   <  g_nogo)
      # Pr(NoGo) evaluated under NoGo-calibration scenario
      nogo_mask <- (PrNoGo_vec_nogo >= g_nogo) & (PrGo_vec_nogo <  g_go)

      PrGo_grid[k_go, k_nogo]   <- mean(go_mask)
      PrNoGo_grid[k_go, k_nogo] <- mean(nogo_mask)
    }
  }

  # ---------------------------------------------------------------------------
  # Stage 3: Select optimal (gamma_go, gamma_nogo)
  #
  # gamma_go  : smallest gamma_go   s.t. max_{gamma_nogo} Pr(Go)   < target_go
  # gamma_nogo: smallest gamma_nogo s.t. max_{gamma_go}   Pr(NoGo) < target_nogo
  # Both worst-case curves are non-increasing in their respective gamma.
  # ---------------------------------------------------------------------------
  max_PrGo_per_g1   <- apply(PrGo_grid,   1L, max, na.rm = TRUE)
  max_PrNoGo_per_g2 <- apply(PrNoGo_grid, 2L, max, na.rm = TRUE)

  go_mask_opt   <- max_PrGo_per_g1   < target_go
  nogo_mask_opt <- max_PrNoGo_per_g2 < target_nogo

  idx_go <- which(go_mask_opt)
  if (length(idx_go) == 0L) {
    gamma_go <- NA_real_
    PrGo_opt <- NA_real_
  } else {
    opt1     <- min(idx_go)
    gamma_go <- gamma_go_grid[opt1]
    PrGo_opt <- max_PrGo_per_g1[opt1]
  }

  idx_nogo <- which(nogo_mask_opt)
  if (length(idx_nogo) == 0L) {
    gamma_nogo <- NA_real_
    PrNoGo_opt <- NA_real_
  } else {
    opt2       <- min(idx_nogo)
    gamma_nogo <- gamma_nogo_grid[opt2]
    PrNoGo_opt <- max_PrNoGo_per_g2[opt2]
  }

  # ---------------------------------------------------------------------------
  # Build and return result
  # ---------------------------------------------------------------------------
  result <- list(
    gamma_go     = gamma_go,
    gamma_nogo   = gamma_nogo,
    PrGo_opt     = PrGo_opt,
    PrNoGo_opt   = PrNoGo_opt,
    grid_results = list(
      gamma_go_grid   = gamma_go_grid,
      gamma_nogo_grid = gamma_nogo_grid,
      PrGo_grid       = PrGo_grid,
      PrNoGo_grid     = PrNoGo_grid
    )
  )

  class(result) <- 'getgamma2cont'
  return(result)
}
