#' Find Optimal Go/NoGo Thresholds for a Single Binary Endpoint
#'
#' Computes the optimal Go threshold \eqn{\gamma_{\mathrm{go}}} and NoGo
#' threshold \eqn{\gamma_{\mathrm{nogo}}} for a single binary endpoint by
#' searching over a grid of candidate values.  The two thresholds are
#' calibrated independently under separate scenarios:
#' \itemize{
#'   \item \eqn{\gamma_{\mathrm{go}}} is the \strong{smallest} value in
#'         \code{gamma_grid} such that the marginal Go probability
#'         \eqn{\Pr(g_{\mathrm{Go}} \ge \gamma_{\mathrm{go}})} is strictly
#'         less than \code{target_go} under the Go-calibration scenario
#'         (\code{pi_t_go}, \code{pi_c_go}); typically the Null scenario.
#'   \item \eqn{\gamma_{\mathrm{nogo}}} is the \strong{smallest} value in
#'         \code{gamma_grid} such that the marginal NoGo probability
#'         \eqn{\Pr(g_{\mathrm{NoGo}} \ge \gamma_{\mathrm{nogo}})} is strictly
#'         less than \code{target_nogo} under the NoGo-calibration
#'         scenario (\code{pi_t_nogo}, \code{pi_c_nogo}); typically the
#'         Alternative scenario.
#' }
#' Here \eqn{g_{\mathrm{Go}} = P(\theta > \theta_{\mathrm{TV}} \mid y_t, y_c)}
#' and \eqn{g_{\mathrm{NoGo}} = P(\theta \le \theta_{\mathrm{MAV}} \mid y_t, y_c)}
#' for \code{prob = 'posterior'}, consistent with the decision rule in
#' \code{\link{pbayesdecisionprob1bin}}.
#'
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or
#'        \code{'external'}.
#' @param theta_TV A numeric scalar in \code{(-1, 1)} giving the Target Value
#'        (TV) threshold for the treatment effect.  Required when
#'        \code{prob = 'posterior'}; set to \code{NULL} otherwise.
#' @param theta_MAV A numeric scalar in \code{(-1, 1)} giving the Minimum
#'        Acceptable Value (MAV) threshold.  Must satisfy
#'        \code{theta_TV > theta_MAV}.  Required when
#'        \code{prob = 'posterior'}; set to \code{NULL} otherwise.
#' @param theta_NULL A numeric scalar in \code{(-1, 1)} giving the null
#'        hypothesis threshold used for the predictive probability.  Required
#'        when \code{prob = 'predictive'}; set to \code{NULL} otherwise.
#' @param pi_t_go A numeric scalar in \code{(0, 1)} giving the true treatment
#'        response rate under the Go-calibration scenario (typically Null).
#' @param pi_c_go A numeric scalar in \code{(0, 1)} giving the true control
#'        response rate under the Go-calibration scenario.  Set to \code{NULL}
#'        for \code{design = 'uncontrolled'}.
#' @param pi_t_nogo A numeric scalar in \code{(0, 1)} giving the true
#'        treatment response rate under the NoGo-calibration scenario
#'        (typically Alternative).
#' @param pi_c_nogo A numeric scalar in \code{(0, 1)} giving the true control
#'        response rate under the NoGo-calibration scenario.  Set to
#'        \code{NULL} for \code{design = 'uncontrolled'}.
#' @param target_go A numeric scalar in \code{(0, 1)} giving the upper bound
#'        on the marginal Go probability under the Go-calibration scenario.
#'        The optimal \eqn{\gamma_{\mathrm{go}}} is the smallest grid value
#'        satisfying \eqn{\Pr(\mathrm{Go}) < \code{target\_go}}.
#' @param target_nogo A numeric scalar in \code{(0, 1)} giving the upper bound
#'        on the marginal NoGo probability under the NoGo-calibration scenario.
#'        The optimal \eqn{\gamma_{\mathrm{nogo}}} is the largest grid value
#'        satisfying \eqn{\Pr(\mathrm{NoGo}) < \code{target\_nogo}}.
#' @param n_t A positive integer giving the number of patients in the
#'        treatment group in the PoC trial.
#' @param n_c A positive integer giving the number of patients in the
#'        control group in the PoC trial.
#' @param a_t A positive numeric scalar giving the first shape parameter of
#'        the Beta prior for the treatment group.
#' @param a_c A positive numeric scalar giving the first shape parameter of
#'        the Beta prior for the control group.
#' @param b_t A positive numeric scalar giving the second shape parameter of
#'        the Beta prior for the treatment group.
#' @param b_c A positive numeric scalar giving the second shape parameter of
#'        the Beta prior for the control group.
#' @param z A non-negative integer giving the hypothetical number of
#'        responders in the control group.  Required when
#'        \code{design = 'uncontrolled'}; set to \code{NULL} otherwise.
#' @param m_t A positive integer giving the future sample size for the
#'        treatment group.  Required when \code{prob = 'predictive'}; set to
#'        \code{NULL} otherwise.
#' @param m_c A positive integer giving the future sample size for the
#'        control group.  Required when \code{prob = 'predictive'}; set to
#'        \code{NULL} otherwise.
#' @param ne_t A positive integer giving the number of patients in the
#'        treatment group of the external data set.  Required when
#'        \code{design = 'external'}; set to \code{NULL} otherwise.
#' @param ne_c A positive integer giving the number of patients in the
#'        control group of the external data set.  Required when
#'        \code{design = 'external'}; set to \code{NULL} otherwise.
#' @param ye_t A non-negative integer giving the number of responders in the
#'        treatment group of the external data set.  Required when
#'        \code{design = 'external'}; set to \code{NULL} otherwise.
#' @param ye_c A non-negative integer giving the number of responders in the
#'        control group of the external data set.  Required when
#'        \code{design = 'external'}; set to \code{NULL} otherwise.
#' @param alpha0e_t A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for the treatment group.  Required when
#'        \code{design = 'external'}; set to \code{NULL} otherwise.
#' @param alpha0e_c A numeric scalar in \code{(0, 1]} giving the power prior
#'        weight for the control group.  Required when
#'        \code{design = 'external'}; set to \code{NULL} otherwise.
#' @param gamma_grid A numeric vector of candidate threshold values in
#'        \code{(0, 1)} to search over.  Defaults to
#'        \code{seq(0.01, 0.99, by = 0.01)}.
#'
#' @return A list of class \code{getgamma1bin} with the following elements:
#' \describe{
#'   \item{gamma_go}{Optimal Go threshold: the smallest value in
#'         \code{gamma_grid} for which \eqn{\Pr(\mathrm{Go}) < \code{target\_go}}
#'         under the Go-calibration scenario.
#'         \code{NA} if no such value exists.}
#'   \item{gamma_nogo}{Optimal NoGo threshold: the smallest value in
#'         \code{gamma_grid} for which \eqn{\Pr(\mathrm{NoGo}) < \code{target\_nogo}}
#'         under the NoGo-calibration scenario.
#'         \code{NA} if no such value exists.}
#'   \item{PrGo_opt}{Marginal \eqn{\Pr(g_{\mathrm{Go}} \ge \gamma_{\mathrm{go}})}
#'         at the optimal \eqn{\gamma_{\mathrm{go}}} under the Go-calibration
#'         scenario.  \code{NA} if \code{gamma_go} is \code{NA}.}
#'   \item{PrNoGo_opt}{Marginal \eqn{\Pr(g_{\mathrm{NoGo}} \ge \gamma_{\mathrm{nogo}})}
#'         at the optimal \eqn{\gamma_{\mathrm{nogo}}} under the NoGo-calibration
#'         scenario.  \code{NA} if \code{gamma_nogo} is \code{NA}.}
#'   \item{grid_results}{A data frame with columns \code{gamma_grid},
#'         \code{PrGo_grid} (marginal Go probability under the Go-calibration
#'         scenario), and \code{PrNoGo_grid} (marginal NoGo probability under
#'         the NoGo-calibration scenario).}
#' }
#'
#' @details
#' The function uses a two-stage precompute-then-sweep strategy:
#' \enumerate{
#'   \item \strong{Precomputation}: All possible outcome pairs
#'         \eqn{(y_t, y_c)} are enumerated.  For each pair,
#'         \code{\link{pbayespostpred1bin}} computes
#'         \eqn{g_{\mathrm{Go}}} (\code{lower.tail = FALSE} at
#'         \code{theta_TV}) and \eqn{g_{\mathrm{NoGo}}} (\code{lower.tail
#'         = TRUE} at \code{theta_MAV}).  This step is independent of
#'         \eqn{\gamma}.
#'   \item \strong{Gamma sweep}: Marginal probabilities are computed as
#'         weighted sums of binary indicators over the grid:
#'         \eqn{\Pr(\mathrm{Go})} uses \code{w_go} (weights under
#'         \code{pi_t_go}, \code{pi_c_go}) and the indicator
#'         \eqn{g_{\mathrm{Go}} \ge \gamma}; \eqn{\Pr(\mathrm{NoGo})}
#'         uses \code{w_nogo} (weights under \code{pi_t_nogo},
#'         \code{pi_c_nogo}) and the indicator
#'         \eqn{g_{\mathrm{NoGo}} \ge \gamma}.
#' }
#' Both \eqn{\Pr(\mathrm{Go})} and \eqn{\Pr(\mathrm{NoGo})} are
#' monotone non-increasing functions of \eqn{\gamma}.  The optimal
#' \eqn{\gamma_{\mathrm{go}}} is the \emph{smallest} grid value
#' crossing below \code{target_go}.  The optimal
#' \eqn{\gamma_{\mathrm{nogo}}} is also the \emph{smallest} grid value
#' crossing below \code{target_nogo}: a smaller \eqn{\gamma_{\mathrm{nogo}}}
#' makes NoGo harder to trigger (more permissive), so this is the least
#' restrictive threshold that still controls the false NoGo rate.
#'
#' @examples
#' # Example 1: Controlled design, posterior probability
#' # gamma_go  : smallest gamma s.t. Pr(Go)   < 0.05 under Null (pi_t = pi_c = 0.15)
#' # gamma_nogo: largest  gamma s.t. Pr(NoGo) < 0.20 under Alt  (pi_t = 0.35, pi_c = 0.15)
#' getgamma1bin(
#'   prob = 'posterior', design = 'controlled',
#'   theta_TV = 0.20, theta_MAV = 0.05, theta_NULL = NULL,
#'   pi_t_go = 0.15, pi_c_go = 0.15,
#'   pi_t_nogo = 0.35, pi_c_nogo = 0.15,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 12L, n_c = 12L,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = NULL, m_t = NULL, m_c = NULL,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL
#' )
#'
#' # Example 2: Uncontrolled design, posterior probability
#' getgamma1bin(
#'   prob = 'posterior', design = 'uncontrolled',
#'   theta_TV = 0.20, theta_MAV = 0.05, theta_NULL = NULL,
#'   pi_t_go = 0.15, pi_c_go = NULL,
#'   pi_t_nogo = 0.35, pi_c_nogo = NULL,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 12L, n_c = 12L,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = 3L, m_t = NULL, m_c = NULL,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL
#' )
#'
#' # Example 3: External design, posterior probability
#' getgamma1bin(
#'   prob = 'posterior', design = 'external',
#'   theta_TV = 0.20, theta_MAV = 0.05, theta_NULL = NULL,
#'   pi_t_go = 0.15, pi_c_go = 0.15,
#'   pi_t_nogo = 0.35, pi_c_nogo = 0.15,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 12L, n_c = 12L,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = NULL, m_t = NULL, m_c = NULL,
#'   ne_t = 15L, ne_c = 15L, ye_t = 6L, ye_c = 4L,
#'   alpha0e_t = 0.5, alpha0e_c = 0.5
#' )
#'
#' # Example 4: Controlled design, predictive probability
#' getgamma1bin(
#'   prob = 'predictive', design = 'controlled',
#'   theta_TV = NULL, theta_MAV = NULL, theta_NULL = 0.10,
#'   pi_t_go = 0.15, pi_c_go = 0.15,
#'   pi_t_nogo = 0.35, pi_c_nogo = 0.15,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 12L, n_c = 12L,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = NULL, m_t = 30L, m_c = 30L,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL
#' )
#'
#' # Example 5: Uncontrolled design, predictive probability
#' getgamma1bin(
#'   prob = 'predictive', design = 'uncontrolled',
#'   theta_TV = NULL, theta_MAV = NULL, theta_NULL = 0.10,
#'   pi_t_go = 0.15, pi_c_go = NULL,
#'   pi_t_nogo = 0.35, pi_c_nogo = NULL,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 12L, n_c = 12L,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = 3L, m_t = 30L, m_c = 30L,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL,
#'   alpha0e_t = NULL, alpha0e_c = NULL
#' )
#'
#' # Example 6: External design, predictive probability
#' getgamma1bin(
#'   prob = 'predictive', design = 'external',
#'   theta_TV = NULL, theta_MAV = NULL, theta_NULL = 0.10,
#'   pi_t_go = 0.15, pi_c_go = 0.15,
#'   pi_t_nogo = 0.35, pi_c_nogo = 0.15,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n_t = 12L, n_c = 12L,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = NULL, m_t = 30L, m_c = 30L,
#'   ne_t = 15L, ne_c = 15L, ye_t = 6L, ye_c = 4L,
#'   alpha0e_t = 0.5, alpha0e_c = 0.5
#' )
#'
#' @importFrom stats dbinom
#' @export
getgamma1bin <- function(prob = 'posterior', design = 'controlled',
                         theta_TV = NULL, theta_MAV = NULL, theta_NULL = NULL,
                         pi_t_go, pi_c_go = NULL,
                         pi_t_nogo, pi_c_nogo = NULL,
                         target_go, target_nogo,
                         n_t, n_c,
                         a_t, a_c, b_t, b_c,
                         z = NULL, m_t = NULL, m_c = NULL,
                         ne_t = NULL, ne_c = NULL,
                         ye_t = NULL, ye_c = NULL,
                         alpha0e_t = NULL, alpha0e_c = NULL,
                         gamma_grid = seq(0.01, 0.99, by = 0.01)) {

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

  if (prob == 'posterior') {
    if (is.null(theta_TV) || is.null(theta_MAV)) {
      stop("'theta_TV' and 'theta_MAV' must be non-NULL when prob = 'posterior'")
    }
    if (!is.numeric(theta_TV) || length(theta_TV) != 1L || is.na(theta_TV) ||
        theta_TV <= -1 || theta_TV >= 1) {
      stop("'theta_TV' must be a single numeric value in (-1, 1)")
    }
    if (!is.numeric(theta_MAV) || length(theta_MAV) != 1L || is.na(theta_MAV) ||
        theta_MAV <= -1 || theta_MAV >= 1) {
      stop("'theta_MAV' must be a single numeric value in (-1, 1)")
    }
    if (theta_TV <= theta_MAV) {
      stop("'theta_TV' must be strictly greater than 'theta_MAV'")
    }
  } else {
    if (is.null(theta_NULL)) {
      stop("'theta_NULL' must be non-NULL when prob = 'predictive'")
    }
    if (!is.numeric(theta_NULL) || length(theta_NULL) != 1L ||
        is.na(theta_NULL) || theta_NULL <= -1 || theta_NULL >= 1) {
      stop("'theta_NULL' must be a single numeric value in (-1, 1)")
    }
    if (is.null(m_t) || is.null(m_c)) {
      stop("'m_t' and 'm_c' must be non-NULL when prob = 'predictive'")
    }
    for (nm in c("m_t", "m_c")) {
      val <- get(nm)
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 1L) {
        stop(paste0("'", nm, "' must be a single positive integer"))
      }
    }
  }

  if (!is.numeric(pi_t_go) || length(pi_t_go) != 1L || is.na(pi_t_go) ||
      pi_t_go <= 0 || pi_t_go >= 1) {
    stop("'pi_t_go' must be a single numeric value in (0, 1)")
  }
  if (!is.numeric(pi_t_nogo) || length(pi_t_nogo) != 1L || is.na(pi_t_nogo) ||
      pi_t_nogo <= 0 || pi_t_nogo >= 1) {
    stop("'pi_t_nogo' must be a single numeric value in (0, 1)")
  }

  if (design != 'uncontrolled') {
    if (is.null(pi_c_go) || !is.numeric(pi_c_go) || length(pi_c_go) != 1L ||
        is.na(pi_c_go) || pi_c_go <= 0 || pi_c_go >= 1) {
      stop("'pi_c_go' must be a single numeric value in (0, 1) for controlled or external design")
    }
    if (is.null(pi_c_nogo) || !is.numeric(pi_c_nogo) || length(pi_c_nogo) != 1L ||
        is.na(pi_c_nogo) || pi_c_nogo <= 0 || pi_c_nogo >= 1) {
      stop("'pi_c_nogo' must be a single numeric value in (0, 1) for controlled or external design")
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

  for (nm in c("a_t", "a_c", "b_t", "b_c")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0) {
      stop(paste0("'", nm, "' must be a single positive numeric value"))
    }
  }

  if (design == 'uncontrolled') {
    if (is.null(z)) {
      stop("'z' must be non-NULL when design = 'uncontrolled'")
    }
    if (!is.numeric(z) || length(z) != 1L || is.na(z) ||
        z != floor(z) || z < 0L) {
      stop("'z' must be a single non-negative integer")
    }
  }

  if (design == 'external') {
    for (nm in c("ne_t", "ne_c", "ye_t", "ye_c")) {
      val <- get(nm)
      if (is.null(val)) {
        stop(paste0("'", nm, "' must be non-NULL when design = 'external'"))
      }
    }
    if (!is.numeric(ne_t) || length(ne_t) != 1L || is.na(ne_t) ||
        ne_t != floor(ne_t) || ne_t < 1L) {
      stop("'ne_t' must be a single positive integer")
    }
    if (!is.numeric(ne_c) || length(ne_c) != 1L || is.na(ne_c) ||
        ne_c != floor(ne_c) || ne_c < 1L) {
      stop("'ne_c' must be a single positive integer")
    }
    if (!is.numeric(ye_t) || length(ye_t) != 1L || is.na(ye_t) ||
        ye_t != floor(ye_t) || ye_t < 0L || ye_t > ne_t) {
      stop("'ye_t' must be a non-negative integer not exceeding 'ne_t'")
    }
    if (!is.numeric(ye_c) || length(ye_c) != 1L || is.na(ye_c) ||
        ye_c != floor(ye_c) || ye_c < 0L || ye_c > ne_c) {
      stop("'ye_c' must be a non-negative integer not exceeding 'ne_c'")
    }
    for (nm in c("alpha0e_t", "alpha0e_c")) {
      val <- get(nm)
      if (is.null(val) || !is.numeric(val) || length(val) != 1L ||
          is.na(val) || val <= 0 || val > 1) {
        stop(paste0("'", nm, "' must be a single numeric value in (0, 1]"))
      }
    }
  }

  if (!is.numeric(gamma_grid) || length(gamma_grid) < 1L ||
      any(is.na(gamma_grid)) || any(gamma_grid <= 0) || any(gamma_grid >= 1)) {
    stop("'gamma_grid' must be a numeric vector with all values in (0, 1)")
  }
  gamma_grid <- sort(unique(gamma_grid))

  # ---------------------------------------------------------------------------
  # Stage 1: Precompute g(y_t, y_c) for all outcome combinations (gamma-free)
  # ---------------------------------------------------------------------------
  theta0_go   <- if (prob == 'posterior') theta_TV  else theta_NULL
  theta0_nogo <- if (prob == 'posterior') theta_MAV else theta_NULL

  if (design == 'uncontrolled') {
    all_y_t <- 0L:n_t

    gPost_Go <- pbayespostpred1bin(
      prob = prob, design = design, theta0 = theta0_go,
      n_t = n_t, n_c = n_c, y_t = all_y_t, y_c = NULL,
      a_t = a_t, a_c = a_c, b_t = b_t, b_c = b_c,
      m_t = m_t, m_c = m_c, z = z,
      ne_t = ne_t, ne_c = ne_c, ye_t = ye_t, ye_c = ye_c,
      alpha0e_t = alpha0e_t, alpha0e_c = alpha0e_c,
      lower.tail = FALSE
    )

    gPost_NoGo <- pbayespostpred1bin(
      prob = prob, design = design, theta0 = theta0_nogo,
      n_t = n_t, n_c = n_c, y_t = all_y_t, y_c = NULL,
      a_t = a_t, a_c = a_c, b_t = b_t, b_c = b_c,
      m_t = m_t, m_c = m_c, z = z,
      ne_t = ne_t, ne_c = ne_c, ye_t = ye_t, ye_c = ye_c,
      alpha0e_t = alpha0e_t, alpha0e_c = alpha0e_c,
      lower.tail = TRUE
    )

    # Weights for Go calibration (typically Null scenario)
    w_go   <- dbinom(all_y_t, n_t, pi_t_go)
    # Weights for NoGo calibration (typically Alternative scenario)
    w_nogo <- dbinom(all_y_t, n_t, pi_t_nogo)

  } else {
    grid    <- expand.grid(y_t = 0L:n_t, y_c = 0L:n_c)
    all_y_t <- grid$y_t
    all_y_c <- grid$y_c

    gPost_Go <- pbayespostpred1bin(
      prob = prob, design = design, theta0 = theta0_go,
      n_t = n_t, n_c = n_c, y_t = all_y_t, y_c = all_y_c,
      a_t = a_t, a_c = a_c, b_t = b_t, b_c = b_c,
      m_t = m_t, m_c = m_c, z = NULL,
      ne_t = ne_t, ne_c = ne_c, ye_t = ye_t, ye_c = ye_c,
      alpha0e_t = alpha0e_t, alpha0e_c = alpha0e_c,
      lower.tail = FALSE
    )

    gPost_NoGo <- pbayespostpred1bin(
      prob = prob, design = design, theta0 = theta0_nogo,
      n_t = n_t, n_c = n_c, y_t = all_y_t, y_c = all_y_c,
      a_t = a_t, a_c = a_c, b_t = b_t, b_c = b_c,
      m_t = m_t, m_c = m_c, z = NULL,
      ne_t = ne_t, ne_c = ne_c, ye_t = ye_t, ye_c = ye_c,
      alpha0e_t = alpha0e_t, alpha0e_c = alpha0e_c,
      lower.tail = TRUE
    )

    # Weights for Go calibration (typically Null scenario)
    w_go   <- dbinom(all_y_t, n_t, pi_t_go)   * dbinom(all_y_c, n_c, pi_c_go)
    # Weights for NoGo calibration (typically Alternative scenario)
    w_nogo <- dbinom(all_y_t, n_t, pi_t_nogo) * dbinom(all_y_c, n_c, pi_c_nogo)
  }

  # ---------------------------------------------------------------------------
  # Stage 2: Sweep gamma_grid using vectorised outer comparison
  # ---------------------------------------------------------------------------
  # Marginal Go   criterion: g_Go   >= gamma  (lower.tail = FALSE -> P(theta > TV))
  # Marginal NoGo criterion: g_NoGo >= gamma  (lower.tail = TRUE  -> P(theta <= MAV))
  # Both Pr(Go) and Pr(NoGo) are monotone non-increasing in gamma.
  go_ind   <- outer(gPost_Go,   gamma_grid, ">=")  # n_outcomes x n_gamma
  nogo_ind <- outer(gPost_NoGo, gamma_grid, ">=")  # n_outcomes x n_gamma

  PrGo_grid   <- colSums(w_go   * go_ind)
  PrNoGo_grid <- colSums(w_nogo * nogo_ind)

  # ---------------------------------------------------------------------------
  # Stage 3: Select optimal gamma values
  # ---------------------------------------------------------------------------
  # gamma_go  : smallest gamma with Pr(Go)   < target_go
  #   -> Pr(Go)   is non-increasing; take the first (smallest) gamma crossing below target
  # gamma_nogo: smallest gamma with Pr(NoGo) < target_nogo
  #   -> Pr(NoGo) is non-increasing; a smaller gamma_nogo means NoGo is harder to trigger,
  #      so the smallest gamma satisfying the constraint is the most permissive threshold
  #      that still controls false NoGo rate below target_nogo
  go_mask   <- PrGo_grid   < target_go
  nogo_mask <- PrNoGo_grid < target_nogo

  idx_go <- which(go_mask)
  if (length(idx_go) == 0L) {
    gamma_go <- NA_real_
    PrGo_opt <- NA_real_
  } else {
    opt1     <- min(idx_go)
    gamma_go <- gamma_grid[opt1]
    PrGo_opt <- PrGo_grid[opt1]
  }

  idx_nogo <- which(nogo_mask)
  if (length(idx_nogo) == 0L) {
    gamma_nogo <- NA_real_
    PrNoGo_opt <- NA_real_
  } else {
    opt2       <- min(idx_nogo)      # smallest gamma with Pr(NoGo) < target_nogo
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
    grid_results = data.frame(
      gamma_grid  = gamma_grid,
      PrGo_grid   = PrGo_grid,
      PrNoGo_grid = PrNoGo_grid
    )
  )

  class(result) <- 'getgamma1bin'
  return(result)
}
