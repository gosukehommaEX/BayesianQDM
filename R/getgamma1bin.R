#' Find Optimal Go/NoGo Thresholds for a Single Binary Endpoint
#'
#' Computes the optimal Go threshold \eqn{\gamma_go} and NoGo threshold
#' \eqn{\gamma_nogo} for a single binary endpoint by searching over a fine grid of
#' candidate values.  The search follows a two-stage approach: posterior or
#' predictive probabilities are precomputed once for every possible outcome pair
#' \eqn{(y_t, y_c)}, and the operating characteristics under each candidate
#' \eqn{\gamma} are obtained by a fast weighted summation without additional
#' probability evaluations.
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
#' @param pi_t A numeric scalar in \code{(0, 1)} giving the true treatment
#'        response rate for the scenario used to search for both
#'        \eqn{\gamma_go} and \eqn{\gamma_nogo}.
#' @param pi_c A numeric scalar in \code{(0, 1)} giving the true control
#'        response rate.  Set to \code{NULL} for
#'        \code{design = 'uncontrolled'}.
#' @param target_go A numeric scalar in \code{(0, 1)} giving the target value
#'        for \eqn{\Pr(\mathrm{Go})} used to determine the optimal
#'        \eqn{\gamma_go}.  The comparison operator applied is specified by
#'        \code{crit_go}.
#' @param target_nogo A numeric scalar in \code{(0, 1)} giving the target
#'        value for \eqn{\Pr(\mathrm{NoGo})} used to determine the optimal
#'        \eqn{\gamma_nogo}.  The comparison operator applied is specified by
#'        \code{crit_nogo}.
#' @param crit_go A character string specifying the comparison operator
#'        applied to \eqn{\Pr(\mathrm{Go})} when searching for
#'        \eqn{\gamma_go}.  Must be one of \code{"<"}, \code{"<="}, \code{">"},
#'        or \code{">="}. Default is \code{"<"}.
#' @param crit_nogo A character string specifying the comparison operator
#'        applied to \eqn{\Pr(\mathrm{NoGo})} when searching for
#'        \eqn{\gamma_nogo}.  Must be one of \code{"<"}, \code{"<="}, \code{">"},
#'        or \code{">="}. Default is \code{"<"}.
#' @param sel_go A character string specifying whether to select the
#'        \code{"smallest"} or \code{"largest"} value in \code{gamma_grid}
#'        among those satisfying the \code{crit_go} criterion.
#'        Default is \code{"smallest"}.
#' @param sel_nogo A character string specifying whether to select the
#'        \code{"smallest"} or \code{"largest"} value in \code{gamma_grid}
#'        among those satisfying the \code{crit_nogo} criterion.
#'        Default is \code{"largest"}.
#' @param n_t A positive integer giving the number of patients in the
#'        treatment group (treatment) in the PoC trial.
#' @param n_c A positive integer giving the number of patients in the
#'        control group (control) in the PoC trial.
#' @param a_t A positive numeric scalar giving the first shape parameter of the
#'        Beta prior for the treatment group.
#' @param a_c A positive numeric scalar giving the first shape parameter of the
#'        Beta prior for the control group.
#' @param b_t A positive numeric scalar giving the second shape parameter of the
#'        Beta prior for the treatment group.
#' @param b_c A positive numeric scalar giving the second shape parameter of the
#'        Beta prior for the control group.
#' @param z A non-negative integer giving the hypothetical number of responders
#'        in the control group.  Required when \code{design = 'uncontrolled'};
#'        set to \code{NULL} otherwise.
#' @param m_t A positive integer giving the future sample size for the treatment
#'        group. Required when \code{prob = 'predictive'}; set to \code{NULL}
#'        otherwise.
#' @param m_c A positive integer giving the future sample size for the control
#'        group. Required when \code{prob = 'predictive'}; set to \code{NULL}
#'        otherwise.
#' @param ne_t A positive integer giving the number of patients in the treatment
#'        group of the external data set.  Required when \code{design = 'external'};
#'        set to \code{NULL} otherwise.
#' @param ne_c A positive integer giving the number of patients in the control
#'        group of the external data set.  Required when \code{design = 'external'};
#'        set to \code{NULL} otherwise.
#' @param ye_t A non-negative integer giving the number of responders in the
#'        treatment group of the external data set.  Required when
#'        \code{design = 'external'}; set to \code{NULL} otherwise.
#' @param ye_c A non-negative integer giving the number of responders in the
#'        control group of the external data set.  Required when
#'        \code{design = 'external'}; set to \code{NULL} otherwise.
#' @param ae_t A numeric scalar in \code{(0, 1]} giving the power prior weight
#'        for the treatment group.  Required when \code{design = 'external'};
#'        set to \code{NULL} otherwise.
#' @param ae_c A numeric scalar in \code{(0, 1]} giving the power prior weight
#'        for the control group.  Required when \code{design = 'external'};
#'        set to \code{NULL} otherwise.
#' @param gamma_grid A numeric vector of candidate threshold values in
#'        \code{(0, 1)} to search over.  Defaults to
#'        \code{seq(0.01, 0.99, by = 0.01)}.
#'
#' @return A list of class \code{getgamma1bin} with the following elements:
#' \describe{
#'   \item{gamma_go}{Optimal Go threshold selected from \code{gamma_grid}
#'         according to \code{crit_go} and \code{sel_go}.
#'         \code{NA} if no value satisfies the criterion.}
#'   \item{gamma_nogo}{Optimal NoGo threshold selected from \code{gamma_grid}
#'         according to \code{crit_nogo} and \code{sel_nogo}.
#'         \code{NA} if no value satisfies the criterion.}
#'   \item{PrGo_at_gamma_go}{Pr(Go) evaluated at the optimal \eqn{\gamma_go}.
#'         \code{NA} if \code{gamma_go} is \code{NA}.}
#'   \item{PrNoGo_at_gamma_nogo}{Pr(NoGo) evaluated at the optimal \eqn{\gamma_nogo}.
#'         \code{NA} if \code{gamma_nogo} is \code{NA}.}
#'   \item{gamma_grid}{The candidate grid used for the search.}
#'   \item{PrGo_grid}{Numeric vector of Pr(Go) values over \code{gamma_grid}.}
#'   \item{PrNoGo_grid}{Numeric vector of Pr(NoGo) values over
#'         \code{gamma_grid}.}
#' }
#'
#' @details
#' The function uses a two-stage precompute-then-sweep strategy for efficiency:
#' \enumerate{
#'   \item \strong{Precomputation}: All possible outcome pairs
#'         \eqn{(y_t, y_c)} are enumerated, and the posterior or predictive
#'         probability \eqn{g(y_t, y_c)} is computed once for every pair
#'         using \code{\link{pbayespostpred1bin}}.  This step is independent
#'         of \eqn{\gamma}.
#'   \item \strong{Gamma sweep}: For each candidate \eqn{\gamma} in
#'         \code{gamma_grid}, Pr(Go) and Pr(NoGo) are computed as weighted
#'         sums of the precomputed indicators, where the weights are binomial
#'         probabilities under \code{pi_t} and \code{pi_c}.  No further
#'         probability evaluations are required at this stage.
#' }
#' For \code{prob = 'posterior'}, two precomputation vectors are used: one
#' evaluated at \code{theta_TV} (Go criterion, \code{lower.tail = FALSE}) and
#' one at \code{theta_MAV} (NoGo criterion, \code{lower.tail = TRUE}).  For
#' \code{prob = 'predictive'}, both criteria use \code{theta_NULL} with the
#' respective tail direction.
#'
#' The optimal \eqn{\gamma_go} is the \code{sel_go} (\code{"smallest"} or
#' \code{"largest"}) value in \code{gamma_grid} for which
#' \eqn{\Pr(\mathrm{Go})} satisfies the \code{crit_go} comparison against
#' \code{target_go}.  Analogously, the optimal \eqn{\gamma_nogo} is the
#' \code{sel_nogo} value satisfying \code{crit_nogo} against
#' \code{target_nogo}.
#'
#' @examples
#' # Example 1: Controlled design, posterior probability
#' # gamma_go: smallest gamma such that Pr(Go) < 0.05
#' # gamma_nogo: largest  gamma such that Pr(NoGo) < 0.20
#' \dontrun{
#' getgamma1bin(
#'   prob = 'posterior', design = 'controlled',
#'   theta_TV = 0.20, theta_MAV = 0.05, theta_NULL = NULL,
#'   pi_t = 0.15, pi_c = 0.15,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n_t = 12L, n_c = 12L,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = NULL, m_t = NULL, m_c = NULL,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, ae_t = NULL, ae_c = NULL
#' )
#' }
#'
#' # Example 2: Uncontrolled design, posterior probability
#' \dontrun{
#' getgamma1bin(
#'   prob = 'posterior', design = 'uncontrolled',
#'   theta_TV = 0.20, theta_MAV = 0.05, theta_NULL = NULL,
#'   pi_t = 0.15, pi_c = NULL,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n_t = 12L, n_c = 12L,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = 3L, m_t = NULL, m_c = NULL,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, ae_t = NULL, ae_c = NULL
#' )
#' }
#'
#' # Example 3: External design, posterior probability
#' \dontrun{
#' getgamma1bin(
#'   prob = 'posterior', design = 'external',
#'   theta_TV = 0.20, theta_MAV = 0.05, theta_NULL = NULL,
#'   pi_t = 0.15, pi_c = 0.15,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n_t = 12L, n_c = 12L,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = NULL, m_t = NULL, m_c = NULL,
#'   ne_t = 15L, ne_c = 15L, ye_t = 6L, ye_c = 4L, ae_t = 0.5, ae_c = 0.5
#' )
#' }
#'
#' # Example 4: Controlled design, predictive probability
#' \dontrun{
#' getgamma1bin(
#'   prob = 'predictive', design = 'controlled',
#'   theta_TV = NULL, theta_MAV = NULL, theta_NULL = 0.10,
#'   pi_t = 0.15, pi_c = 0.15,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n_t = 12L, n_c = 12L,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = NULL, m_t = 30L, m_c = 30L,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, ae_t = NULL, ae_c = NULL
#' )
#' }
#'
#' # Example 5: Uncontrolled design, predictive probability
#' \dontrun{
#' getgamma1bin(
#'   prob = 'predictive', design = 'uncontrolled',
#'   theta_TV = NULL, theta_MAV = NULL, theta_NULL = 0.10,
#'   pi_t = 0.15, pi_c = NULL,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n_t = 12L, n_c = 12L,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = 3L, m_t = 30L, m_c = 30L,
#'   ne_t = NULL, ne_c = NULL, ye_t = NULL, ye_c = NULL, ae_t = NULL, ae_c = NULL
#' )
#' }
#'
#' # Example 6: External design, predictive probability
#' \dontrun{
#' getgamma1bin(
#'   prob = 'predictive', design = 'external',
#'   theta_TV = NULL, theta_MAV = NULL, theta_NULL = 0.10,
#'   pi_t = 0.15, pi_c = 0.15,
#'   target_go = 0.05, target_nogo = 0.20,
#'   crit_go = '<', crit_nogo = '<',
#'   sel_go = 'smallest', sel_nogo = 'largest',
#'   n_t = 12L, n_c = 12L,
#'   a_t = 0.5, a_c = 0.5, b_t = 0.5, b_c = 0.5,
#'   z = NULL, m_t = 30L, m_c = 30L,
#'   ne_t = 15L, ne_c = 15L, ye_t = 6L, ye_c = 4L, ae_t = 0.5, ae_c = 0.5
#' )
#' }
#'
#'
#' @importFrom stats dbinom
#' @export
getgamma1bin <- function(prob = 'posterior', design = 'controlled',
                         theta_TV = NULL, theta_MAV = NULL, theta_NULL = NULL,
                         pi_t, pi_c = NULL,
                         target_go, target_nogo,
                         crit_go  = '<',        crit_nogo  = '<',
                         sel_go   = 'smallest', sel_nogo   = 'largest',
                         n_t, n_c,
                         a_t, a_c, b_t, b_c,
                         z = NULL, m_t = NULL, m_c = NULL,
                         ne_t = NULL, ne_c = NULL,
                         ye_t = NULL, ye_c = NULL,
                         ae_t = NULL, ae_c = NULL,
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

  if (!is.numeric(pi_t) || length(pi_t) != 1L || is.na(pi_t) ||
      pi_t <= 0 || pi_t >= 1) {
    stop("'pi_t' must be a single numeric value in (0, 1)")
  }

  if (design != 'uncontrolled') {
    if (is.null(pi_c) || !is.numeric(pi_c) || length(pi_c) != 1L ||
        is.na(pi_c) || pi_c <= 0 || pi_c >= 1) {
      stop("'pi_c' must be a single numeric value in (0, 1) for controlled or external design")
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
    for (nm in c("ae_t", "ae_c")) {
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
  n_gamma    <- length(gamma_grid)

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
      ne_t = ne_t, ne_c = ne_c, ye_t = ye_t, ye_c = ye_c, ae_t = ae_t, ae_c = ae_c,
      lower.tail = FALSE
    )

    gPost_NoGo <- pbayespostpred1bin(
      prob = prob, design = design, theta0 = theta0_nogo,
      n_t = n_t, n_c = n_c, y_t = all_y_t, y_c = NULL,
      a_t = a_t, a_c = a_c, b_t = b_t, b_c = b_c,
      m_t = m_t, m_c = m_c, z = z,
      ne_t = ne_t, ne_c = ne_c, ye_t = ye_t, ye_c = ye_c, ae_t = ae_t, ae_c = ae_c,
      lower.tail = TRUE
    )

    w <- dbinom(all_y_t, n_t, pi_t)

  } else {
    grid    <- expand.grid(y_t = 0L:n_t, y_c = 0L:n_c)
    all_y_t <- grid$y_t
    all_y_c <- grid$y_c

    gPost_Go <- pbayespostpred1bin(
      prob = prob, design = design, theta0 = theta0_go,
      n_t = n_t, n_c = n_c, y_t = all_y_t, y_c = all_y_c,
      a_t = a_t, a_c = a_c, b_t = b_t, b_c = b_c,
      m_t = m_t, m_c = m_c, z = NULL,
      ne_t = ne_t, ne_c = ne_c, ye_t = ye_t, ye_c = ye_c, ae_t = ae_t, ae_c = ae_c,
      lower.tail = FALSE
    )

    gPost_NoGo <- pbayespostpred1bin(
      prob = prob, design = design, theta0 = theta0_nogo,
      n_t = n_t, n_c = n_c, y_t = all_y_t, y_c = all_y_c,
      a_t = a_t, a_c = a_c, b_t = b_t, b_c = b_c,
      m_t = m_t, m_c = m_c, z = NULL,
      ne_t = ne_t, ne_c = ne_c, ye_t = ye_t, ye_c = ye_c, ae_t = ae_t, ae_c = ae_c,
      lower.tail = TRUE
    )

    w <- dbinom(all_y_t, n_t, pi_t) * dbinom(all_y_c, n_c, pi_c)
  }

  # ---------------------------------------------------------------------------
  # Stage 2: Sweep gamma_grid using vectorised outer comparison
  # ---------------------------------------------------------------------------
  # go_ind[i, k]   = 1 if gPost_Go[i]   >= gamma_grid[k]
  # nogo_ind[i, k] = 1 if gPost_NoGo[i] <= gamma_grid[k]
  go_ind   <- outer(gPost_Go,   gamma_grid, ">=")  # n_outcomes x n_gamma
  nogo_ind <- outer(gPost_NoGo, gamma_grid, "<=")  # n_outcomes x n_gamma

  PrGo_grid   <- colSums(w * go_ind)
  PrNoGo_grid <- colSums(w * nogo_ind)

  # ---------------------------------------------------------------------------
  # Stage 3: Apply comparison operators and select optimal gamma values
  # ---------------------------------------------------------------------------
  # Helper: apply comparison operator element-wise
  .compare <- function(x, op, val) {
    switch(op,
           '<'  = x <  val,
           '<=' = x <= val,
           '>'  = x >  val,
           '>=' = x >= val)
  }

  # Helper: select smallest or largest qualifying index
  .select_idx <- function(mask, sel) {
    idx <- which(mask)
    if (length(idx) == 0L) return(NA_integer_)
    if (sel == 'smallest') min(idx) else max(idx)
  }

  # gamma_go
  opt1 <- .select_idx(.compare(PrGo_grid,   crit_go,   target_go),   sel_go)
  if (is.na(opt1)) {
    gamma_go       <- NA_real_
    PrGo_at_gamma_go <- NA_real_
  } else {
    gamma_go       <- gamma_grid[opt1]
    PrGo_at_gamma_go <- PrGo_grid[opt1]
  }

  # gamma_nogo
  opt2 <- .select_idx(.compare(PrNoGo_grid, crit_nogo, target_nogo), sel_nogo)
  if (is.na(opt2)) {
    gamma_nogo        <- NA_real_
    PrNoGo_at_gamma_nogo  <- NA_real_
  } else {
    gamma_nogo        <- gamma_grid[opt2]
    PrNoGo_at_gamma_nogo  <- PrNoGo_grid[opt2]
  }

  # ---------------------------------------------------------------------------
  # Build and return result
  # ---------------------------------------------------------------------------
  result <- list(
    gamma_go           = gamma_go,
    gamma_nogo           = gamma_nogo,
    PrGo_at_gamma_go   = PrGo_at_gamma_go,
    PrNoGo_at_gamma_nogo = PrNoGo_at_gamma_nogo,
    gamma_grid       = gamma_grid,
    PrGo_grid        = PrGo_grid,
    PrNoGo_grid      = PrNoGo_grid
  )

  class(result) <- 'getgamma1bin'
  return(result)
}
