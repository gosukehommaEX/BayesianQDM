#' Find Optimal Go/NoGo Thresholds for a Single Binary Endpoint
#'
#' Computes the optimal Go threshold \eqn{\gamma_1} and NoGo threshold
#' \eqn{\gamma_2} for a single binary endpoint by searching over a fine grid of
#' candidate values.  The search follows a two-stage approach: posterior or
#' predictive probabilities are precomputed once for every possible outcome pair
#' \eqn{(y_1, y_2)}, and the operating characteristics under each candidate
#' \eqn{\gamma} are obtained by a fast weighted summation without additional
#' probability evaluations.
#'
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or
#'        \code{'external'}.
#' @param theta.TV A numeric scalar in \code{(-1, 1)} giving the Target Value
#'        (TV) threshold for the treatment effect.  Required when
#'        \code{prob = 'posterior'}; set to \code{NULL} otherwise.
#' @param theta.MAV A numeric scalar in \code{(-1, 1)} giving the Minimum
#'        Acceptable Value (MAV) threshold.  Must satisfy
#'        \code{theta.TV > theta.MAV}.  Required when
#'        \code{prob = 'posterior'}; set to \code{NULL} otherwise.
#' @param theta.NULL A numeric scalar in \code{(-1, 1)} giving the null
#'        hypothesis threshold used for the predictive probability.  Required
#'        when \code{prob = 'predictive'}; set to \code{NULL} otherwise.
#' @param pi1 A numeric scalar in \code{(0, 1)} giving the true treatment
#'        response rate for the scenario used to search for both
#'        \eqn{\gamma_1} and \eqn{\gamma_2}.
#' @param pi2 A numeric scalar in \code{(0, 1)} giving the true control
#'        response rate.  Set to \code{NULL} for
#'        \code{design = 'uncontrolled'}.
#' @param target_go A numeric scalar in \code{(0, 1)} giving the upper bound
#'        on \eqn{\Pr(\mathrm{Go})}.  The optimal \eqn{\gamma_1} is the
#'        smallest value in \code{gamma_grid} such that
#'        \eqn{\Pr(\mathrm{Go} \mid \pi_1, \pi_2)} is strictly less than
#'        \code{target_go}.
#' @param target_nogo A numeric scalar in \code{(0, 1)} giving the upper bound
#'        on \eqn{\Pr(\mathrm{NoGo})}.  The optimal \eqn{\gamma_2} is the
#'        largest value in \code{gamma_grid} such that
#'        \eqn{\Pr(\mathrm{NoGo} \mid \pi_1, \pi_2)} is strictly less than
#'        \code{target_nogo}.
#' @param n1 A positive integer giving the number of patients in group 1
#'        (treatment) in the PoC trial.
#' @param n2 A positive integer giving the number of patients in group 2
#'        (control) in the PoC trial.
#' @param a1 A positive numeric scalar giving the first shape parameter of the
#'        Beta prior for group 1.
#' @param a2 A positive numeric scalar giving the first shape parameter of the
#'        Beta prior for group 2.
#' @param b1 A positive numeric scalar giving the second shape parameter of the
#'        Beta prior for group 1.
#' @param b2 A positive numeric scalar giving the second shape parameter of the
#'        Beta prior for group 2.
#' @param z A non-negative integer giving the hypothetical number of responders
#'        in the control group.  Required when \code{design = 'uncontrolled'};
#'        set to \code{NULL} otherwise.
#' @param m1 A positive integer giving the future sample size for group 1.
#'        Required when \code{prob = 'predictive'}; set to \code{NULL}
#'        otherwise.
#' @param m2 A positive integer giving the future sample size for group 2.
#'        Required when \code{prob = 'predictive'}; set to \code{NULL}
#'        otherwise.
#' @param ne1 A positive integer giving the number of patients in group 1 of
#'        the external data set.  Required when \code{design = 'external'};
#'        set to \code{NULL} otherwise.
#' @param ne2 A positive integer giving the number of patients in group 2 of
#'        the external data set.  Required when \code{design = 'external'};
#'        set to \code{NULL} otherwise.
#' @param ye1 A non-negative integer giving the number of responders in group 1
#'        of the external data set.  Required when \code{design = 'external'};
#'        set to \code{NULL} otherwise.
#' @param ye2 A non-negative integer giving the number of responders in group 2
#'        of the external data set.  Required when \code{design = 'external'};
#'        set to \code{NULL} otherwise.
#' @param ae1 A numeric scalar in \code{(0, 1]} giving the power prior weight
#'        for group 1.  Required when \code{design = 'external'}; set to
#'        \code{NULL} otherwise.
#' @param ae2 A numeric scalar in \code{(0, 1]} giving the power prior weight
#'        for group 2.  Required when \code{design = 'external'}; set to
#'        \code{NULL} otherwise.
#' @param gamma_grid A numeric vector of candidate threshold values in
#'        \code{(0, 1)} to search over.  Defaults to
#'        \code{seq(0.01, 0.99, by = 0.01)}.
#'
#' @return A list of class \code{getgamma1bin} with the following elements:
#' \describe{
#'   \item{gamma1}{Optimal Go threshold: the smallest value in
#'         \code{gamma_grid} such that \eqn{\Pr(\mathrm{Go})} is strictly less
#'         than \code{target_go}.  \code{NA} if no such value exists.}
#'   \item{gamma2}{Optimal NoGo threshold: the largest value in
#'         \code{gamma_grid} such that \eqn{\Pr(\mathrm{NoGo})} is strictly
#'         less than \code{target_nogo}.  \code{NA} if no such value exists.}
#'   \item{PrGo_at_gamma1}{Pr(Go) evaluated at the optimal \eqn{\gamma_1}.
#'         \code{NA} if \code{gamma1} is \code{NA}.}
#'   \item{PrNoGo_at_gamma2}{Pr(NoGo) evaluated at the optimal \eqn{\gamma_2}.
#'         \code{NA} if \code{gamma2} is \code{NA}.}
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
#'         \eqn{(y_1, y_2)} are enumerated, and the posterior or predictive
#'         probability \eqn{g(y_1, y_2)} is computed once for every pair
#'         using \code{\link{pbayespostpred1bin}}.  This step is independent
#'         of \eqn{\gamma}.
#'   \item \strong{Gamma sweep}: For each candidate \eqn{\gamma} in
#'         \code{gamma_grid}, Pr(Go) and Pr(NoGo) are computed as weighted
#'         sums of the precomputed indicators, where the weights are binomial
#'         probabilities under \code{pi1} and \code{pi2}.  No further
#'         probability evaluations are required at this stage.
#' }
#' For \code{prob = 'posterior'}, two precomputation vectors are used: one
#' evaluated at \code{theta.TV} (Go criterion, \code{lower.tail = FALSE}) and
#' one at \code{theta.MAV} (NoGo criterion, \code{lower.tail = TRUE}).  For
#' \code{prob = 'predictive'}, both criteria use \code{theta.NULL} with the
#' respective tail direction.
#'
#' @examples
#' # Example 1: Controlled design, posterior probability
#' getgamma1bin(
#'   prob = 'posterior', design = 'controlled',
#'   theta.TV = 0.20, theta.MAV = 0.05, theta.NULL = NULL,
#'   pi1 = 0.15, pi2 = 0.15,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n1 = 12L, n2 = 12L,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   z = NULL, m1 = NULL, m2 = NULL,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
#' )
#'
#' # Example 2: Uncontrolled design, posterior probability
#' getgamma1bin(
#'   prob = 'posterior', design = 'uncontrolled',
#'   theta.TV = 0.20, theta.MAV = 0.05, theta.NULL = NULL,
#'   pi1 = 0.15, pi2 = NULL,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n1 = 12L, n2 = 12L,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   z = 3L, m1 = NULL, m2 = NULL,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
#' )
#'
#' # Example 3: Controlled design, predictive probability
#' getgamma1bin(
#'   prob = 'predictive', design = 'controlled',
#'   theta.TV = NULL, theta.MAV = NULL, theta.NULL = 0.10,
#'   pi1 = 0.15, pi2 = 0.15,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n1 = 12L, n2 = 12L,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   z = NULL, m1 = 30L, m2 = 30L,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
#' )
#'
#' # Example 4: External design, posterior probability
#' getgamma1bin(
#'   prob = 'posterior', design = 'external',
#'   theta.TV = 0.20, theta.MAV = 0.05, theta.NULL = NULL,
#'   pi1 = 0.15, pi2 = 0.15,
#'   target_go = 0.05, target_nogo = 0.20,
#'   n1 = 12L, n2 = 12L,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   z = NULL, m1 = NULL, m2 = NULL,
#'   ne1 = 15L, ne2 = 15L, ye1 = 6L, ye2 = 4L, ae1 = 0.5, ae2 = 0.5
#' )
#'
#' @importFrom stats dbinom
#' @export
getgamma1bin <- function(prob = 'posterior', design = 'controlled',
                         theta.TV = NULL, theta.MAV = NULL, theta.NULL = NULL,
                         pi1, pi2 = NULL,
                         target_go, target_nogo,
                         n1, n2,
                         a1, a2, b1, b2,
                         z = NULL, m1 = NULL, m2 = NULL,
                         ne1 = NULL, ne2 = NULL,
                         ye1 = NULL, ye2 = NULL,
                         ae1 = NULL, ae2 = NULL,
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
    if (is.null(theta.TV) || is.null(theta.MAV)) {
      stop("'theta.TV' and 'theta.MAV' must be non-NULL when prob = 'posterior'")
    }
    if (!is.numeric(theta.TV) || length(theta.TV) != 1L || is.na(theta.TV) ||
        theta.TV <= -1 || theta.TV >= 1) {
      stop("'theta.TV' must be a single numeric value in (-1, 1)")
    }
    if (!is.numeric(theta.MAV) || length(theta.MAV) != 1L || is.na(theta.MAV) ||
        theta.MAV <= -1 || theta.MAV >= 1) {
      stop("'theta.MAV' must be a single numeric value in (-1, 1)")
    }
    if (theta.TV <= theta.MAV) {
      stop("'theta.TV' must be strictly greater than 'theta.MAV'")
    }
  } else {
    if (is.null(theta.NULL)) {
      stop("'theta.NULL' must be non-NULL when prob = 'predictive'")
    }
    if (!is.numeric(theta.NULL) || length(theta.NULL) != 1L ||
        is.na(theta.NULL) || theta.NULL <= -1 || theta.NULL >= 1) {
      stop("'theta.NULL' must be a single numeric value in (-1, 1)")
    }
    if (is.null(m1) || is.null(m2)) {
      stop("'m1' and 'm2' must be non-NULL when prob = 'predictive'")
    }
    for (nm in c("m1", "m2")) {
      val <- get(nm)
      if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
          val != floor(val) || val < 1L) {
        stop(paste0("'", nm, "' must be a single positive integer"))
      }
    }
  }

  if (!is.numeric(pi1) || length(pi1) != 1L || is.na(pi1) ||
      pi1 <= 0 || pi1 >= 1) {
    stop("'pi1' must be a single numeric value in (0, 1)")
  }

  if (design != 'uncontrolled') {
    if (is.null(pi2) || !is.numeric(pi2) || length(pi2) != 1L ||
        is.na(pi2) || pi2 <= 0 || pi2 >= 1) {
      stop("'pi2' must be a single numeric value in (0, 1) for controlled or external design")
    }
  }

  for (nm in c("target_go", "target_nogo")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val <= 0 || val >= 1) {
      stop(paste0("'", nm, "' must be a single numeric value in (0, 1)"))
    }
  }

  for (nm in c("n1", "n2")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val != floor(val) || val < 1L) {
      stop(paste0("'", nm, "' must be a single positive integer"))
    }
  }

  for (nm in c("a1", "a2", "b1", "b2")) {
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
    for (nm in c("ne1", "ne2", "ye1", "ye2")) {
      val <- get(nm)
      if (is.null(val)) {
        stop(paste0("'", nm, "' must be non-NULL when design = 'external'"))
      }
    }
    if (!is.numeric(ne1) || length(ne1) != 1L || is.na(ne1) ||
        ne1 != floor(ne1) || ne1 < 1L) {
      stop("'ne1' must be a single positive integer")
    }
    if (!is.numeric(ne2) || length(ne2) != 1L || is.na(ne2) ||
        ne2 != floor(ne2) || ne2 < 1L) {
      stop("'ne2' must be a single positive integer")
    }
    if (!is.numeric(ye1) || length(ye1) != 1L || is.na(ye1) ||
        ye1 != floor(ye1) || ye1 < 0L || ye1 > ne1) {
      stop("'ye1' must be a non-negative integer not exceeding 'ne1'")
    }
    if (!is.numeric(ye2) || length(ye2) != 1L || is.na(ye2) ||
        ye2 != floor(ye2) || ye2 < 0L || ye2 > ne2) {
      stop("'ye2' must be a non-negative integer not exceeding 'ne2'")
    }
    for (nm in c("ae1", "ae2")) {
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
  # Stage 1: Precompute g(y1, y2) for all outcome combinations (gamma-free)
  # ---------------------------------------------------------------------------
  theta0_go   <- if (prob == 'posterior') theta.TV  else theta.NULL
  theta0_nogo <- if (prob == 'posterior') theta.MAV else theta.NULL

  if (design == 'uncontrolled') {
    # Enumerate y1 only; y2 is fixed at z inside pbayespostpred1bin
    all_y1 <- 0L:n1

    gPost_Go <- pbayespostpred1bin(
      prob = prob, design = design, theta0 = theta0_go,
      n1 = n1, n2 = n2, y1 = all_y1, y2 = NULL,
      a1 = a1, a2 = a2, b1 = b1, b2 = b2,
      m1 = m1, m2 = m2, z = z,
      ne1 = ne1, ne2 = ne2, ye1 = ye1, ye2 = ye2, ae1 = ae1, ae2 = ae2,
      lower.tail = FALSE
    )

    gPost_NoGo <- pbayespostpred1bin(
      prob = prob, design = design, theta0 = theta0_nogo,
      n1 = n1, n2 = n2, y1 = all_y1, y2 = NULL,
      a1 = a1, a2 = a2, b1 = b1, b2 = b2,
      m1 = m1, m2 = m2, z = z,
      ne1 = ne1, ne2 = ne2, ye1 = ye1, ye2 = ye2, ae1 = ae1, ae2 = ae2,
      lower.tail = TRUE
    )

    # Binomial weights under the specified scenario
    w <- dbinom(all_y1, n1, pi1)

  } else {
    # Enumerate all (y1, y2) pairs for controlled and external designs
    grid    <- expand.grid(y1 = 0L:n1, y2 = 0L:n2)
    all_y1e <- grid$y1
    all_y2e <- grid$y2

    gPost_Go <- pbayespostpred1bin(
      prob = prob, design = design, theta0 = theta0_go,
      n1 = n1, n2 = n2, y1 = all_y1e, y2 = all_y2e,
      a1 = a1, a2 = a2, b1 = b1, b2 = b2,
      m1 = m1, m2 = m2, z = NULL,
      ne1 = ne1, ne2 = ne2, ye1 = ye1, ye2 = ye2, ae1 = ae1, ae2 = ae2,
      lower.tail = FALSE
    )

    gPost_NoGo <- pbayespostpred1bin(
      prob = prob, design = design, theta0 = theta0_nogo,
      n1 = n1, n2 = n2, y1 = all_y1e, y2 = all_y2e,
      a1 = a1, a2 = a2, b1 = b1, b2 = b2,
      m1 = m1, m2 = m2, z = NULL,
      ne1 = ne1, ne2 = ne2, ye1 = ye1, ye2 = ye2, ae1 = ae1, ae2 = ae2,
      lower.tail = TRUE
    )

    # Binomial weights under the specified scenario
    w <- dbinom(all_y1e, n1, pi1) * dbinom(all_y2e, n2, pi2)
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
  # Stage 3: Find optimal gamma1 and gamma2
  # ---------------------------------------------------------------------------

  # gamma1: smallest gamma such that Pr(Go) < target_go
  idx1 <- which(PrGo_grid < target_go)
  if (length(idx1) == 0L) {
    gamma1         <- NA_real_
    PrGo_at_gamma1 <- NA_real_
  } else {
    gamma1         <- gamma_grid[min(idx1)]
    PrGo_at_gamma1 <- PrGo_grid[min(idx1)]
  }

  # gamma2: largest gamma such that Pr(NoGo) < target_nogo
  idx2 <- which(PrNoGo_grid < target_nogo)
  if (length(idx2) == 0L) {
    gamma2            <- NA_real_
    PrNoGo_at_gamma2  <- NA_real_
  } else {
    gamma2            <- gamma_grid[max(idx2)]
    PrNoGo_at_gamma2  <- PrNoGo_grid[max(idx2)]
  }

  # ---------------------------------------------------------------------------
  # Build and return result
  # ---------------------------------------------------------------------------
  result <- list(
    gamma1           = gamma1,
    gamma2           = gamma2,
    PrGo_at_gamma1   = PrGo_at_gamma1,
    PrNoGo_at_gamma2 = PrNoGo_at_gamma2,
    gamma_grid       = gamma_grid,
    PrGo_grid        = PrGo_grid,
    PrNoGo_grid      = PrNoGo_grid
  )

  class(result) <- 'getgamma1bin'
  return(result)
}
