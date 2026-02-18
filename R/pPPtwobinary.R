#' Bayesian Posterior or Posterior Predictive Probability for a Clinical Trial
#' with Two Binary Endpoints
#'
#' Computes Bayesian posterior probability or posterior predictive probability
#' for clinical trials with two binary endpoints under a Dirichlet-multinomial
#' conjugate model. The function returns probabilities for nine decision regions
#' (posterior) or four decision regions (predictive) defined by target values
#' (TV) and minimum acceptable values (MAV) for both endpoints.
#' External data can be incorporated through power priors.
#'
#' @param prob A character string specifying the probability type.
#'        Must be \code{'posterior'} or \code{'predictive'}.
#' @param design A character string specifying the trial design.
#'        Must be \code{'controlled'}, \code{'uncontrolled'}, or \code{'external'}.
#' @param theta.TV1 A numeric scalar giving the target value (TV) threshold for
#'        Endpoint 1 (upper boundary of the Gray region).
#' @param theta.MAV1 A numeric scalar giving the minimum acceptable value (MAV)
#'        threshold for Endpoint 1 (lower boundary of the Gray region).
#'        For \code{prob = 'posterior'}, must satisfy \code{theta.TV1 > theta.MAV1}.
#'        For \code{prob = 'predictive'}, must equal \code{theta.TV1}.
#' @param theta.TV2 A numeric scalar giving the target value (TV) threshold for
#'        Endpoint 2.
#' @param theta.MAV2 A numeric scalar giving the minimum acceptable value (MAV)
#'        threshold for Endpoint 2.
#'        For \code{prob = 'posterior'}, must satisfy \code{theta.TV2 > theta.MAV2}.
#'        For \code{prob = 'predictive'}, must equal \code{theta.TV2}.
#' @param x1_00 A non-negative integer giving the count of (0,0) responses in
#'        group 1 (Endpoint 1 = 0, Endpoint 2 = 0).
#' @param x1_01 A non-negative integer giving the count of (0,1) responses in
#'        group 1 (Endpoint 1 = 0, Endpoint 2 = 1).
#' @param x1_10 A non-negative integer giving the count of (1,0) responses in
#'        group 1 (Endpoint 1 = 1, Endpoint 2 = 0).
#' @param x1_11 A non-negative integer giving the count of (1,1) responses in
#'        group 1 (Endpoint 1 = 1, Endpoint 2 = 1).
#' @param x2_00 A non-negative integer giving the count of (0,0) responses in
#'        group 2. For \code{design = 'uncontrolled'}, these are hypothetical
#'        control counts.
#' @param x2_01 A non-negative integer giving the count of (0,1) responses in
#'        group 2.
#' @param x2_10 A non-negative integer giving the count of (1,0) responses in
#'        group 2.
#' @param x2_11 A non-negative integer giving the count of (1,1) responses in
#'        group 2.
#' @param a1_00 A positive numeric scalar giving the Dirichlet prior parameter
#'        for the (0,0) response pattern in group 1.
#' @param a1_01 A positive numeric scalar giving the Dirichlet prior parameter
#'        for the (0,1) response pattern in group 1.
#' @param a1_10 A positive numeric scalar giving the Dirichlet prior parameter
#'        for the (1,0) response pattern in group 1.
#' @param a1_11 A positive numeric scalar giving the Dirichlet prior parameter
#'        for the (1,1) response pattern in group 1.
#' @param a2_00 A positive numeric scalar giving the Dirichlet prior parameter
#'        for the (0,0) response pattern in group 2.
#' @param a2_01 A positive numeric scalar giving the Dirichlet prior parameter
#'        for the (0,1) response pattern in group 2.
#' @param a2_10 A positive numeric scalar giving the Dirichlet prior parameter
#'        for the (1,0) response pattern in group 2.
#' @param a2_11 A positive numeric scalar giving the Dirichlet prior parameter
#'        for the (1,1) response pattern in group 2.
#' @param m1 A positive integer giving the future sample size for group 1.
#'        Required when \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param m2 A positive integer giving the future sample size for group 2.
#'        Required when \code{prob = 'predictive'}; otherwise set to \code{NULL}.
#' @param xe1_00 A non-negative integer giving the external group 1 count for
#'        pattern (0,0). Required when \code{design = 'external'} and external
#'        treatment data are used; otherwise set to \code{NULL}.
#' @param xe1_01 A non-negative integer giving the external group 1 count for
#'        pattern (0,1). Required for external treatment data; otherwise \code{NULL}.
#' @param xe1_10 A non-negative integer giving the external group 1 count for
#'        pattern (1,0). Required for external treatment data; otherwise \code{NULL}.
#' @param xe1_11 A non-negative integer giving the external group 1 count for
#'        pattern (1,1). Required for external treatment data; otherwise \code{NULL}.
#' @param xe2_00 A non-negative integer giving the external group 2 count for
#'        pattern (0,0). Required when \code{design = 'external'} and external
#'        control data are used; otherwise set to \code{NULL}.
#' @param xe2_01 A non-negative integer giving the external group 2 count for
#'        pattern (0,1). Required for external control data; otherwise \code{NULL}.
#' @param xe2_10 A non-negative integer giving the external group 2 count for
#'        pattern (1,0). Required for external control data; otherwise \code{NULL}.
#' @param xe2_11 A non-negative integer giving the external group 2 count for
#'        pattern (1,1). Required for external control data; otherwise \code{NULL}.
#' @param ae1 A numeric scalar in \code{(0, 1]} giving the power prior weight for
#'        group 1. Required when external treatment data are used; otherwise
#'        \code{NULL}.
#' @param ae2 A numeric scalar in \code{(0, 1]} giving the power prior weight for
#'        group 2. Required when external control data are used; otherwise
#'        \code{NULL}.
#' @param nMC A positive integer giving the number of Monte Carlo draws used to
#'        estimate region probabilities. Default is 10000. Larger values reduce
#'        Monte Carlo error.
#'
#' @return A named numeric vector of region probabilities.
#'         For \code{prob = 'posterior'}: length 9, named \code{R1}, ..., \code{R9},
#'         corresponding to regions defined by TV and MAV thresholds for both
#'         endpoints (column-major order: Endpoint 2 varies fastest).
#'         For \code{prob = 'predictive'}: length 4, named \code{R1}, ..., \code{R4}.
#'         All elements are non-negative and sum to 1.
#'
#' @details
#' The model uses a Dirichlet-multinomial conjugate analysis. The four response
#' categories are ordered as (0,0), (0,1), (1,0), (1,1).
#'
#' Prior: \eqn{p_k \sim \mathrm{Dirichlet}(\alpha_{k,00},\, \alpha_{k,01},\,
#' \alpha_{k,10},\, \alpha_{k,11})}.
#'
#' Posterior after observing \eqn{x_k}: \eqn{p_k \mid x_k \sim
#' \mathrm{Dirichlet}(\alpha_{k,00} + x_{k,00},\; \ldots,\;
#' \alpha_{k,11} + x_{k,11})}.
#'
#' Marginal response rates:
#' \eqn{\pi_{k1} = p_{k,10} + p_{k,11}} (Endpoint 1) and
#' \eqn{\pi_{k2} = p_{k,01} + p_{k,11}} (Endpoint 2).
#'
#' Treatment effects: \eqn{\theta_1 = \pi_{11} - \pi_{21}},
#' \eqn{\theta_2 = \pi_{12} - \pi_{22}}.
#'
#' For \code{design = 'external'}, the power prior augments the Dirichlet prior:
#' \deqn{\alpha_{k,lm}^* = \alpha_{k,lm} + ae_k \cdot xe_{k,lm}.}
#'
#' For \code{prob = 'posterior'} the nine regions are (column-major order with
#' Endpoint 2 varying fastest):
#' \itemize{
#'   \item R1: \eqn{\theta_1 > TV_1} AND \eqn{\theta_2 > TV_2}
#'   \item R2: \eqn{\theta_1 > TV_1} AND \eqn{TV_2 \ge \theta_2 > MAV_2}
#'   \item R3: \eqn{\theta_1 > TV_1} AND \eqn{\theta_2 \le MAV_2}
#'   \item R4: \eqn{TV_1 \ge \theta_1 > MAV_1} AND \eqn{\theta_2 > TV_2}
#'   \item R5: \eqn{TV_1 \ge \theta_1 > MAV_1} AND \eqn{TV_2 \ge \theta_2 > MAV_2}
#'   \item R6: \eqn{TV_1 \ge \theta_1 > MAV_1} AND \eqn{\theta_2 \le MAV_2}
#'   \item R7: \eqn{\theta_1 \le MAV_1} AND \eqn{\theta_2 > TV_2}
#'   \item R8: \eqn{\theta_1 \le MAV_1} AND \eqn{TV_2 \ge \theta_2 > MAV_2}
#'   \item R9: \eqn{\theta_1 \le MAV_1} AND \eqn{\theta_2 \le MAV_2}
#' }
#'
#' @examples
#' # Example 1: Posterior probability for controlled design
#' pPPtwobinary(
#'   prob = 'posterior', design = 'controlled',
#'   theta.TV1 = 0.20, theta.MAV1 = 0.10,
#'   theta.TV2 = 0.20, theta.MAV2 = 0.10,
#'   x1_00 = 5, x1_01 = 3, x1_10 = 4, x1_11 = 8,
#'   x2_00 = 8, x2_01 = 4, x2_10 = 5, x2_11 = 3,
#'   a1_00 = 0.5, a1_01 = 0.5, a1_10 = 0.5, a1_11 = 0.5,
#'   a2_00 = 0.5, a2_01 = 0.5, a2_10 = 0.5, a2_11 = 0.5,
#'   m1 = NULL, m2 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 10000
#' )
#'
#' # Example 2: Posterior predictive probability for controlled design
#' pPPtwobinary(
#'   prob = 'predictive', design = 'controlled',
#'   theta.TV1 = 0.15, theta.MAV1 = 0.15,
#'   theta.TV2 = 0.15, theta.MAV2 = 0.15,
#'   x1_00 = 3, x1_01 = 2, x1_10 = 3, x1_11 = 7,
#'   x2_00 = 6, x2_01 = 3, x2_10 = 4, x2_11 = 2,
#'   a1_00 = 1, a1_01 = 1, a1_10 = 1, a1_11 = 1,
#'   a2_00 = 1, a2_01 = 1, a2_10 = 1, a2_11 = 1,
#'   m1 = 50, m2 = 50,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 10000
#' )
#'
#' # Example 3: External control design with power prior
#' pPPtwobinary(
#'   prob = 'posterior', design = 'external',
#'   theta.TV1 = 0.15, theta.MAV1 = 0.08,
#'   theta.TV2 = 0.15, theta.MAV2 = 0.08,
#'   x1_00 = 3, x1_01 = 2, x1_10 = 3, x1_11 = 7,
#'   x2_00 = 5, x2_01 = 3, x2_10 = 4, x2_11 = 3,
#'   a1_00 = 0.5, a1_01 = 0.5, a1_10 = 0.5, a1_11 = 0.5,
#'   a2_00 = 0.5, a2_01 = 0.5, a2_10 = 0.5, a2_11 = 0.5,
#'   m1 = NULL, m2 = NULL,
#'   xe1_00 = 2, xe1_01 = 2, xe1_10 = 3, xe1_11 = 8,
#'   xe2_00 = 8, xe2_01 = 4, xe2_10 = 6, xe2_11 = 2,
#'   ae1 = 0.5, ae2 = 0.5,
#'   nMC = 10000
#' )
#'
#' # Example 4: Uncontrolled design with hypothetical control counts
#' pPPtwobinary(
#'   prob = 'posterior', design = 'uncontrolled',
#'   theta.TV1 = 0.20, theta.MAV1 = 0.10,
#'   theta.TV2 = 0.20, theta.MAV2 = 0.10,
#'   x1_00 = 5, x1_01 = 3, x1_10 = 4, x1_11 = 8,
#'   x2_00 = 8, x2_01 = 4, x2_10 = 5, x2_11 = 3,
#'   a1_00 = 0.5, a1_01 = 0.5, a1_10 = 0.5, a1_11 = 0.5,
#'   a2_00 = 0.5, a2_01 = 0.5, a2_10 = 0.5, a2_11 = 0.5,
#'   m1 = NULL, m2 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 10000
#' )
#'
#' @importFrom stats rbinom
#' @export
pPPtwobinary <- function(prob = 'posterior', design = 'controlled',
                         theta.TV1, theta.MAV1, theta.TV2, theta.MAV2,
                         x1_00, x1_01, x1_10, x1_11,
                         x2_00, x2_01, x2_10, x2_11,
                         a1_00, a1_01, a1_10, a1_11,
                         a2_00, a2_01, a2_10, a2_11,
                         m1, m2,
                         xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
                         xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
                         ae1 = NULL, ae2 = NULL,
                         nMC = 10000) {

  # --- Input validation ---
  if (!is.character(prob) || length(prob) != 1L ||
      !prob %in% c('posterior', 'predictive')) {
    stop("'prob' must be either 'posterior' or 'predictive'")
  }

  if (!is.character(design) || length(design) != 1L ||
      !design %in% c('controlled', 'uncontrolled', 'external')) {
    stop("'design' must be 'controlled', 'uncontrolled', or 'external'")
  }

  if (!is.numeric(nMC) || length(nMC) != 1L || is.na(nMC) ||
      nMC != floor(nMC) || nMC < 1L) {
    stop("'nMC' must be a single positive integer")
  }

  # Validate threshold parameters
  for (nm in c("theta.TV1", "theta.MAV1", "theta.TV2", "theta.MAV2")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val)) {
      stop(paste0("'", nm, "' must be a single numeric value"))
    }
  }

  if (prob == 'posterior') {
    if (theta.TV1 <= theta.MAV1) {
      stop("'theta.TV1' must be strictly greater than 'theta.MAV1'")
    }
    if (theta.TV2 <= theta.MAV2) {
      stop("'theta.TV2' must be strictly greater than 'theta.MAV2'")
    }
  } else {
    if (theta.TV1 != theta.MAV1) {
      stop("'theta.TV1' must equal 'theta.MAV1' when prob = 'predictive'")
    }
    if (theta.TV2 != theta.MAV2) {
      stop("'theta.TV2' must equal 'theta.MAV2' when prob = 'predictive'")
    }
  }

  # Validate observed count parameters (must be non-negative integers)
  for (nm in c("x1_00", "x1_01", "x1_10", "x1_11",
               "x2_00", "x2_01", "x2_10", "x2_11")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) ||
        val != floor(val) || val < 0L) {
      stop(paste0("'", nm, "' must be a single non-negative integer"))
    }
  }

  # Validate prior parameters (must be positive)
  for (nm in c("a1_00", "a1_01", "a1_10", "a1_11",
               "a2_00", "a2_01", "a2_10", "a2_11")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || val <= 0) {
      stop(paste0("'", nm, "' must be a single positive numeric value"))
    }
  }

  # Validate prob-specific parameters
  if (prob == 'predictive') {
    if (is.null(m1) || is.null(m2)) {
      stop("'m1' and 'm2' must be non-NULL when prob = 'predictive'")
    }
    if (!is.numeric(m1) || length(m1) != 1L || is.na(m1) ||
        m1 != floor(m1) || m1 < 1L) {
      stop("'m1' must be a single positive integer")
    }
    if (!is.numeric(m2) || length(m2) != 1L || is.na(m2) ||
        m2 != floor(m2) || m2 < 1L) {
      stop("'m2' must be a single positive integer")
    }
  }

  # Validate design-specific parameters
  if (design == 'external') {
    has_external1 <- !is.null(xe1_00) && !is.null(xe1_01) &&
      !is.null(xe1_10) && !is.null(xe1_11) && !is.null(ae1)
    has_external2 <- !is.null(xe2_00) && !is.null(xe2_01) &&
      !is.null(xe2_10) && !is.null(xe2_11) && !is.null(ae2)

    if (!has_external1 && !has_external2) {
      stop("For design = 'external', at least one complete set of external data (xe1_* + ae1 or xe2_* + ae2) must be provided")
    }

    if (!is.null(ae1)) {
      if (!is.numeric(ae1) || length(ae1) != 1L || is.na(ae1) ||
          ae1 <= 0 || ae1 > 1) {
        stop("'ae1' must be a single numeric value in (0, 1]")
      }
    }
    if (!is.null(ae2)) {
      if (!is.numeric(ae2) || length(ae2) != 1L || is.na(ae2) ||
          ae2 <= 0 || ae2 > 1) {
        stop("'ae2' must be a single numeric value in (0, 1]")
      }
    }
  }

  # --- Compute posterior Dirichlet parameters ---
  # Group 1: incorporate external treatment data if available
  xe1_contrib <- if (!is.null(ae1) && design == 'external') ae1 else 0
  alpha1_00_post <- a1_00 + x1_00 + xe1_contrib * ifelse(!is.null(xe1_00), xe1_00, 0)
  alpha1_01_post <- a1_01 + x1_01 + xe1_contrib * ifelse(!is.null(xe1_01), xe1_01, 0)
  alpha1_10_post <- a1_10 + x1_10 + xe1_contrib * ifelse(!is.null(xe1_10), xe1_10, 0)
  alpha1_11_post <- a1_11 + x1_11 + xe1_contrib * ifelse(!is.null(xe1_11), xe1_11, 0)

  # Group 2: incorporate external control data if available
  xe2_contrib <- if (!is.null(ae2) && design == 'external') ae2 else 0
  alpha2_00_post <- a2_00 + x2_00 + xe2_contrib * ifelse(!is.null(xe2_00), xe2_00, 0)
  alpha2_01_post <- a2_01 + x2_01 + xe2_contrib * ifelse(!is.null(xe2_01), xe2_01, 0)
  alpha2_10_post <- a2_10 + x2_10 + xe2_contrib * ifelse(!is.null(xe2_10), xe2_10, 0)
  alpha2_11_post <- a2_11 + x2_11 + xe2_contrib * ifelse(!is.null(xe2_11), xe2_11, 0)

  # --- Monte Carlo estimation ---
  if (prob == 'posterior') {

    # Sample from posterior Dirichlet distributions
    p1_samples <- rdirichlet(nMC, c(alpha1_00_post, alpha1_01_post,
                                    alpha1_10_post, alpha1_11_post))
    p2_samples <- rdirichlet(nMC, c(alpha2_00_post, alpha2_01_post,
                                    alpha2_10_post, alpha2_11_post))

    # Marginal response rates:
    #   pi_k1 = p_k10 + p_k11  (Endpoint 1 success)
    #   pi_k2 = p_k01 + p_k11  (Endpoint 2 success)
    theta1 <- (p1_samples[, 3] + p1_samples[, 4]) -
      (p2_samples[, 3] + p2_samples[, 4])
    theta2 <- (p1_samples[, 2] + p1_samples[, 4]) -
      (p2_samples[, 2] + p2_samples[, 4])

    # Assign each draw to one of 9 regions (column-major: Endpoint 2 varies fastest)
    R1 <- 1L * (theta1 >  theta.TV1) +
      2L * ((theta.TV1 >= theta1) & (theta1 > theta.MAV1)) +
      3L * (theta.MAV1 >= theta1)
    R2 <- 1L * (theta2 >  theta.TV2) +
      2L * ((theta.TV2 >= theta2) & (theta2 > theta.MAV2)) +
      3L * (theta.MAV2 >= theta2)
    R <- (R1 - 1L) * 3L + R2

    Pr_R       <- tabulate(R, nbins = 9L) / nMC
    names(Pr_R) <- paste0("R", 1:9)

  } else {

    # Sample from posterior Dirichlet distributions
    p1 <- rdirichlet(nMC, c(alpha1_00_post, alpha1_01_post,
                            alpha1_10_post, alpha1_11_post))
    p2 <- rdirichlet(nMC, c(alpha2_00_post, alpha2_01_post,
                            alpha2_10_post, alpha2_11_post))

    # Generate future multinomial counts via sequential binomial draws
    # (exact Dirichlet-Multinomial sampling: draw from categories 1-3 in turn,
    #  remainder goes to category 4)
    x1_future <- matrix(0L, nMC, 4L)
    x2_future <- matrix(0L, nMC, 4L)
    rem1       <- rep(m1, nMC)
    rem2       <- rep(m2, nMC)
    used_p1    <- rep(0, nMC)
    used_p2    <- rep(0, nMC)

    for (j in seq_len(3L)) {
      denom1 <- pmax(1 - used_p1, 0)
      denom2 <- pmax(1 - used_p2, 0)
      ratio1 <- pmin(pmax(ifelse(denom1 > 0, p1[, j] / denom1, 0), 0), 1)
      ratio2 <- pmin(pmax(ifelse(denom2 > 0, p2[, j] / denom2, 0), 0), 1)

      draw1 <- rbinom(nMC, rem1, ratio1)
      draw2 <- rbinom(nMC, rem2, ratio2)

      x1_future[, j] <- draw1
      x2_future[, j] <- draw2
      rem1            <- rem1 - draw1
      rem2            <- rem2 - draw2
      used_p1         <- used_p1 + p1[, j]
      used_p2         <- used_p2 + p2[, j]
    }
    x1_future[, 4L] <- rem1
    x2_future[, 4L] <- rem2

    # Marginal response proportions from future data
    theta1 <- (x1_future[, 3L] + x1_future[, 4L]) / m1 -
      (x2_future[, 3L] + x2_future[, 4L]) / m2
    theta2 <- (x1_future[, 2L] + x1_future[, 4L]) / m1 -
      (x2_future[, 2L] + x2_future[, 4L]) / m2

    # Assign each draw to one of 4 regions
    R1 <- 1L * (theta1 >  theta.TV1) + 2L * (theta.MAV1 >= theta1)
    R2 <- 1L * (theta2 >  theta.TV2) + 2L * (theta.MAV2 >= theta2)
    R  <- (R1 - 1L) * 2L + R2

    Pr_R       <- tabulate(R, nbins = 4L) / nMC
    names(Pr_R) <- paste0("R", 1:4)
  }

  return(Pr_R)
}
