#' Calculate Bayesian Posterior Probability or Posterior Predictive Probability
#' for a Clinical Trial with Two Binary Endpoints
#'
#' This function computes Bayesian posterior probability or posterior predictive
#' probability for clinical trials with two binary endpoints. The function supports
#' controlled and external control designs, using Dirichlet-multinomial conjugate
#' priors. External data can be incorporated through power priors. The function
#' returns probabilities for nine decision regions defined by target values (TV)
#' and minimum acceptable values (MAV) for both endpoints.
#'
#' @param prob A character string specifying the type of probability to calculate.
#'        Options are \code{'posterior'} for posterior probability or \code{'predictive'}
#'        for posterior predictive probability.
#' @param design A character string specifying the type of trial design. Options are
#'        \code{'controlled'} for randomized controlled trials, \code{'uncontrolled'}
#'        for single-arm studies, or \code{'external'} for designs incorporating external data.
#' @param theta.TV1 A numeric value representing the target value (TV) threshold for
#'        Endpoint 1. This defines the upper boundary for the "Gray" region.
#' @param theta.MAV1 A numeric value representing the minimum acceptable value (MAV)
#'        threshold for Endpoint 1. This defines the lower boundary for the "Gray" region.
#' @param theta.TV2 A numeric value representing the target value (TV) threshold for
#'        Endpoint 2. This defines the upper boundary for the "Gray" region.
#' @param theta.MAV2 A numeric value representing the minimum acceptable value (MAV)
#'        threshold for Endpoint 2. This defines the lower boundary for the "Gray" region.
#' @param x1_00 A non-negative integer representing the count of (0,0) responses
#'        in group 1 for the PoC trial (Endpoint 1 = 0, Endpoint 2 = 0).
#' @param x1_01 A non-negative integer representing the count of (0,1) responses
#'        in group 1 for the PoC trial (Endpoint 1 = 0, Endpoint 2 = 1).
#' @param x1_10 A non-negative integer representing the count of (1,0) responses
#'        in group 1 for the PoC trial (Endpoint 1 = 1, Endpoint 2 = 0).
#' @param x1_11 A non-negative integer representing the count of (1,1) responses
#'        in group 1 for the PoC trial (Endpoint 1 = 1, Endpoint 2 = 1).
#' @param x2_00 A non-negative integer representing the count of (0,0) responses in group 2.
#'        For \code{design = 'controlled'} or \code{'external'}: observed control data.
#'        For \code{design = 'uncontrolled'}: hypothetical control response count based on
#'        historical data or prior knowledge.
#' @param x2_01 A non-negative integer representing the count of (0,1) responses in group 2.
#'        For \code{design = 'controlled'} or \code{'external'}: observed control data.
#'        For \code{design = 'uncontrolled'}: hypothetical control response count.
#' @param x2_10 A non-negative integer representing the count of (1,0) responses in group 2.
#'        For \code{design = 'controlled'} or \code{'external'}: observed control data.
#'        For \code{design = 'uncontrolled'}: hypothetical control response count.
#' @param x2_11 A non-negative integer representing the count of (1,1) responses in group 2.
#'        For \code{design = 'controlled'} or \code{'external'}: observed control data.
#'        For \code{design = 'uncontrolled'}: hypothetical control response count.
#' @param a1_00 A positive numeric value for the Dirichlet prior parameter α1_00.
#' @param a1_01 A positive numeric value for the Dirichlet prior parameter α1_01.
#' @param a1_10 A positive numeric value for the Dirichlet prior parameter α1_10.
#' @param a1_11 A positive numeric value for the Dirichlet prior parameter α1_11.
#' @param a2_00 A positive numeric value for the Dirichlet prior parameter α2_00.
#'        For \code{design = 'uncontrolled'}: typically set equal to a1_00.
#' @param a2_01 A positive numeric value for the Dirichlet prior parameter α2_01.
#'        For \code{design = 'uncontrolled'}: typically set equal to a1_01.
#' @param a2_10 A positive numeric value for the Dirichlet prior parameter α2_10.
#'        For \code{design = 'uncontrolled'}: typically set equal to a1_10.
#' @param a2_11 A positive numeric value for the Dirichlet prior parameter α2_11.
#'        For \code{design = 'uncontrolled'}: typically set equal to a1_11.
#' @param m1 A positive integer representing the number of patients in group 1 for
#'        the future trial (required if \code{prob = 'predictive'}).
#' @param m2 A positive integer representing the number of patients in group 2 for
#'        the future trial (required if \code{prob = 'predictive'}).
#' @param xe1_00 A non-negative integer for external group 1 count (0,0)
#'        (required if \code{design = 'external'} and external treatment data used).
#' @param xe1_01 A non-negative integer for external group 1 count (0,1)
#'        (required if \code{design = 'external'} and external treatment data used).
#' @param xe1_10 A non-negative integer for external group 1 count (1,0)
#'        (required if \code{design = 'external'} and external treatment data used).
#' @param xe1_11 A non-negative integer for external group 1 count (1,1)
#'        (required if \code{design = 'external'} and external treatment data used).
#' @param xe2_00 A non-negative integer for external group 2 count (0,0)
#'        (required if \code{design = 'external'} and external control data used).
#' @param xe2_01 A non-negative integer for external group 2 count (0,1)
#'        (required if \code{design = 'external'} and external control data used).
#' @param xe2_10 A non-negative integer for external group 2 count (1,0)
#'        (required if \code{design = 'external'} and external control data used).
#' @param xe2_11 A non-negative integer for external group 2 count (1,1)
#'        (required if \code{design = 'external'} and external control data used).
#' @param ae1 A numeric value in (0, 1] representing the power prior scale parameter
#'        for group 1 (required if \code{design = 'external'} and external treatment data used).
#'        Controls the degree of borrowing: 0 = no borrowing, 1 = full borrowing.
#' @param ae2 A numeric value in (0, 1] representing the power prior scale parameter
#'        for group 2 (required if \code{design = 'external'} and external control data used).
#'        Controls the degree of borrowing: 0 = no borrowing, 1 = full borrowing.
#' @param nMC A positive integer representing the number of Monte Carlo samples
#'        for probability calculation (default = 10000). Larger values provide
#'        more accurate estimates but require more computation time.
#'
#' @return A named numeric vector of length 9 containing the posterior or posterior
#'         predictive probabilities for each of the nine decision regions (R1, ..., R9).
#'         The regions are defined by the treatment effects θ1 = π11 - π21 and
#'         θ2 = π12 - π22 relative to the thresholds:
#'         \itemize{
#'           \item R1: θ1 > TV1 AND θ2 > TV2 (both endpoints exceed target)
#'           \item R2: θ1 > TV1 AND TV2 ≥ θ2 > MAV2
#'           \item R3: θ1 > TV1 AND θ2 ≤ MAV2
#'           \item R4: TV1 ≥ θ1 > MAV1 AND θ2 > TV2
#'           \item R5: TV1 ≥ θ1 > MAV1 AND TV2 ≥ θ2 > MAV2 (both in Gray zone)
#'           \item R6: TV1 ≥ θ1 > MAV1 AND θ2 ≤ MAV2
#'           \item R7: θ1 ≤ MAV1 AND θ2 > TV2
#'           \item R8: θ1 ≤ MAV1 AND TV2 ≥ θ2 > MAV2
#'           \item R9: θ1 ≤ MAV1 AND θ2 ≤ MAV2 (both endpoints below minimum)
#'         }
#'
#' @details
#' The function computes probabilities based on Dirichlet-multinomial conjugate analysis:
#' \itemize{
#'   \item **Data structure**: Two binary endpoints (y_k1, y_k2) are aggregated into
#'         four response categories: (0,0), (0,1), (1,0), and (1,1)
#'   \item **Prior**: p_k ~ Dirichlet(αk_00, αk_01, αk_10, αk_11)
#'   \item **Posterior**: p_k | x_k ~ Dirichlet(αk_00 + xk_00, αk_01 + xk_01,
#'         αk_10 + xk_10, αk_11 + xk_11)
#'   \item **Treatment effects**: θ1 = π11 - π21, θ2 = π12 - π22, where
#'         πk1 = pk_10 + pk_11 (marginal probability for Endpoint 1),
#'         πk2 = pk_01 + pk_11 (marginal probability for Endpoint 2)
#' }
#'
#' **Note on correlation**: The two endpoints are correlated within each arm because
#' they share the p_11 component. This correlation is naturally captured through the
#' Dirichlet distribution and Monte Carlo sampling.
#'
#' For **external control designs**, power priors are used to incorporate historical
#' data. The effective prior becomes:
#' Dirichlet(αk_00 + aek × xek_00, αk_01 + aek × xek_01,
#'           αk_10 + aek × xek_10, αk_11 + aek × xek_11)
#'
#' **Posterior predictive distribution**:
#' tilde.x_k | x_k ~ Dirichlet-Multinomial(mk, αk_00 + xk_00, αk_01 + xk_01,
#'                                      αk_10 + xk_10, αk_11 + xk_11)
#'
#' **Monte Carlo estimation**: Due to the correlation between endpoints and the
#' complexity of the joint distribution, probabilities are estimated using Monte
#' Carlo simulation. Larger nMC values provide more accurate estimates.
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
#' # Example 4: Uncontrolled design (single-arm study)
#' # Use hypothetical control data with same priors as treatment
#' pPPtwobinary(
#'   prob = 'posterior', design = 'uncontrolled',
#'   theta.TV1 = 0.20, theta.MAV1 = 0.10,
#'   theta.TV2 = 0.20, theta.MAV2 = 0.10,
#'   x1_00 = 5, x1_01 = 3, x1_10 = 4, x1_11 = 8,
#'   x2_00 = 8, x2_01 = 4, x2_10 = 5, x2_11 = 3,  # Hypothetical control counts
#'   a1_00 = 0.5, a1_01 = 0.5, a1_10 = 0.5, a1_11 = 0.5,
#'   a2_00 = 0.5, a2_01 = 0.5, a2_10 = 0.5, a2_11 = 0.5,  # Same as a1_xx
#'   m1 = NULL, m2 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 10000
#' )
#'
#' # Example 5: Visualizing the 9 regions
#' \dontrun{
#' result <- pPPtwobinary(
#'   prob = 'posterior', design = 'controlled',
#'   theta.TV1 = 0.20, theta.MAV1 = 0.10,
#'   theta.TV2 = 0.20, theta.MAV2 = 0.10,
#'   x1_00 = 4, x1_01 = 3, x1_10 = 5, x1_11 = 13,
#'   x2_00 = 10, x2_01 = 5, x2_10 = 7, x2_11 = 3,
#'   a1_00 = 0.25, a1_01 = 0.25, a1_10 = 0.25, a1_11 = 0.25,
#'   a2_00 = 0.25, a2_01 = 0.25, a2_10 = 0.25, a2_11 = 0.25,
#'   m1 = NULL, m2 = NULL,
#'   xe1_00 = NULL, xe1_01 = NULL, xe1_10 = NULL, xe1_11 = NULL,
#'   xe2_00 = NULL, xe2_01 = NULL, xe2_10 = NULL, xe2_11 = NULL,
#'   ae1 = NULL, ae2 = NULL,
#'   nMC = 10000
#' )
#'
#' # Display as 3x3 matrix
#' matrix(result, nrow = 3, byrow = FALSE,
#'        dimnames = list(
#'          c("θ2 > TV2", "TV2 ≥ θ2 > MAV2", "MAV2 ≥ θ2"),
#'          c("θ1 > TV1", "TV1 ≥ θ1 > MAV1", "MAV1 ≥ θ1")
#'        ))
#' }
#'
#' @importFrom stats rmultinom
#' @export
pPPtwobinary <- function(prob = 'posterior', design = 'controlled',
                         theta.TV1, theta.MAV1, theta.TV2, theta.MAV2,
                         x1_00, x1_01, x1_10, x1_11,
                         x2_00, x2_01, x2_10, x2_11,
                         a1_00, a1_01, a1_10, a1_11,
                         a2_00, a2_01, a2_10, a2_11,
                         m1, m2,
                         xe1_00, xe1_01, xe1_10, xe1_11,
                         xe2_00, xe2_01, xe2_10, xe2_11,
                         ae1, ae2,
                         nMC = 10000) {

  # Validate prob parameter
  if(!prob %in% c('posterior', 'predictive')) {
    stop("prob must be either 'posterior' or 'predictive'")
  }

  # Validate design parameter
  if(!design %in% c('controlled', 'external', 'uncontrolled')) {
    stop("design must be 'controlled', 'external', or 'uncontrolled'")
  }

  # Validate threshold ordering
  if(prob == 'posterior') {
    if(theta.TV1 <= theta.MAV1) {
      stop("theta.TV1 must be greater than theta.MAV1")
    }
    if(theta.TV2 <= theta.MAV2) {
      stop("theta.TV2 must be greater than theta.MAV2")
    }
  } else if(prob =='predictive') {
    if(theta.TV1 != theta.MAV1) {
      stop("theta.TV1 must be same to the theta.MAV1")
    }
    if(theta.TV2 != theta.MAV2) {
      stop("theta.TV2 must be same to the theta.MAV2")
    }
  }

  # Validate parameter sets for posterior predictive probability
  if((prob == 'predictive') & (is.null(m1) | is.null(m2))) {
    stop('If you calculate posterior predictive probability, m1 and m2 must be non-null')
  }

  # Validate parameter sets for external design
  if(design == 'external') {
    # Check if at least one external data source is provided
    has_external1 <- !is.null(xe1_00) & !is.null(xe1_01) & !is.null(xe1_10) & !is.null(xe1_11) & !is.null(ae1)
    has_external2 <- !is.null(xe2_00) & !is.null(xe2_01) & !is.null(xe2_10) & !is.null(xe2_11) & !is.null(ae2)

    if(!has_external1 & !has_external2) {
      stop('If you use external design, at least one set of external data must be provided')
    }
  }

  # Calculate posterior Dirichlet parameters for group 1 (treatment)
  # Incorporate external data if design is 'external'
  alpha1_00_post <- a1_00 + x1_00 + (design == 'external') * ifelse(!is.null(ae1), xe1_00 * ae1, 0)
  alpha1_01_post <- a1_01 + x1_01 + (design == 'external') * ifelse(!is.null(ae1), xe1_01 * ae1, 0)
  alpha1_10_post <- a1_10 + x1_10 + (design == 'external') * ifelse(!is.null(ae1), xe1_10 * ae1, 0)
  alpha1_11_post <- a1_11 + x1_11 + (design == 'external') * ifelse(!is.null(ae1), xe1_11 * ae1, 0)

  # Calculate posterior Dirichlet parameters for group 2 (control)
  # Require observed control data for controlled and external designs
  if(design != 'uncontrolled') {
    if(is.null(x2_00) | is.null(x2_01) | is.null(x2_10) | is.null(x2_11)) {
      stop('For controlled and external designs, x2_00, x2_01, x2_10, and x2_11 must be non-null')
    }
  } else {
    # For uncontrolled design, x2 parameters represent hypothetical control data
    if(is.null(x2_00) | is.null(x2_01) | is.null(x2_10) | is.null(x2_11)) {
      stop('For uncontrolled design, x2_00, x2_01, x2_10, and x2_11 must be provided as hypothetical control counts')
    }
  }

  alpha2_00_post <- a2_00 + x2_00 + (design == 'external') * ifelse(!is.null(ae2), xe2_00 * ae2, 0)
  alpha2_01_post <- a2_01 + x2_01 + (design == 'external') * ifelse(!is.null(ae2), xe2_01 * ae2, 0)
  alpha2_10_post <- a2_10 + x2_10 + (design == 'external') * ifelse(!is.null(ae2), xe2_10 * ae2, 0)
  alpha2_11_post <- a2_11 + x2_11 + (design == 'external') * ifelse(!is.null(ae2), xe2_11 * ae2, 0)

  # Monte Carlo simulation for probability calculation
  if(prob == 'posterior') {
    # Sample from posterior Dirichlet distributions
    # Group 1 (treatment)
    p1_samples <- rdirichlet(nMC, c(alpha1_00_post, alpha1_01_post, alpha1_10_post, alpha1_11_post))

    # Group 2 (control)
    p2_samples <- rdirichlet(nMC, c(alpha2_00_post, alpha2_01_post, alpha2_10_post, alpha2_11_post))

    # Calculate marginal response rates for each endpoint
    # pi_k1 = p_k10 + p_k11 (Endpoint 1: Y_k1 = 1)
    # pi_k2 = p_k01 + p_k11 (Endpoint 2: Y_k2 = 1)
    pi1_1 <- p1_samples[, 3] + p1_samples[, 4]  # p1_10 + p1_11
    pi1_2 <- p1_samples[, 2] + p1_samples[, 4]  # p1_01 + p1_11
    pi2_1 <- p2_samples[, 3] + p2_samples[, 4]  # p2_10 + p2_11
    pi2_2 <- p2_samples[, 2] + p2_samples[, 4]  # p2_01 + p2_11

    # Calculate treatment effects
    theta1 <- pi1_1 - pi2_1
    theta2 <- pi1_2 - pi2_2

  } else if(prob == 'predictive') {
    # Sample from posterior predictive Dirichlet-Multinomial distributions
    # This requires sampling from Dirichlet then Multinomial

    # Sample p from posterior Dirichlet
    p1 <- rdirichlet(nMC, c(alpha1_00_post, alpha1_01_post, alpha1_10_post, alpha1_11_post))
    p2 <- rdirichlet(nMC, c(alpha2_00_post, alpha2_01_post, alpha2_10_post, alpha2_11_post))

    x1_future <- matrix(0L, nMC, 4)
    x2_future <- matrix(0L, nMC, 4)
    rem1 <- m1
    rem2 <- m2
    used_prob1 <- rep(0, nMC)
    used_prob2 <- rep(0, nMC)

    for (j in seq_len(3)) {
      denom1 <- pmax(1 - used_prob1, 0)
      denom2 <- pmax(1 - used_prob2, 0)

      ratio1 <- ifelse(denom1 > 0, p1[, j] / denom1, 0)
      ratio1 <- pmin(pmax(ratio1, 0), 1)

      ratio2 <- ifelse(denom2 > 0, p2[, j] / denom2, 0)
      ratio2 <- pmin(pmax(ratio2, 0), 1)

      draw1 <- rbinom(nMC, rem1, ratio1)
      draw2 <- rbinom(nMC, rem2, ratio2)

      x1_future[, j] <- draw1
      x2_future[, j] <- draw2

      rem1 <- rem1 - draw1
      rem2 <- rem2 - draw2

      used_prob1 <- used_prob1 + p1[, j]
      used_prob2 <- used_prob2 + p2[, j]
    }
    x1_future[, 4] <- rem1
    x2_future[, 4] <- rem2

    # Calculate marginal response rates from future data
    pi1_1_pred <- (x1_future[, 3] + x1_future[, 4]) / m1
    pi1_2_pred <- (x1_future[, 2] + x1_future[, 4]) / m1
    pi2_1_pred <- (x2_future[, 3] + x2_future[, 4]) / m2
    pi2_2_pred <- (x2_future[, 2] + x2_future[, 4]) / m2

    # Calculate treatment effects
    theta1 <- pi1_1_pred - pi2_1_pred
    theta2 <- pi1_2_pred - pi2_2_pred
  }

  if(prob == 'posterior') {
    # Calculate region indices (1-9) based on treatment effects and thresholds
    # Region index for Endpoint 1
    R1 <- 1 * (theta1 > theta.TV1) +
      2 * ((theta.TV1 >= theta1) & (theta1 > theta.MAV1)) +
      3 * (theta.MAV1 >= theta1)

    # Region index for Endpoint 2
    R2 <- 1 * (theta2 > theta.TV2) +
      2 * ((theta.TV2 >= theta2) & (theta2 > theta.MAV2)) +
      3 * (theta.MAV2 >= theta2)

    # Combined region index (1-9)
    R <- (R1 - 1) * 3 + R2

    # Calculate region probabilities Pr(R_l | D)
    Pr_R <- tabulate(R, nbins = 9) / nMC
    names(Pr_R) <- paste0("R", 1:9)
  } else if(prob =='predictive') {
    # Region index for Endpoint 1
    R1 <- 1 * (theta1 > theta.TV1) + 2 * (theta.MAV1 >= theta1)

    # Region index for Endpoint 2
    R2 <- 1 * (theta2 > theta.TV2) + 2 * (theta.MAV2 >= theta2)

    # Combined region index (1-4)
    R <- (R1 - 1) * 2 + R2

    # Calculate region probabilities Pr(R_l | D)
    Pr_R <- tabulate(R, nbins = 4) / nMC
    names(Pr_R) <- paste0("R", 1:4)
  }

  return(Pr_R)
}
