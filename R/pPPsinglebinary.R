#' Calculate Bayesian Posterior Probability or Posterior Predictive Probability
#' for a Clinical Trial with a Single Binary Endpoint
#'
#' This function computes Bayesian posterior probability or posterior predictive
#' probability for binary outcome clinical trials. The function supports controlled,
#' uncontrolled, and external control designs, using beta-binomial conjugate priors.
#' External data can be incorporated through power priors. This function supports
#' vector inputs for y1 and y2 parameters for efficient batch processing.
#'
#' @param prob A character string specifying the type of probability to calculate.
#'        Options are \code{'posterior'} for posterior probability or \code{'predictive'}
#'        for posterior predictive probability.
#' @param design A character string specifying the type of trial design. Options are
#'        \code{'controlled'} for randomized controlled trials, \code{'uncontrolled'}
#'        for single-arm studies, or \code{'external'} for designs incorporating external data.
#' @param theta0 A numeric value representing the pre-specified threshold value for
#'        the treatment effect (difference in response rates).
#' @param n1 A positive integer representing the number of patients in group 1
#'        (treatment) for the proof-of-concept (PoC) trial.
#' @param n2 A positive integer representing the number of patients in group 2.
#'        For \code{design = 'controlled'} or \code{'external'}: sample size of the
#'        control group. For \code{design = 'uncontrolled'}: hypothetical control
#'        sample size (typically set equal to n1).
#' @param y1 A non-negative integer or vector of integers representing the observed
#'        number of responders in group 1 for the PoC trial (must satisfy 0 <= y1 <= n1).
#' @param y2 A non-negative integer or vector of integers representing the number
#'        of responders in group 2. For \code{design = 'controlled'} or \code{'external'}:
#'        observed number of responders in the control group (must satisfy 0 <= y2 <= n2).
#'        For \code{design = 'uncontrolled'}: set to NULL and use z parameter instead.
#'        If vectors, must have the same length as y1.
#' @param a1 A positive numeric value representing the first shape parameter (alpha)
#'        of the prior beta distribution for group 1.
#' @param a2 A positive numeric value representing the first shape parameter (alpha)
#'        of the prior beta distribution for group 2 (control or hypothetical control).
#' @param b1 A positive numeric value representing the second shape parameter (beta)
#'        of the prior beta distribution for group 1.
#' @param b2 A positive numeric value representing the second shape parameter (beta)
#'        of the prior beta distribution for group 2 (control or hypothetical control).
#' @param m1 A positive integer representing the number of patients in group 1 for
#'        the future trial (required if \code{prob = 'predictive'}, otherwise set to NULL).
#' @param m2 A positive integer representing the number of patients in group 2 for
#'        the future trial (required if \code{prob = 'predictive'}, otherwise set to NULL).
#' @param z A non-negative integer representing the hypothetical control responder count
#'        (required if \code{design = 'uncontrolled'}, otherwise set to NULL).
#'        This specifies the expected number of responders in a hypothetical control
#'        group based on historical data or prior knowledge. When used, y2 should be NULL.
#' @param ne1 A positive integer representing the number of patients in group 1 for
#'        the external data (required if \code{design = 'external'}, otherwise set to NULL).
#' @param ne2 A positive integer representing the number of patients in group 2 for
#'        the external data (required if \code{design = 'external'}, otherwise set to NULL).
#' @param ye1 A non-negative integer representing the observed number of responders
#'        in group 1 for the external data (required if \code{design = 'external'},
#'        otherwise set to NULL).
#' @param ye2 A non-negative integer representing the observed number of responders
#'        in group 2 for the external data (required if \code{design = 'external'},
#'        otherwise set to NULL).
#' @param ae1 A numeric value in (0, 1] representing the power prior scale parameter
#'        for group 1 (required if \code{design = 'external'}, otherwise set to NULL).
#'        Controls the degree of borrowing: 0 = no borrowing, 1 = full borrowing.
#' @param ae2 A numeric value in (0, 1] representing the power prior scale parameter
#'        for group 2 (required if \code{design = 'external'}, otherwise set to NULL).
#'        Controls the degree of borrowing: 0 = no borrowing, 1 = full borrowing.
#' @param lower.tail A logical value; if TRUE (default), probabilities are
#'        P(treatment effect <= theta0), otherwise P(treatment effect > theta0).
#'
#' @return A numeric value or vector in \code{[0, 1]} representing the Bayesian posterior
#'         probability or posterior predictive probability that the treatment effect
#'         exceeds (or is below) the threshold theta0. If y1 and y2 are vectors,
#'         returns a vector of the same length.
#'
#' @details
#' The function computes probabilities based on beta-binomial conjugate analysis:
#' \itemize{
#'   \item **Prior**: The proportion of responders (response_j) for each group (j = 1, 2)
#'         follows a beta distribution: response_j ~ Beta(aj, bj)
#'   \item **Posterior**: After observing data, the posterior distribution is also
#'         beta: response_j | yj ~ Beta(aj + yj, bj + nj - yj)
#'   \item **Posterior predictive**: The predictive distribution of future trial
#'         data follows a beta-binomial distribution
#' }
#'
#' **Design-specific handling**:
#' \itemize{
#'   \item **Controlled design**: Uses observed control data (y2, n2) directly
#'   \item **Uncontrolled design**: Uses hypothetical control value z. The calculation
#'         treats z as the expected number of responders in a hypothetical control group
#'         with sample size n2. Parameter y2 should be NULL when using uncontrolled design.
#'   \item **External design**: Incorporates external data through power priors.
#'         The effective prior becomes Beta(aj + aej times yej, bj + aej times (nej - yej)),
#'         where aej controls the degree of borrowing.
#' }
#'
#' **Probability types**:
#' \itemize{
#'   \item **Posterior probability**: P(response1 - response2 > theta0 | current data)
#'   \item **Posterior predictive probability**: P(future trial success | current data)
#' }
#'
#' @examples
#' # Example 1: Controlled design - posterior probability (scalar)
#' pPPsinglebinary(
#'   prob = 'posterior', design = 'controlled', theta0 = 0.15,
#'   n1 = 12, n2 = 15, y1 = 7, y2 = 5,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = NULL, m2 = NULL, z = NULL,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 2: Vector input for y1 and y2
#' pPPsinglebinary(
#'   prob = 'posterior', design = 'controlled', theta0 = 0.15,
#'   n1 = 12, n2 = 15, y1 = c(5, 6, 7, 8), y2 = c(3, 4, 5, 6),
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = NULL, m2 = NULL, z = NULL,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 3: Uncontrolled design - using hypothetical control with z parameter
#' pPPsinglebinary(
#'   prob = 'posterior', design = 'uncontrolled', theta0 = 0.20,
#'   n1 = 20, n2 = 20, y1 = 12, y2 = NULL,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = NULL, m2 = NULL, z = 3,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 4: Posterior predictive probability for controlled design
#' pPPsinglebinary(
#'   prob = 'predictive', design = 'controlled', theta0 = 0.10,
#'   n1 = 12, n2 = 15, y1 = 7, y2 = 5,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = 30, m2 = 30, z = NULL,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 5: External design with power prior (50 percent borrowing for control)
#' pPPsinglebinary(
#'   prob = 'posterior', design = 'external', theta0 = 0.15,
#'   n1 = 12, n2 = 15, y1 = 7, y2 = 9,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = NULL, m2 = NULL, z = NULL,
#'   ne1 = 12, ne2 = 12, ye1 = 6, ye2 = 6, ae1 = 0.5, ae2 = 0.5,
#'   lower.tail = FALSE
#' )
#'
#' @importFrom stats integrate
#' @export
pPPsinglebinary <- function(prob = 'posterior', design = 'controlled', theta0,
                            n1, n2, y1, y2, a1, a2, b1, b2,
                            m1, m2, z = NULL, ne1, ne2, ye1, ye2, ae1, ae2, lower.tail = TRUE) {
  # Validate parameter sets for posterior predictive probability
  if((prob == 'predictive') & (sum(sapply(list(m1, m2), is.null)) > 0)) {
    stop('If you calculate posterior predictive probability, m1 and m2 must be non-null')
  }

  # Validate parameter sets for external design
  if((design == 'external') & (sum(sapply(list(ne1, ne2, ye1, ye2, ae1, ae2), is.null)) > 0)) {
    stop('If you use external data, ne1, ne2, ye1, ye2, ae1, and ae2 must be non-null')
  }

  # Validate parameter sets for uncontrolled design
  if((design == 'uncontrolled') & (is.null(z))) {
    stop('If you consider uncontrolled design, z (hypothetical control responder count) must be non-null')
  }

  # Ensure y1 is a numeric vector
  y1 <- as.numeric(y1)

  # Resolve y2: for uncontrolled design, replicate z to match length of y1
  if(is.null(y2)) {
    if(design == 'uncontrolled' & !is.null(z)) {
      y2 <- rep(z, length(y1))
    } else {
      stop('y2 must be provided or z must be non-null for uncontrolled design')
    }
  }
  y2 <- as.numeric(y2)

  # Validate that y1 and y2 have the same length
  if(length(y1) != length(y2)) {
    stop('y1 and y2 must have the same length')
  }

  # Compute posterior shape parameters for all (y1, y2) pairs simultaneously
  # (vectorized over the full input vectors - no loop required)
  if(design == 'external') {
    # Incorporate external data through power priors
    s11 <- y1 + a1 + ae1 * ye1
    s21 <- n1 - y1 + b1 + ae1 * (ne1 - ye1)
    s12 <- y2 + a2 + ae2 * ye2
    s22 <- n2 - y2 + b2 + ae2 * (ne2 - ye2)
  } else {
    # Controlled or uncontrolled design (y2 already set to z for uncontrolled)
    s11 <- y1 + a1
    s21 <- n1 - y1 + b1
    s12 <- y2 + a2
    s22 <- n2 - y2 + b2
  }

  # Calculate probabilities using mapply to apply scalar functions across all pairs
  if(prob == 'posterior') {
    # Posterior probability P(response1 - response2 > theta0 | data)
    # p2betadiff requires scalar parameters, so mapply is used
    results <- mapply(
      p2betadiff,
      alpha1 = s11, alpha2 = s12, beta1 = s21, beta2 = s22,
      MoreArgs = list(q = theta0, lower.tail = lower.tail)
    )
  } else if(prob == 'predictive') {
    # Posterior predictive probability P(future trial success | current data)
    # p2betabinomdiff requires scalar parameters, so mapply is used
    results <- mapply(
      p2betabinomdiff,
      alpha1 = s11, alpha2 = s12, beta1 = s21, beta2 = s22,
      MoreArgs = list(q = theta0, m1 = m1, m2 = m2, lower.tail = lower.tail)
    )
  }

  return(results)
}
