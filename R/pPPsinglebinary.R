#' Calculate Bayesian Posterior Probability or Posterior Predictive Probability
#' for a Clinical Trial with a Single Binary Endpoint
#'
#' This function computes Bayesian posterior probability or posterior predictive
#' probability for binary outcome clinical trials. The function supports controlled,
#' uncontrolled, and external control designs, using beta-binomial conjugate priors.
#' External data can be incorporated through power priors.
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
#' @param y1 A non-negative integer representing the observed number of responders
#'        in group 1 for the PoC trial (must satisfy 0 ≤ y1 ≤ n1).
#' @param y2 A non-negative integer representing the number of responders in group 2.
#'        For \code{design = 'controlled'} or \code{'external'}: observed number of
#'        responders in the control group (must satisfy 0 ≤ y2 ≤ n2).
#'        For \code{design = 'uncontrolled'}: hypothetical control responder count
#'        based on historical data or prior knowledge. This is mathematically equivalent
#'        to the \code{z} parameter used in pGNGsinglebinary().
#' @param a1 A positive numeric value representing the first shape parameter (α)
#'        of the prior beta distribution for group 1.
#' @param a2 A positive numeric value representing the first shape parameter (α)
#'        of the prior beta distribution for group 2 (control or hypothetical control).
#' @param b1 A positive numeric value representing the second shape parameter (β)
#'        of the prior beta distribution for group 1.
#' @param b2 A positive numeric value representing the second shape parameter (β)
#'        of the prior beta distribution for group 2 (control or hypothetical control).
#' @param m1 A positive integer representing the number of patients in group 1 for
#'        the future trial (required if \code{prob = 'predictive'}, otherwise set to NULL).
#' @param m2 A positive integer representing the number of patients in group 2 for
#'        the future trial (required if \code{prob = 'predictive'}, otherwise set to NULL).
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
#'        P(treatment effect ≤ theta0), otherwise P(treatment effect > theta0).
#'
#' @return A numeric value in \code{[0, 1]} representing the Bayesian posterior probability
#'         or posterior predictive probability that the treatment effect exceeds
#'         (or is below) the threshold theta0.
#'
#' @details
#' The function computes probabilities based on beta-binomial conjugate analysis:
#' \itemize{
#'   \item **Prior**: The proportion of responders (πj) for each group (j = 1, 2)
#'         follows a beta distribution: πj ~ Beta(aj, bj)
#'   \item **Posterior**: After observing data, the posterior distribution is also
#'         beta: πj | yj ~ Beta(aj + yj, bj + nj - yj)
#'   \item **Posterior predictive**: The predictive distribution of future trial
#'         data follows a beta-binomial distribution
#' }
#'
#' **Design-specific handling**:
#' \itemize{
#'   \item **Controlled design**: Uses observed control data (y2, n2) directly
#'   \item **Uncontrolled design**: Uses hypothetical control values (y2, n2) based on
#'         historical data. The calculation is mathematically identical to controlled
#'         design, but y2 represents an expected/hypothetical value rather than observed data.
#'   \item **External design**: Incorporates external data through power priors.
#'         The effective prior becomes Beta(aj + aej × yej, bj + aej × (nej - yej)),
#'         where aej controls the degree of borrowing.
#' }
#'
#' **Probability types**:
#' \itemize{
#'   \item **Posterior probability**: P(π1 - π2 > theta0 | current data)
#'   \item **Posterior predictive probability**: P(future trial success | current data)
#' }
#'
#' @examples
#' # Example 1: Controlled design - posterior probability
#' pPPsinglebinary(
#'   prob = 'posterior', design = 'controlled', theta0 = 0.15,
#'   n1 = 12, n2 = 15, y1 = 7, y2 = 5,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = NULL, m2 = NULL,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 2: Uncontrolled design - using hypothetical control
#' # y2 = 3 represents expected responders in hypothetical control (n2 = 20)
#' pPPsinglebinary(
#'   prob = 'posterior', design = 'uncontrolled', theta0 = 0.20,
#'   n1 = 20, n2 = 20, y1 = 12, y2 = 3,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = NULL, m2 = NULL,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 3: Posterior predictive probability for controlled design
#' pPPsinglebinary(
#'   prob = 'predictive', design = 'controlled', theta0 = 0.10,
#'   n1 = 12, n2 = 15, y1 = 7, y2 = 5,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = 30, m2 = 30,
#'   ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 4: External design with power prior (50% borrowing for control)
#' pPPsinglebinary(
#'   prob = 'posterior', design = 'external', theta0 = 0.15,
#'   n1 = 12, n2 = 15, y1 = 7, y2 = 9,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = NULL, m2 = NULL,
#'   ne1 = 12, ne2 = 12, ye1 = 6, ye2 = 6, ae1 = 0.5, ae2 = 0.5,
#'   lower.tail = FALSE
#' )
#'
#' @export
pPPsinglebinary <- function(prob = 'posterior', design = 'controlled', theta0,
                            n1, n2, y1, y2, a1, a2, b1, b2,
                            m1, m2, ne1, ne2, ye1, ye2, ae1, ae2, lower.tail = TRUE) {
  # Validate parameter sets for posterior predictive probability
  if((prob == 'predictive') & (sum(sapply(list(m1, m2), is.null)) > 0)) {
    stop('If you calculate posterior predictive probability, m1 and m2 must be non-null')
  }

  # Validate parameter sets for external design
  if((design == 'external') & (sum(sapply(list(ne1, ne2, ye1, ye2, ae1, ae2), is.null)) > 0)) {
    stop('If you use external data, ne1, ne2, ye1, ye2, ae1, and ae2 must be non-null')
  }

  # Calculate posterior shape parameters for group 1 (treatment)
  # For controlled/uncontrolled: s11 = y1 + a1
  # For external: s11 = y1 + a1 + ae1 * ye1 (incorporate external treatment data)
  s11 <- y1 + a1 + (design == 'external') * c(ye1 * ae1, 0)[1]
  s21 <- n1 - y1 + b1 + (design == 'external') * c((ne1 - ye1) * ae1, 0)[1]

  # Calculate posterior shape parameters for group 2 (control or hypothetical control)
  # For controlled: y2 = observed control responders
  # For uncontrolled: y2 = hypothetical control responders (equivalent to z in pGNGsinglebinary)
  # For external: additionally incorporate external control data via power prior
  s12 <- y2 + a2 + (design == 'external') * c(ye2 * ae2, 0)[1]
  s22 <- n2 - y2 + b2 + (design == 'external') * c((ne2 - ye2) * ae2, 0)[1]

  # Calculate probability based on specified type
  if(prob == 'posterior') {
    # Calculate posterior probability P(π1 - π2 > theta0 | data)
    # using the difference of two beta distributions
    # Works identically for controlled, uncontrolled, and external designs
    results <- p2betadiff(theta0, s11, s12, s21, s22, lower.tail)
  } else if(prob == 'predictive') {
    # Calculate posterior predictive probability P(future success | current data)
    # using the difference of two beta-binomial distributions
    results <- p2betabinomdiff(theta0, m1, m2, s11, s12, s21, s22, lower.tail)
  }

  return(results)
}
