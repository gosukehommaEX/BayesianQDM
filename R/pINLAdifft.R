#' Cumulative Distribution Function of the Difference of Two t-Distributed Variables by INLA
#'
#' This function calculates the cumulative distribution function (CDF) of the difference
#' between two independent t-distributed random variables using the Integrated Nested Laplace
#' Approximation (INLA) with power prior for external data. Specifically, it computes
#' P(T1 - T2 > q) where T1 and T2 follow t-distributions estimated from current and external data.
#'
#' @param nINLAsample A positive integer representing the number of iterations for INLA sampling.
#' @param q A numeric value representing the quantile threshold.
#' @param mu.n1 A numeric value representing the mean value of normal distribution for group 1 in PoC and external trials.
#' @param mu.n2 A numeric value representing the mean value of normal distribution for group 2 in PoC and external trials.
#' @param sd.n1 A positive numeric value representing the standard deviation of normal distribution for group 1 in PoC and external trials.
#' @param sd.n2 A positive numeric value representing the standard deviation of normal distribution for group 2 in PoC and external trials.
#' @param n1 A positive integer representing the sample size for group 1 in PoC trial.
#' @param n2 A positive integer representing the sample size for group 2 in PoC trial.
#' @param ne1 A positive integer representing the sample size for group 1 in external trial (can be NULL if no external treatment data).
#' @param ne2 A positive integer representing the sample size for group 2 in external trial (can be NULL if no external control data).
#' @param alpha01 A positive numeric value representing the scale parameter of the power prior for group 1 (can be NULL if no external treatment data).
#' @param alpha02 A positive numeric value representing the scale parameter of the power prior for group 2 (can be NULL if no external control data).
#'
#' @return A numeric value representing P(T1 - T2 > q), the probability that the difference
#'         between the two t-distributed variables exceeds the quantile q.
#'
#' @details
#' This function uses the Integrated Nested Laplace Approximation (INLA) to estimate posterior
#' distributions when external data are available. The method involves:
#' \itemize{
#'   \item Incorporating external data through power priors with scaling parameters α01 and α02
#'   \item Using INLA for fast approximate Bayesian inference
#'   \item Generating posterior samples to compute the probability of interest
#'   \item Supporting flexible combinations of external treatment and/or control data
#' }
#'
#' The power prior approach allows borrowing strength from historical data while controlling
#' the degree of borrowing through the power parameters. INLA provides computationally efficient
#' approximations compared to full MCMC approaches.
#'
#' @examples
#' \dontrun{
#' # Both external treatment and control data
#' pINLAdifft(nINLAsample = 1e+4, q = 4, mu.n1 = 5, mu.n2 = 0, sd.n1 = 1, sd.n2 = 1,
#'           n1 = 12, n2 = 12, ne1 = 24, ne2 = 24, alpha01 = 0.5, alpha02 = 0.5)
#'
#' # External control data only
#' pINLAdifft(nINLAsample = 1e+4, q = 4, mu.n1 = 5, mu.n2 = 0, sd.n1 = 1, sd.n2 = 1,
#'           n1 = 12, n2 = 12, ne1 = NULL, ne2 = 24, alpha01 = NULL, alpha02 = 0.5)
#' }
#'
#' @importFrom stats rnorm
#' @importFrom INLA inla inla.posterior.sample
#' @export
pINLAdifft <- function(nINLAsample, q, mu.n1, mu.n2, sd.n1, sd.n2, n1, n2, ne1, ne2, alpha01, alpha02) {
  # Check parameter sets for external data availability
  if((is.null(ne1) | is.null(alpha01)) & (is.null(ne2) | is.null(alpha02))) {
    stop('At least one complete pair (ne1 & alpha01) or (ne2 & alpha02) must be provided for external data analysis')
  }

  # Generate PoC trial data from the specified distributions
  y1 <- rnorm(n1, mu.n1, sd.n1)
  y2 <- rnorm(n2, mu.n2, sd.n2)

  if((!is.null(ne1) & !is.null(alpha01)) & (!is.null(ne2) & !is.null(alpha02))) {
    ## Both external treatment and control data are available
    # Generate external data
    ye1 <- rnorm(ne1, mu.n1, sd.n1)
    ye2 <- rnorm(ne2, mu.n2, sd.n2)
    # Combine current and external data with appropriate weights
    data <- data.frame(
      y = c(y1, y2, ye1, ye2),
      weight = c(rep(1, n2 + n1), rep(c(alpha01, alpha02), c(ne1, ne2))),
      group = factor(rep(c(1, 2, 1, 2), c(n1, n2, ne1, ne2)), levels = c(2, 1))
    )
  } else if((!is.null(ne1) & !is.null(alpha01)) & (is.null(ne2) | is.null(alpha02))) {
    ## External treatment data only
    # Generate external treatment data
    ye1 <- rnorm(ne1, mu.n1, sd.n1)
    # Combine current and external treatment data
    data <- data.frame(
      y = c(y1, y2, ye1),
      weight = c(rep(1, n2 + n1), rep(alpha01, ne1)),
      group = factor(rep(c(1, 2, 1), c(n1, n2, ne1)), levels = c(2, 1))
    )
  } else if((is.null(ne1) | is.null(alpha01)) & (!is.null(ne2) & !is.null(alpha02))) {
    ## External control data only
    # Generate external control data
    ye2 <- rnorm(ne2, mu.n2, sd.n2)
    # Combine current and external control data
    data <- data.frame(
      y = c(y1, y2, ye2),
      weight = c(rep(1, n2 + n1), rep(alpha02, ne2)),
      group = factor(rep(c(1, 2, 2), c(n1, n2, ne2)), levels = c(2, 1))
    )
  }

  # Fit Bayesian model using INLA with power prior weighting
  result <- INLA::inla(
    formula = y ~ group,
    data = data,
    family = 'gaussian',
    weights = data$weight,
    control.family = list(hyper = list(prec = list(initial = 0, fixed = FALSE))),
    control.compute = list(config = TRUE)
  )

  # Generate posterior samples from the fitted model
  samples <- INLA::inla.posterior.sample(nINLAsample, result)

  # Calculate the probability that the treatment effect exceeds the threshold q
  results <- mean(sapply(samples, function(x) x$latent['group1:1', 1]) > q)

  return(results)
}
