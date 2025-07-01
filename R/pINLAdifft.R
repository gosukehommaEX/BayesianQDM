#' Cumulative distribution function (CDF) of the difference of two t-distributed outcomes by INLA.
#'
#' This function calculates CDF(s) of difference between two variables following t-distribution
#' using the integrated nested Laplace approximation (INLA) with power prior for external data.
#'
#' @param nINLAsample A number of iterations for INLA sampling.
#' @param q Quantile value.
#' @param mu.n1 Mean value of normal distribution for group 1 in PoC and external trials.
#' @param mu.n2 Mean value of normal distribution for group 2 in PoC and external trials.
#' @param sd.n1 Standard deviation value of normal distribution for group 1 in PoC and external trials.
#' @param sd.n2 Standard deviation value of normal distribution for group 2 in PoC and external trials.
#' @param n1 Sample size for group 1 in PoC trial.
#' @param n2 Sample size for group 2 in PoC trial.
#' @param ne1 Sample size for group 1 in external trial (can be NULL if no external treatment data).
#' @param ne2 Sample size for group 2 in external trial (can be NULL if no external control data).
#' @param alpha01 Scale parameter of the power prior for group 1 (can be NULL if no external treatment data).
#' @param alpha02 Scale parameter of the power prior for group 2 (can be NULL if no external control data).
#'
#' @return Cumulative probability(s)
#'
#' @examples
#' # Both external treatment and control
#' pINLAdifft(nINLAsample = 1e+4, q = 4, mu.n1 = 5, mu.n2 = 0, sd.n1 = 1, sd.n2 = 1,
#'           n1 = 12, n2 = 12, ne1 = 24, ne2 = 24, alpha01 = 0.5, alpha02 = 0.5)
#'
#' # External control only
#' pINLAdifft(nINLAsample = 1e+4, q = 4, mu.n1 = 5, mu.n2 = 0, sd.n1 = 1, sd.n2 = 1,
#'           n1 = 12, n2 = 12, ne1 = NULL, ne2 = 24, alpha01 = NULL, alpha02 = 0.5)
#'
#' @importFrom stats rnorm
#' @importFrom INLA inla inla.posterior.sample
#' @export
pINLAdifft = function(nINLAsample, q, mu.n1, mu.n2, sd.n1, sd.n2, n1, n2, ne1, ne2, alpha01, alpha02) {
  # Check parameter sets
  if((is.null(ne1) | is.null(alpha01)) & (is.null(ne2) | is.null(alpha02))) {
    stop('At least one complete pair (ne1 & alpha01) or (ne2 & alpha02) must be provided for external data analysis')
  }

  # Generate PoC trial data
  y1 <- rnorm(n1, mu.n1, sd.n1)
  y2 <- rnorm(n2, mu.n2, sd.n2)

  if((!is.null(ne1) & !is.null(alpha01)) & (!is.null(ne2) & !is.null(alpha02))) {
    ## Both external treatment and control
    # External data
    ye1 <- rnorm(ne1, mu.n1, sd.n1)
    ye2 <- rnorm(ne2, mu.n2, sd.n2)
    # Dataset
    data <- data.frame(
      y = c(y1, y2, ye1, ye2),
      weight = c(rep(1, n2 + n1), rep(c(alpha01, alpha02), c(ne1, ne2))),
      group = factor(rep(c(1, 2, 1, 2), c(n1, n2, ne1, ne2)), levels = c(2, 1))
    )
  } else if((!is.null(ne1) & !is.null(alpha01)) & (is.null(ne2) | is.null(alpha02))) {
    ## External treatment only
    # External data
    ye1 <- rnorm(ne1, mu.n1, sd.n1)
    # Dataset
    data <- data.frame(
      y = c(y1, y2, ye1),
      weight = c(rep(1, n2 + n1), rep(alpha01, ne1)),
      group = factor(rep(c(1, 2, 1), c(n1, n2, ne1)), levels = c(2, 1))
    )
  } else if((is.null(ne1) | is.null(alpha01)) & (!is.null(ne2) & !is.null(alpha02))) {
    ## External control only
    # External data
    ye2 <- rnorm(ne2, mu.n2, sd.n2)
    # Dataset
    data <- data.frame(
      y = c(y1, y2, ye2),
      weight = c(rep(1, n2 + n1), rep(alpha02, ne2)),
      group = factor(rep(c(1, 2, 2), c(n1, n2, ne2)), levels = c(2, 1))
    )
  }
  # Get posterior of treatment effect (group coefficient)
  result <- inla(
    formula = y ~ group,
    data = data,
    family = 'gaussian',
    weights = weight,
    control.family = list(hyper = list(prec = list(initial = 0, fixed = FALSE))),
    control.compute = list(config = TRUE)
  )
  samples <- inla.posterior.sample(nINLAsample, result)
  results <- mean(sapply(samples, function(x) x$latent['group1:1', 1]) > q)
  return(results)
}
