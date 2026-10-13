#' Beta-binomial probability mass function
#'
#' Computes the probability mass function of the Beta-binomial distribution
#' on the log scale to avoid overflow and underflow for large \code{m}.
#'
#' @param k Integer vector of counts from 0 to \code{m}.
#' @param m Number of trials.
#' @param a,b Positive shape parameters.
#'
#' @return A numeric vector of probabilities.
#' @noRd
dbetabinom <- function(k, m, a, b) {
  exp(lchoose(m, k) + lbeta(k + a, m - k + b) - lbeta(a, b))
}
