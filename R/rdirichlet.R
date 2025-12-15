#' Generate Random Samples from Dirichlet Distribution
#'
#' This function generates random samples from a Dirichlet distribution using
#' the property that if X_i ~ Gamma(α_i, 1) independently, then
#' (X_1/S, ..., X_k/S) ~ Dirichlet(α_1, ..., α_k), where S = sum(X_i).
#'
#' @param n A positive integer representing the number of random vectors to generate.
#' @param alpha A numeric vector of positive shape parameters (α_1, ..., α_k) for
#'        the Dirichlet distribution. All elements must be positive.
#'
#' @return A numeric matrix of size \code{n × length(alpha)}, where each row is
#'         a random sample from the Dirichlet distribution. Each row sums to 1,
#'         and all elements are in the interval \code{[0, 1]}.
#'
#' @details
#' The Dirichlet distribution is a multivariate generalization of the Beta distribution
#' and is commonly used as a prior distribution in Bayesian statistics, particularly
#' for multinomial proportions.
#'
#' The probability density function is:
#' \deqn{f(x_1, ..., x_k) = \frac{\Gamma(\sum α_i)}{\prod \Gamma(α_i)} \prod x_i^{α_i - 1}}
#'
#' where \eqn{x_i > 0}, \eqn{\sum x_i = 1}, and \eqn{α_i > 0}.
#'
#' **Properties**:
#' \itemize{
#'   \item Each marginal distribution is Beta: \eqn{X_i ~ Beta(α_i, \sum_{j≠i} α_j)}
#'   \item \eqn{E[X_i] = α_i / \sum α_j}
#'   \item Components are negatively correlated (unless one component dominates)
#' }
#'
#' **Implementation**: The function uses the Gamma representation method:
#' 1. Generate independent \eqn{Y_i ~ Gamma(α_i, 1)}
#' 2. Normalize: \eqn{X_i = Y_i / \sum Y_j}
#'
#' @examples
#' # Example 1: Generate 5 samples from Dirichlet(1, 1, 1) - uniform on simplex
#' samples <- rdirichlet(5, c(1, 1, 1))
#' print(samples)
#' # Check that rows sum to 1
#' rowSums(samples)
#'
#' # Example 2: Generate samples with unequal concentrations
#' samples <- rdirichlet(1000, c(2, 5, 3))
#' # Expected values should be approximately c(0.2, 0.5, 0.3)
#' colMeans(samples)
#'
#' # Example 3: Sparse Dirichlet (small alpha values)
#' samples <- rdirichlet(100, c(0.1, 0.1, 0.1, 0.1))
#' # Should produce sparse vectors (most weight on one component)
#' head(samples)
#'
#' # Example 4: Concentrated Dirichlet (large alpha values)
#' samples <- rdirichlet(100, c(100, 100, 100))
#' # Should produce samples concentrated around (1/3, 1/3, 1/3)
#' colMeans(samples)
#'
#' # Example 5: Use for multinomial prior in Bayesian inference
#' # Prior: Dirichlet(0.5, 0.5, 0.5, 0.5) - Jeffreys prior for 4 categories
#' prior_alpha <- c(0.5, 0.5, 0.5, 0.5)
#' # Observed data: counts in 4 categories
#' observed_counts <- c(10, 5, 8, 7)
#' # Posterior: Dirichlet(0.5+10, 0.5+5, 0.5+8, 0.5+7)
#' posterior_samples <- rdirichlet(1000, prior_alpha + observed_counts)
#' # Posterior mean
#' colMeans(posterior_samples)
#'
#' @importFrom stats rgamma
#' @export
rdirichlet <- function(n, alpha) {
  # Input validation
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer")
  }

  if (!is.numeric(alpha) || any(alpha <= 0)) {
    stop("alpha must be a numeric vector with all positive values")
  }

  if (length(alpha) < 2) {
    stop("alpha must have at least 2 elements")
  }

  # Generate random numbers from Gamma distribution
  # Each row corresponds to one Dirichlet sample
  k <- length(alpha)
  gamma_samples <- matrix(
    rgamma(n * k, shape = alpha, rate = 1),
    nrow = n,
    byrow = TRUE
  )

  # Normalize each row to sum to 1
  row_sums <- rowSums(gamma_samples)
  dirichlet_samples <- gamma_samples / row_sums

  return(dirichlet_samples)
}
