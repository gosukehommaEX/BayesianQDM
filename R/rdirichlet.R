#' Generate Random Samples from a Dirichlet Distribution
#'
#' Generates random samples from a Dirichlet distribution using the Gamma
#' representation: if \eqn{Y_i \sim \mathrm{Gamma}(\alpha_i, 1)} independently
#' for \eqn{i = 1, \ldots, K}, then
#' \eqn{(Y_1 / S, \ldots, Y_K / S) \sim \mathrm{Dirichlet}(\alpha_1, \ldots, \alpha_K)},
#' where \eqn{S = \sum_{i=1}^{K} Y_i}.
#'
#' @param n A positive integer specifying the number of random vectors to
#'        generate.
#' @param alpha A numeric vector of length \eqn{K \ge 2} containing positive
#'        concentration parameters of the Dirichlet distribution.
#'        All elements must be strictly positive.
#'
#' @return A numeric matrix of dimensions \code{n x K} where each row is one
#'         random draw from the Dirichlet distribution, with all elements in
#'         \code{[0, 1]} and each row summing to 1.
#'         When \code{n = 1}, a numeric vector of length \code{K} is returned.
#'
#' @details
#' The Dirichlet distribution is a multivariate generalisation of the Beta
#' distribution and is commonly used as a conjugate prior for multinomial
#' proportions in Bayesian statistics.
#'
#' The probability density function is:
#' \deqn{f(x_1, \ldots, x_K) =
#'   \frac{\Gamma\!\left(\sum_{i=1}^{K} \alpha_i\right)}
#'        {\prod_{i=1}^{K} \Gamma(\alpha_i)}
#'   \prod_{i=1}^{K} x_i^{\alpha_i - 1}}
#' where \eqn{x_i > 0} and \eqn{\sum_{i=1}^{K} x_i = 1}.
#'
#' Key properties:
#' \itemize{
#'   \item Each marginal follows a Beta distribution:
#'         \eqn{X_i \sim \mathrm{Beta}\!\left(\alpha_i,\,
#'         \sum_{l \neq i} \alpha_l\right)}.
#'   \item \eqn{E[X_i] = \alpha_i / \sum_{l=1}^{K} \alpha_l}.
#'   \item Components are negatively correlated unless one component dominates.
#' }
#'
#' Implementation steps:
#' \enumerate{
#'   \item Generate independent \eqn{Y_i \sim \mathrm{Gamma}(\alpha_i, 1)} for
#'         each \eqn{i = 1, \ldots, K}.
#'   \item Normalise: \eqn{X_i = Y_i / \sum_{l=1}^{K} Y_l}.
#' }
#'
#' @examples
#' # Example 1: Generate 5 samples from Dirichlet(1, 1, 1) - uniform on simplex
#' samples <- rdirichlet(5, c(1, 1, 1))
#' print(samples)
#' rowSums(samples)  # Each row should sum to 1
#'
#' # Example 2: Generate samples with unequal concentrations
#' samples <- rdirichlet(1000, c(2, 5, 3))
#' colMeans(samples)  # Expected values: approximately c(0.2, 0.5, 0.3)
#'
#' # Example 3: Sparse Dirichlet (small alpha values)
#' samples <- rdirichlet(100, c(0.1, 0.1, 0.1, 0.1))
#' head(samples)  # Most weight concentrated on one component
#'
#' # Example 4: Concentrated Dirichlet (large alpha values)
#' samples <- rdirichlet(100, c(100, 100, 100))
#' colMeans(samples)  # Concentrated around c(1/3, 1/3, 1/3)
#'
#' # Example 5: Bayesian update with Jeffreys prior for 4 categories
#' prior_alpha     <- c(0.5, 0.5, 0.5, 0.5)
#' observed_counts <- c(10, 5, 8, 7)
#' posterior_samples <- rdirichlet(1000, prior_alpha + observed_counts)
#' colMeans(posterior_samples)  # Posterior mean
#'
#' @importFrom stats rgamma
#' @export
rdirichlet <- function(n, alpha) {

  # --- Input validation ---
  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      n != floor(n) || n < 1L) {
    stop("'n' must be a single positive integer")
  }

  if (!is.numeric(alpha) || length(alpha) < 2L || any(is.na(alpha)) ||
      any(alpha <= 0)) {
    stop("'alpha' must be a numeric vector of length >= 2 with all elements strictly positive")
  }

  # --- Generate Dirichlet samples via the Gamma representation ---
  # Build an (n x k) matrix where column j contains n independent
  # Gamma(alpha[j], 1) draws. Using rep(alpha, each = n) ensures that
  # the j-th block of n values uses shape alpha[j]; filling column-wise
  # (byrow = FALSE, the R default) aligns each block with column j.
  k <- length(alpha)
  gamma_samples <- matrix(
    rgamma(n * k, shape = rep(alpha, each = n), rate = 1),
    nrow = n,
    ncol = k
  )

  # Normalise each row so that it sums to 1
  dirichlet_samples <- gamma_samples / rowSums(gamma_samples)

  # Return a plain vector when only one sample is requested
  if (n == 1L) {
    dirichlet_samples <- as.numeric(dirichlet_samples)
  }

  return(dirichlet_samples)
}
