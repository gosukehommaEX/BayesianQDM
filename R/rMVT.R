#' Sample from a Bivariate Non-Standardized Multivariate t-Distribution
#'
#' Generates random samples from a bivariate non-standardized multivariate
#' t-distribution \eqn{t_\nu(\mu, V)}, defined as the distribution of
#' \eqn{\mu + V^{1/2} Z / \sqrt{W / \nu}}, where \eqn{Z \sim N_2(0, I_2)} and
#' \eqn{W \sim \chi^2(\nu)} are independent. This parameterisation corresponds
#' to the marginal posterior distribution of the mean vector under both the
#' vague (Jeffreys) prior and the Normal-Inverse-Wishart (NIW) conjugate prior
#' for bivariate normal data (see Details).
#'
#' @param n A positive integer. Number of random samples to draw.
#' @param df A positive numeric value. Degrees of freedom \eqn{\nu > 0}.
#' @param mu A numeric vector of length 2. Location parameter
#'        \eqn{\mu = (\mu_1, \mu_2)^T}.
#' @param V A 2x2 positive-definite numeric matrix. Scale matrix. The
#'        covariance of the distribution is \eqn{V \nu / (\nu - 2)} for
#'        \eqn{\nu > 2}.
#'
#' @return A numeric matrix of dimension \code{n x 2}, where each row is one
#'         independent draw from \eqn{t_\nu(\mu, V)}.
#'
#' @details
#' The sampling algorithm uses the normal-chi-squared mixture representation:
#' \enumerate{
#'   \item Draw \eqn{Z^{(u)} \sim N_2(0, I_2)} for \eqn{u = 1, \ldots, n}.
#'   \item Draw \eqn{W^{(u)} \sim \chi^2(\nu)} independently.
#'   \item Return \eqn{\mu + R^T Z^{(u)} / \sqrt{W^{(u)} / \nu}}, where
#'         \eqn{R} is the upper Cholesky factor of \eqn{V} (i.e. \eqn{V = R^T R}).
#' }
#'
#' **Connection to posterior distributions:**
#'
#' Under the vague (Jeffreys) prior \eqn{P(\mu_k, \Sigma_k) \propto |\Sigma_k|^{-3/2}},
#' the marginal posterior of \eqn{\mu_k} is
#' \deqn{\mu_k | Y_k \sim t_{n_k - 2}\!\left(\bar{y}_k,\; \frac{S_k}{n_k(n_k - 2)}\right)}
#' so \eqn{\nu = n_k - 2} and \eqn{V = S_k / [n_k(n_k - 2)]}.
#'
#' Under the NIW prior \eqn{\mu_k | \Sigma_k \sim N_2(\mu_{0k}, \Sigma_k / \kappa_{0k})},
#' \eqn{\Sigma_k \sim \mathrm{Inv\text{-}Wishart}_{\nu_{0k}}(\Lambda_{0k}^{-1})},
#' the marginal posterior of \eqn{\mu_k} is
#' \deqn{\mu_k | Y_k \sim t_{\nu_{nk} - 1}\!\left(\mu_{nk},\; \frac{\Lambda_{nk}}{\kappa_{nk}(\nu_{nk} - 1)}\right)}
#' with updated hyperparameters
#' \deqn{\mu_{nk} = \frac{\kappa_{0k} \mu_{0k} + n_k \bar{y}_k}{\kappa_{nk}},\quad
#'       \kappa_{nk} = \kappa_{0k} + n_k,\quad \nu_{nk} = \nu_{0k} + n_k,}
#' \deqn{\Lambda_{nk} = \Lambda_{0k} + S_k + \frac{\kappa_{0k} n_k}{\kappa_{0k} + n_k}
#'       (\bar{y}_k - \mu_{0k})(\bar{y}_k - \mu_{0k})^T.}
#'
#' A small ridge \eqn{10^{-10} I_2} is added to \eqn{V} before Cholesky
#' decomposition to guard against numerical non-positive-definiteness.
#'
#' @examples
#' # Draw 1000 samples from t_5(mu, V) and check sample mean
#' set.seed(1)
#' mu <- c(2.0, -1.0)
#' V  <- matrix(c(4.0, 1.2, 1.2, 1.0), 2, 2)
#' samples <- rMVT(n = 1000L, df = 5, mu = mu, V = V)
#' colMeans(samples)  # should be close to c(2, -1)
#'
#' # Posterior sampling under vague prior (n = 20 observations)
#' set.seed(2)
#' n   <- 20L
#' ybar <- c(1.5, 0.8)
#' S    <- matrix(c(18.0, 3.6, 3.6, 9.0), 2, 2)   # sum-of-squares matrix
#' df_post <- n - 2L                                  # vague prior df
#' V_post  <- S / (n * (n - 2L))                     # vague prior scale
#' mu_post <- rMVT(n = 500L, df = df_post, mu = ybar, V = V_post)
#' colMeans(mu_post)  # should be close to ybar = c(1.5, 0.8)
#'
#' @importFrom stats rnorm rchisq
#' @export
rMVT <- function(n, df, mu, V) {

  # --- Input validation ---
  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      n != floor(n) || n < 1L)
    stop("'n' must be a single positive integer")

  if (!is.numeric(df) || length(df) != 1L || is.na(df) || df <= 0)
    stop("'df' must be a single positive numeric value")

  if (!is.numeric(mu) || length(mu) != 2L || any(is.na(mu)))
    stop("'mu' must be a numeric vector of length 2")

  if (!is.matrix(V) || !is.numeric(V) || nrow(V) != 2L || ncol(V) != 2L ||
      any(is.na(V)))
    stop("'V' must be a 2x2 numeric matrix")

  n <- as.integer(n)

  # --- Cholesky factor of V (with ridge for numerical stability) ---
  R <- tryCatch(
    chol(V + diag(1e-10, 2L)),
    error = function(e)
      stop("'V' is not positive-definite: ", conditionMessage(e))
  )

  # --- Normal-chi-squared mixture ---
  # Z[u, ] ~ N_2(0, I_2), independently for each u
  Z <- matrix(rnorm(n * 2L), nrow = n, ncol = 2L)

  # W[u] ~ chi^2(df), independently for each u
  W <- rchisq(n, df = df)

  # Scale: divide each row of Z %*% R by sqrt(W / df)
  # Z %*% R produces N_2(0, V) draws; dividing by sqrt(W/df) inflates to t_df(0, V)
  scale_factor <- sqrt(W / df)   # length n
  samples <- Z %*% R / scale_factor + matrix(mu, nrow = n, ncol = 2L, byrow = TRUE)

  samples
}
