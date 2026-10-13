#' Random generation from the Dirichlet distribution
#'
#' Draws random vectors from the Dirichlet distribution by normalizing
#' independent Gamma variables, using the random number generator of R.
#'
#' @param n Number of draws (a whole number of at least 1).
#' @param alpha Numeric vector of at least two positive concentration
#'   parameters.
#'
#' @return A numeric matrix with \code{n} rows and \code{length(alpha)}
#'   columns whose rows sum to 1 (a matrix is returned even when
#'   \code{n = 1}). Column names are taken from \code{names(alpha)}.
#'
#' @family distributions
#' @examples
#' set.seed(1)
#' rdirichlet(3, c(1.5, 0.5, 2, 1))
#' @export
rdirichlet <- function(n, alpha) {
  check_numeric(n, "n", lower = 1, integer = TRUE)
  check_numeric(alpha, "alpha", len = NULL, lower = 0, lower_open = TRUE)
  if (length(alpha) < 2L) {
    stop("'alpha' must have length 2 or more.", call. = FALSE)
  }
  n <- as.integer(round(n))
  k <- length(alpha)
  # Column j holds n Gamma(alpha[j], 1) draws
  g <- matrix(stats::rgamma(n * k, shape = rep(alpha, each = n)), nrow = n,
              ncol = k)
  out <- g / rowSums(g)
  colnames(out) <- names(alpha)
  out
}
