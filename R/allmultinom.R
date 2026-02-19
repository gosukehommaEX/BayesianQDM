#' Enumerate All Multinomial Count Vectors for a Bivariate Binary Outcome
#'
#' Generates the complete set of non-negative integer vectors
#' \eqn{(x_{00}, x_{01}, x_{10}, x_{11})} that satisfy
#' \eqn{x_{00} + x_{01} + x_{10} + x_{11} = n}, where \eqn{n} is the arm
#' sample size.  Each row of the returned matrix corresponds to one possible
#' realisation of the aggregated bivariate binary counts defined in
#' Appendix A1.1 of the reference.  This enumeration is a prerequisite for
#' exact scenario assessment (operating characteristics) in
#' \code{pbayesdecisionprob2bin}.
#'
#' @param n A single non-negative integer giving the arm sample size.
#'
#' @return An integer matrix with \eqn{\binom{n+3}{3}} rows and 4 columns
#'         named \code{x00}, \code{x01}, \code{x10}, \code{x11}.
#'         Each row is a distinct non-negative integer solution to
#'         \eqn{x_{00} + x_{01} + x_{10} + x_{11} = n}.
#'         Rows are ordered by \code{x00} (ascending), then \code{x10}
#'         (ascending), then \code{x01} (ascending).
#'
#' @details
#' The number of non-negative integer solutions to
#' \eqn{x_{00} + x_{01} + x_{10} + x_{11} = n} is the stars-and-bars count
#' \deqn{\binom{n+3}{3} = \frac{(n+1)(n+2)(n+3)}{6}.}
#' For example, \eqn{n = 10} yields 286 rows and \eqn{n = 20} yields 1771
#' rows.  The matrix is pre-allocated to avoid repeated memory reallocation,
#' making the function efficient for the sample sizes typical in rare-disease
#' proof-of-concept studies.
#'
#' This function is an internal computational building block used by
#' \code{pbayesdecisionprob2bin} to enumerate the sample space over which multinomial
#' probabilities and decision indicators are summed when computing exact
#' operating characteristics.
#'
#' @examples
#' # Example 1: n = 2 (smallest non-trivial case)
#' allmultinom(2)
#'
#' # Example 2: n = 10 (typical rare-disease PoC arm size)
#' mat <- allmultinom(10)
#' nrow(mat)          # Should be choose(13, 3) = 286
#' all(rowSums(mat) == 10)  # Every row must sum to n
#'
#' # Example 3: n = 0 (edge case - only the all-zero row)
#' allmultinom(0)
#'
#' # Example 4: Verify column names
#' colnames(allmultinom(5))
#'
#' # Example 5: Row counts match the stars-and-bars formula
#' n <- 15L
#' mat <- allmultinom(n)
#' nrow(mat) == choose(n + 3L, 3L)  # Should be TRUE
#'
#' @export
allmultinom <- function(n) {

  # --- Input validation ---
  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      n != floor(n) || n < 0L) {
    stop("'n' must be a single non-negative integer")
  }

  n <- as.integer(n)

  # --- Pre-allocate output matrix ---
  # Number of rows = C(n+3, 3) = (n+1)(n+2)(n+3) / 6
  n_rows <- as.integer((n + 1L) * (n + 2L) * (n + 3L) / 6L)
  out    <- matrix(0L, nrow = n_rows, ncol = 4L)
  colnames(out) <- c("x00", "x01", "x10", "x11")

  # --- Enumerate all solutions via nested loops ---
  # Outer loop: x00 from 0 to n
  # Middle loop: x10 from 0 to (n - x00)
  # Inner loop: x01 from 0 to (n - x00 - x10)
  # x11 is the remainder: n - x00 - x10 - x01
  k <- 0L
  for (x00 in 0L:n) {
    rem1 <- n - x00
    for (x10 in 0L:rem1) {
      rem2 <- rem1 - x10

      # Fill a block of (rem2 + 1) rows at once using vectorised assignment
      m   <- rem2 + 1L
      idx <- (k + 1L):(k + m)

      out[idx, 1L] <- x00           # x00 constant within block
      out[idx, 2L] <- 0L:rem2       # x01 varies from 0 to rem2
      out[idx, 3L] <- x10           # x10 constant within block
      out[idx, 4L] <- rem2 - 0L:rem2  # x11 = rem2 - x01

      k <- k + m
    }
  }

  return(out)
}
