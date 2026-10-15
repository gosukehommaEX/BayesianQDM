#' Multinomial probabilities of count vectors
#'
#' @param x Integer matrix of counts with one row per outcome; all rows have
#'   the same total.
#' @param p Numeric vector of cell probabilities (summing to 1); a vector of
#'   missing values gives missing probabilities.
#'
#' @return A numeric vector with one probability per row of \code{x}. Cells
#'   with probability 0 contribute a factor 1 if their count is 0 and make the
#'   probability 0 otherwise.
#' @noRd
dmultinom_counts <- function(x, p) {
  if (anyNA(p)) {
    return(rep(NA_real_, nrow(x)))
  }
  n <- sum(x[1L, ])
  terms <- x * rep(log(p), each = nrow(x))
  terms[x == 0L] <- 0
  exp(lgamma(n + 1) - rowSums(lgamma(x + 1)) + rowSums(terms))
}
