#' Warn about positive probabilities of Miss decisions
#'
#' @param miss Numeric vector of probabilities of Miss decisions.
#'
#' @return \code{NULL}, invisibly.
#' @noRd
warn_miss <- function(miss) {
  k <- sum(miss > 0)
  if (k > 0L) {
    warning(sprintf(paste("The Go and NoGo criteria are both met (Miss) with",
                          "positive probability in %d scenario(s); see the",
                          "column 'Miss'. A Miss decision is impossible if",
                          "gamma_go + gamma_nogo > 1."), k), call. = FALSE)
  }
  invisible(NULL)
}
