#' Warn about positive probabilities of Miss decisions
#'
#' @param miss Numeric vector of probabilities of Miss decisions.
#' @param rule Rule of the criterion (\code{"joint"}, \code{"and"} or
#'   \code{"or"}), which determines when Miss decisions are impossible.
#'
#' @return \code{NULL}, invisibly.
#' @noRd
warn_miss <- function(miss, rule = "joint") {
  k <- sum(miss > 0, na.rm = TRUE)
  if (k > 0L) {
    hint <- switch(
      rule,
      joint = "A Miss decision is impossible if gamma_go + gamma_nogo > 1.",
      and = paste("A Miss decision is impossible if gamma_go + gamma_nogo >",
                  "1 for at least one endpoint."),
      or = paste("With rule = \"or\", a Miss decision can occur for any",
                 "thresholds (Go on one endpoint and NoGo on the other).")
    )
    warning(sprintf(paste("The Go and NoGo criteria are both met (Miss) with",
                          "positive probability in %d scenario(s); see the",
                          "column 'Miss'. %s"), k, hint), call. = FALSE)
  }
  invisible(NULL)
}
