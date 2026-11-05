#' Print the operating characteristics of several designs
#'
#' @param x An object of class \code{"qdm_comparison"}.
#' @param digits Number of significant digits of the probabilities.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' qdm_compare(list(n12 = des, n24 = qdm_sample_size(des, 24, 24)),
#'             data.frame(pi_t = c(0.3, 0.6), pi_c = 0.3))
#' @export
print.qdm_comparison <- function(x, digits = 4L, ...) {
  method <- attr(x, "method")
  n <- nlevels(x$design)
  cat(sprintf("Operating characteristics of %d %s\n", n,
              if (n == 1L) "design" else "designs"))
  if (!is.null(method)) {
    if (length(unique(method)) == 1L) {
      cat(sprintf("Method: %s\n", method[[1L]]))
    } else {
      cat(sprintf("  %s: %s\n", names(method), method), sep = "")
    }
  }
  y <- x
  attr(y, "designs") <- NULL
  attr(y, "method") <- NULL
  class(y) <- "data.frame"
  print(y, digits = digits, row.names = FALSE)
  invisible(x)
}
