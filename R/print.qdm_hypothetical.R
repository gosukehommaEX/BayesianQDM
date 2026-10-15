#' Print a hypothetical control specification
#'
#' @param x An object of class \code{"qdm_hypothetical"}.
#' @param ... Passed to \code{\link{format.qdm_hypothetical}}.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' print(hypothetical(data = data_binary(n = 15, y = 5)))
#' @export
print.qdm_hypothetical <- function(x, ...) {
  cat("Hypothetical control:", format(x, ...), "\n")
  invisible(x)
}
