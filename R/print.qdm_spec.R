#' Print a design specification
#'
#' @param x An object of class \code{"qdm_spec"}, created by
#'   \code{\link{qdm_endpoint}} and the other design steps.
#' @param ... Passed to \code{\link{format.qdm_spec}}.
#'
#' @return \code{x}, invisibly.
#'
#' @examples
#' qdm_endpoint("binary") |>
#'   qdm_sample_size(n_t = 12, n_c = 12) |>
#'   print()
#' @export
print.qdm_spec <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}
