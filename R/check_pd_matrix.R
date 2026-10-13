#' Validate a symmetric positive definite matrix
#'
#' Checks that an argument is a finite numeric square matrix of the given
#' dimension that is symmetric and positive definite.
#'
#' @param x Object to check.
#' @param name Name of the argument, used in error messages.
#' @param dim Required number of rows and columns.
#'
#' @return \code{x}, invisibly.
#' @noRd
check_pd_matrix <- function(x, name, dim = 2L) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) != dim || ncol(x) != dim) {
    stop(sprintf("'%s' must be a numeric %d x %d matrix.", name, dim, dim),
         call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("'%s' must contain finite values only.", name),
         call. = FALSE)
  }
  if (max(abs(x - t(x))) > 1e-8 * max(1, max(abs(x)))) {
    stop(sprintf("'%s' must be symmetric.", name), call. = FALSE)
  }
  ev <- eigen((x + t(x)) / 2, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0) {
    stop(sprintf("'%s' must be positive definite.", name), call. = FALSE)
  }
  invisible(x)
}
