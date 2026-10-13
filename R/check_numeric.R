#' Validate a numeric argument
#'
#' Checks the type, length, finiteness, integer property and range of a
#' numeric argument, and stops with an informative message if a check fails.
#'
#' @param x Object to check.
#' @param name Name of the argument, used in error messages.
#' @param len Required length, or \code{NULL} to accept any positive length.
#' @param lower,upper Lower and upper bounds for all elements.
#' @param lower_open,upper_open Logical; if \code{TRUE}, the bound itself is
#'   not allowed.
#' @param integer Logical; if \code{TRUE}, all elements must be whole numbers.
#' @param allow_inf Logical; if \code{TRUE}, infinite values are accepted.
#'
#' @return \code{x}, invisibly.
#' @noRd
check_numeric <- function(x, name, len = 1L, lower = -Inf, upper = Inf,
                          lower_open = FALSE, upper_open = FALSE,
                          integer = FALSE, allow_inf = FALSE) {
  if (!is.numeric(x) || length(x) == 0L) {
    stop(sprintf("'%s' must be numeric.", name), call. = FALSE)
  }
  if (!is.null(len) && length(x) != len) {
    stop(sprintf("'%s' must have length %d.", name, len), call. = FALSE)
  }
  if (anyNA(x)) {
    stop(sprintf("'%s' must not contain missing values.", name),
         call. = FALSE)
  }
  if (!allow_inf && any(is.infinite(x))) {
    stop(sprintf("'%s' must be finite.", name), call. = FALSE)
  }
  if (integer && any(is.finite(x) & abs(x - round(x)) > 1e-8)) {
    stop(sprintf("'%s' must contain whole numbers only.", name),
         call. = FALSE)
  }
  bad_lower <- if (lower_open) any(x <= lower) else any(x < lower)
  bad_upper <- if (upper_open) any(x >= upper) else any(x > upper)
  if (bad_lower || bad_upper) {
    parts <- character(0)
    if (is.finite(lower)) {
      parts <- c(parts, paste(if (lower_open) "greater than" else "at least",
                              format(lower)))
    }
    if (is.finite(upper)) {
      parts <- c(parts, paste(if (upper_open) "less than" else "at most",
                              format(upper)))
    }
    stop(sprintf("'%s' must be %s.", name, paste(parts, collapse = " and ")),
         call. = FALSE)
  }
  invisible(x)
}
