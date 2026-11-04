#' Format a design specification
#'
#' @param x An object of class \code{"qdm_spec"}, created by
#'   \code{\link{qdm_endpoint}} and the other design steps.
#' @param digits Number of significant digits for numeric values.
#' @param ... Not used.
#'
#' @return A character vector with one element per line of the description.
#'   For a specification that is not complete, it lists the arguments set and
#'   those still needed. For a complete specification, it gives the
#'   description of the design that \code{\link{qdm_build}} creates (see
#'   \code{\link{format.qdm_design}}), or the error that building it gives.
#'
#' @examples
#' spec <- qdm_endpoint("binary") |>
#'   qdm_sample_size(n_t = 12, n_c = 12)
#' format(spec)
#' @export
format.qdm_spec <- function(x, digits = 4L, ...) {
  endpoint <- switch(x$endpoint,
                     binary = "one binary endpoint",
                     continuous = "one continuous endpoint",
                     binary2 = "two binary endpoints",
                     continuous2 = "two continuous endpoints",
                     tte = "one time-to-event endpoint")
  unset <- spec_missing(x)
  if (length(unset) > 0L) {
    fun <- get(design_function(x$endpoint), mode = "function")
    set <- intersect(names(formals(fun)), names(x$args))
    return(c(sprintf("Design specification: %s (incomplete)", endpoint),
             sprintf("  Arguments set: %s",
                     if (length(set) == 0L) "none" else
                       paste(set, collapse = ", ")),
             sprintf("  Still needed: %s", paste(unset, collapse = ", "))))
  }
  des <- tryCatch(qdm_build(x), error = function(e) e)
  if (inherits(des, "error")) {
    return(c(sprintf("Design specification: %s (not valid)", endpoint),
             paste0("  ", conditionMessage(des))))
  }
  c("Design specification (complete):", format(des, digits = digits))
}
