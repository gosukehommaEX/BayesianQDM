#' Blocks of a report
#'
#' @param x A design, a design specification, or a result of
#'   \code{qdm_oc()}, \code{qdm_calibrate()}, \code{qdm_find_n()},
#'   \code{qdm_boundary()} or \code{qdm_decide()}.
#' @param digits Number of decimal places of the probabilities.
#'
#' @return A list of blocks (see \code{report_block()}): those of the design,
#'   followed by those of the result.
#' @noRd
report_blocks <- function(x, digits) {
  if (inherits(x, "qdm_spec")) {
    x <- qdm_build(x)
  }
  if (inherits(x, "qdm_design")) {
    return(report_design(x))
  }
  if (inherits(x, "qdm_oc")) {
    return(c(report_design(attr(x, "design")), report_oc(x, digits)))
  }
  if (inherits(x, "qdm_calibration")) {
    if (is.null(x$design)) {
      stop(paste("The calibration did not find both thresholds, so there is",
                 "no design to report."), call. = FALSE)
    }
    return(c(report_design(x$design), report_calibration(x, digits)))
  }
  if (inherits(x, "qdm_n_search")) {
    if (is.null(x$design)) {
      stop(paste("The search did not find a sample size that meets the",
                 "targets, so there is no design to report."), call. = FALSE)
    }
    return(c(report_design(x$design), report_n_search(x, digits)))
  }
  if (inherits(x, "qdm_boundary")) {
    return(c(report_design(attr(x, "design")), report_boundary(x)))
  }
  if (inherits(x, "qdm_decision")) {
    return(c(report_design(x$design), report_decision(x, digits)))
  }
  stop(paste("'x' must be a design, a design specification, or a result of",
             "qdm_oc(), qdm_calibrate(), qdm_find_n(), qdm_boundary() or",
             "qdm_decide()."), call. = FALSE)
}
