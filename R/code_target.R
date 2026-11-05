#' Design or specification whose R code is written
#'
#' @param x The argument \code{x} of \code{show_code()}.
#'
#' @return A design specification or a design: \code{x} itself, or the design
#'   of a result of an analysis.
#' @noRd
code_target <- function(x) {
  if (inherits(x, c("qdm_spec", "qdm_design"))) {
    return(x)
  }
  if (inherits(x, c("qdm_oc", "qdm_boundary"))) {
    design <- attr(x, "design")
  } else if (inherits(x, c("qdm_calibration", "qdm_n_search",
                           "qdm_decision"))) {
    design <- x$design
  } else {
    stop(paste("'x' must be a design, a design specification, or a result",
               "of qdm_oc(), qdm_calibrate(), qdm_find_n(), qdm_boundary()",
               "or qdm_decide()."), call. = FALSE)
  }
  if (!inherits(design, "qdm_design")) {
    stop(paste("'x' has no design (a calibration or a search of the sample",
               "size that did not find one)."), call. = FALSE)
  }
  design
}
