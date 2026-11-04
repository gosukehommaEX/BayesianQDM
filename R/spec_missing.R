#' Required arguments not yet set in a design specification
#'
#' @param x A design specification (class \code{"qdm_spec"}).
#'
#' @return A character vector describing the required arguments that are not
#'   set and the functions that set them, empty if none is missing.
#' @noRd
spec_missing <- function(x) {
  required <- c(n_t = "the sample sizes (qdm_sample_size())",
                n_c = "the sample sizes (qdm_sample_size())",
                prior_t = "the prior distributions (qdm_prior())",
                criterion = "the decision criterion (qdm_criterion())")
  unique(unname(required[setdiff(names(required), names(x$args))]))
}
