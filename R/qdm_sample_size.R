#' Set the sample sizes of a design specification
#'
#' Sets the arguments \code{n_t} and \code{n_c} of a design specification
#' started by \code{\link{qdm_endpoint}}. They are checked when the design is
#' built, by the \code{design_*()} function of the endpoint.
#'
#' @param x A design specification started by \code{\link{qdm_endpoint}}.
#' @param n_t Number of subjects in the treatment group.
#' @param n_c Number of subjects in the concurrent control group: 0 for a
#'   study without a concurrent control group (external control data only, a
#'   hypothetical control or a fixed control).
#'
#' @return \code{x} with the sample sizes set.
#'
#' @family design steps
#' @examples
#' qdm_endpoint("binary") |>
#'   qdm_sample_size(n_t = 12, n_c = 12)
#' @export
qdm_sample_size <- function(x, n_t, n_c) {
  spec_set(x, list(n_t = n_t, n_c = n_c), "qdm_sample_size")
}
