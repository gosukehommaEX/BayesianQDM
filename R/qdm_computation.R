#' Set the computation of the probabilities of a design specification
#'
#' Sets the arguments of the \code{design_*()} function that choose how the
#' probabilities of the regions are computed. Each endpoint has only some of
#' them: \code{method} for one or two continuous endpoints, \code{nMC} for
#' one continuous endpoint, two binary endpoints, two continuous endpoints and
#' a time-to-event endpoint, and \code{enumerate} for two binary endpoints. An
#' argument that the \code{design_*()} function of the endpoint does not have
#' is an error. The values are checked when the design is built, by the
#' \code{design_*()} function. All the arguments of the step are set together:
#' an argument that is not given returns to the default, so
#' \code{qdm_computation(x, nMC = 5000)} also returns \code{method} to its
#' default.
#'
#' @param x A design specification started by \code{\link{qdm_endpoint}},
#'   or a design, which is converted to the specification that builds it.
#' @param method Method of computation of the probabilities:
#'   \code{"NI"}, \code{"MM"} or \code{"MC"} for one continuous endpoint
#'   (see \code{\link{design_continuous}}), \code{"MM"} or \code{"MC"} for
#'   two continuous endpoints (see \code{\link{design_continuous2}}), or
#'   \code{NULL} for the default.
#' @param nMC Number of Monte Carlo draws or simulated future studies, or
#'   \code{NULL} for the default.
#' @param enumerate For two binary endpoints, whether the operating
#'   characteristics enumerate all outcomes (see
#'   \code{\link{design_binary2}}), or \code{NULL} for the default.
#'
#' @return A design specification: \code{x}, or the specification of the
#'   design \code{x}, with the arguments set.
#'
#' @family design steps
#' @examples
#' qdm_endpoint("continuous") |>
#'   qdm_computation(method = "MM")
#' @export
qdm_computation <- function(x, method = NULL, nMC = NULL, enumerate = NULL) {
  spec_set(x, list(method = method, nMC = nMC, enumerate = enumerate),
           "qdm_computation")
}
