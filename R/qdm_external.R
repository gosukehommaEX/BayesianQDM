#' Set the external data of a design specification
#'
#' Sets the arguments \code{external_t} and \code{external_c} of a design
#' specification started by \code{\link{qdm_endpoint}}: external data of the
#' treatment and control groups, borrowed through a power prior. With
#' \code{n_c = 0} (see \code{\link{qdm_sample_size}}) and external control
#' data, the control group is the external control data only. The arguments
#' are checked when the design is built, by the \code{design_*()} function
#' of the endpoint.
#'
#' @param x A design specification started by \code{\link{qdm_endpoint}}.
#' @param external_t,external_c External data of the treatment and control
#'   groups, created by \code{\link{borrow}} with the \code{data_*()}
#'   function of the endpoint, or \code{NULL} (no external data).
#'
#' @return \code{x} with the external data set.
#'
#' @family design steps
#' @examples
#' qdm_endpoint("binary") |>
#'   qdm_external(external_c = borrow(data_binary(n = 20, y = 8), a0 = 0.5))
#' @export
qdm_external <- function(x, external_t = NULL, external_c = NULL) {
  spec_set(x, list(external_t = external_t, external_c = external_c),
           "qdm_external")
}
