#' Start a design specification for an endpoint
#'
#' Starts the specification of a design that is built step by step, for
#' example with the native pipe \code{|>} of R. The other steps set the
#' arguments of the \code{design_*()} function of the endpoint:
#' \describe{
#'   \item{\code{\link{qdm_sample_size}}}{\code{n_t} and \code{n_c};}
#'   \item{\code{\link{qdm_prior}}}{\code{prior_t} and \code{prior_c};}
#'   \item{\code{\link{qdm_criterion}}}{\code{criterion};}
#'   \item{\code{\link{qdm_external}}}{\code{external_t} and
#'     \code{external_c};}
#'   \item{\code{\link{qdm_control}}}{\code{hypothetical_c} and
#'     \code{fixed_c};}
#'   \item{\code{\link{qdm_tte}}}{\code{milestone}, \code{effect},
#'     \code{conduct} and \code{conduct_future} (time-to-event endpoint);}
#'   \item{\code{\link{qdm_computation}}}{\code{method}, \code{nMC} and
#'     \code{enumerate}.}
#' }
#' \code{\link{qdm_build}} creates the design with the \code{design_*()}
#' function of the endpoint.
#'
#' The steps only store their arguments, so they can be given in any order.
#' A step given again sets all its arguments anew: an argument that is not
#' given, or is given as \code{NULL}, returns to the default of the
#' \code{design_*()} function. The arguments are checked when the design is
#' built, by the \code{design_*()} function, so the design and its error
#' messages are those of that function. \code{\link{qdm_decide}},
#' \code{\link{qdm_oc}}, \code{\link{qdm_compare}},
#' \code{\link{qdm_calibrate}}, \code{\link{qdm_find_n}},
#' \code{\link{qdm_boundary}} and \code{\link{qdm_scenarios}} accept a
#' specification in place of a design and build it, and
#' \code{\link{show_report}} and \code{\link{show_code}} accept it as well.
#'
#' A design created by a \code{design_*()} function can also be given to the
#' steps and to \code{\link{qdm_build}}: it is converted to the
#' specification that builds it (the arguments that are \code{NULL} in the
#' design are not set), so that one part of a design can be changed with a
#' step, for example
#' \code{design_binary(...) |> qdm_sample_size(n_t = 20, n_c = 20)}.
#' As with a step given again, all the arguments of the step are set: for
#' example, \code{qdm_external(des, external_c = ...)} also removes the
#' external data of the treatment group of \code{des}.
#'
#' @param endpoint Type of endpoint: \code{"binary"} (one binary endpoint,
#'   \code{\link{design_binary}}), \code{"continuous"} (one continuous
#'   endpoint, \code{\link{design_continuous}}), \code{"binary2"} (two binary
#'   endpoints, \code{\link{design_binary2}}), \code{"continuous2"} (two
#'   continuous endpoints, \code{\link{design_continuous2}}) or \code{"tte"}
#'   (one time-to-event endpoint, \code{\link{design_tte}}).
#'
#' @return An object of class \code{"qdm_spec"}, a list with the elements
#'   \code{endpoint} and \code{args} (the arguments of the \code{design_*()}
#'   function set so far, none at the start).
#'
#' @family design steps
#' @examples
#' spec <- qdm_endpoint("binary") |>
#'   qdm_sample_size(n_t = 12, n_c = 12) |>
#'   qdm_prior(prior_t = prior_beta(0.5, 0.5),
#'             prior_c = prior_beta(0.5, 0.5)) |>
#'   qdm_criterion(criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                     gamma_go = 0.5, gamma_nogo = 0.8))
#' spec
#'
#' # The same design as with design_binary()
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' identical(qdm_build(spec), des)
#'
#' # The specification is built by the functions that use a design
#' spec |>
#'   qdm_decide(data_t = data_binary(12, 9), data_c = data_binary(12, 4))
#' @export
qdm_endpoint <- function(endpoint) {
  endpoints <- c("binary", "continuous", "binary2", "continuous2", "tte")
  if (missing(endpoint) || !is.character(endpoint) ||
      length(endpoint) != 1L || !(endpoint %in% endpoints)) {
    stop(sprintf("'endpoint' must be one of %s.",
                 paste0("\"", endpoints, "\"", collapse = ", ")),
         call. = FALSE)
  }
  structure(list(endpoint = endpoint, args = list()), class = "qdm_spec")
}
