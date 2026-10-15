#' Check the prior distribution of the control group of a design
#'
#' The prior distribution of the control group is not used for a fixed
#' control and, for continuous endpoints, for a hypothetical control; it must
#' then be \code{NULL}, so that a specification is not silently ignored.
#' Otherwise it is checked by \code{check}.
#'
#' @param prior_c The argument \code{prior_c} of the design.
#' @param control Form of the control group (see \code{control_form()}).
#' @param check Function that checks the prior, called as
#'   \code{check(prior_c, "prior_c")}.
#' @param unused Forms of the control group for which the prior is not used.
#'
#' @return \code{prior_c}, invisibly.
#' @noRd
check_prior_control <- function(prior_c, control, check, unused = "fixed") {
  if (control %in% unused) {
    if (!is.null(prior_c)) {
      reason <- if (control == "fixed") {
        "whose value is given by 'fixed_c'"
      } else {
        paste("whose distribution is determined by 'hypothetical_c' and the",
              "posterior distribution of the treatment group")
      }
      stop(sprintf("'prior_c' must be NULL for a %s control, %s.", control,
                   reason), call. = FALSE)
    }
  } else {
    check(prior_c, "prior_c")
  }
  invisible(prior_c)
}
