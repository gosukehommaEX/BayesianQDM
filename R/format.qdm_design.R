#' Format a design
#'
#' @param x An object of class \code{"qdm_design"}.
#' @param digits Number of significant digits for numeric values.
#' @param ... Not used.
#'
#' @return A character vector with one element per line of the description.
#'
#' @examples
#' des <- design_binary(
#'   n_t = 12, n_c = 12,
#'   prior_t = prior_beta(0.5, 0.5), prior_c = prior_beta(0.5, 0.5),
#'   criterion = criterion_posterior(theta_tv = 0.2, theta_mav = 0.1,
#'                                   gamma_go = 0.5, gamma_nogo = 0.8)
#' )
#' format(des)
#' @export
format.qdm_design <- function(x, digits = 4L, ...) {
  endpoint <- switch(x$endpoint,
                     binary = "one binary endpoint",
                     continuous = "one continuous endpoint",
                     binary2 = "two binary endpoints",
                     stop("Unknown design.", call. = FALSE))
  control <- switch(x$control,
                    concurrent = sprintf("concurrent, n = %d", x$n_c),
                    external = "external data only",
                    hypothetical = "hypothetical")
  out <- c(sprintf("Design: %s", endpoint),
           sprintf("  Treatment: n = %d, prior %s", x$n_t,
                   format(x$prior_t, digits = digits)))
  if (!is.null(x$external_t)) {
    out <- c(out, sprintf("    External data: %s",
                          format(x$external_t, digits = digits)))
  }
  out <- c(out, if (is.null(x$prior_c)) {
    sprintf("  Control: %s", control)
  } else {
    sprintf("  Control: %s, prior %s", control,
            format(x$prior_c, digits = digits))
  })
  if (!is.null(x$external_c)) {
    out <- c(out, sprintf("    External data: %s",
                          format(x$external_c, digits = digits)))
  }
  if (!is.null(x$hypothetical_c)) {
    out <- c(out, sprintf("    Hypothetical control: %s",
                          format(x$hypothetical_c, digits = digits)))
  }
  out <- c(out, format(x$criterion, digits = digits))
  if (!is.null(x$method)) {
    out <- c(out, sprintf("Probabilities: %s", method_label(x$method,
                                                             x$nMC)))
  }
  if (x$endpoint == "binary2") {
    out <- c(out,
             sprintf("Probabilities: %s",
                     if (x$criterion$type == "posterior") {
                       method_label("MC", x$nMC)
                     } else {
                       "exact computation"
                     }),
             sprintf("Operating characteristics: %s",
                     if (x$enumerate) "enumeration of all outcomes" else
                       "simulation"))
  }
  out
}
