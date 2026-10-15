#' Format the parameters of a posterior distribution
#'
#' @param design An object of class \code{"qdm_design"}.
#' @param post Named numeric vector of posterior parameters, or the fixed
#'   control value(s) named \code{value} (one endpoint) or \code{value1} and
#'   \code{value2} (two endpoints).
#' @param digits Number of significant digits.
#'
#' @return A character string such as \code{"Beta(a = 7.5, b = 5.5)"} for a
#'   binary endpoint, or the t distribution of the mean, such as
#'   \code{"mean ~ t(df = 11, location = 18, scale = 7.217)"}, for a
#'   continuous endpoint, or \code{"Dirichlet(alpha = (2.25, 1.25, 3.25,
#'   4.25))"} for two binary endpoints, or
#'   \code{"hazard ~ Gamma(shape = 5.1, rate = 41)"} for a time-to-event
#'   endpoint, or \code{"fixed value 0.2"} for a fixed control
#'   (\code{"fixed median survival time 6"} for the hazard ratio).
#' @noRd
format_posterior <- function(design, post, digits = 4L) {
  if (startsWith(names(post)[1L], "value")) {
    if (identical(design$effect, "hazard_ratio")) {
      return(sprintf("fixed median survival time %s",
                     format_num(post, digits)))
    }
    return(if (length(post) == 1L) {
      sprintf("fixed value %s", format_num(post, digits))
    } else {
      sprintf("fixed values (%s)", format_num(post, digits))
    })
  }
  switch(design$endpoint,
         binary = sprintf("Beta(a = %s, b = %s)",
                          format_num(post[["a"]], digits),
                          format_num(post[["b"]], digits)),
         continuous = sprintf(
           "mean ~ t(df = %s, location = %s, scale = %s)",
           format_num(post[["nu"]], digits), format_num(post[["mu"]], digits),
           format_num(post[["sigma"]] / sqrt(post[["kappa"]]), digits)
         ),
         binary2 = sprintf("Dirichlet(alpha = (%s))",
                           format_num(post, digits)),
         continuous2 = {
           f <- 1 / (post[["kappa"]] * (post[["nu"]] - 1))
           sprintf("mean ~ t(df = %s, location = (%s), scale = (%s; %s))",
                   format_num(post[["nu"]] - 1, digits),
                   format_num(post[c("mu1", "mu2")], digits),
                   format_num(f * post[c("L11", "L12")], digits),
                   format_num(f * post[c("L12", "L22")], digits))
         },
         tte = sprintf("hazard ~ Gamma(shape = %s, rate = %s)",
                       format_num(post[["shape"]], digits),
                       format_num(post[["rate"]], digits)))
}
