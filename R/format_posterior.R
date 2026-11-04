#' Format the parameters of a posterior distribution
#'
#' @param design An object of class \code{"qdm_design"}.
#' @param post Named numeric vector of posterior parameters, or the fixed
#'   control value(s) named \code{value} (one endpoint) or \code{value1} and
#'   \code{value2} (two endpoints).
#' @param digits Number of significant digits.
#' @param num Function that formats a numeric vector as one character
#'   string: by default \code{format_num()} with \code{digits}, which gives
#'   the numbers of a vector the same number of decimals; the reports of
#'   \code{show_report()} use \code{report_num()}.
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
format_posterior <- function(design, post, digits = 4L,
                             num = function(v) format_num(v, digits)) {
  if (startsWith(names(post)[1L], "value")) {
    if (identical(design$effect, "hazard_ratio")) {
      return(sprintf("fixed median survival time %s",
                     num(post)))
    }
    return(if (length(post) == 1L) {
      sprintf("fixed value %s", num(post))
    } else {
      sprintf("fixed values (%s)", num(post))
    })
  }
  switch(design$endpoint,
         binary = sprintf("Beta(a = %s, b = %s)",
                          num(post[["a"]]),
                          num(post[["b"]])),
         continuous = sprintf(
           "mean ~ t(df = %s, location = %s, scale = %s)",
           num(post[["nu"]]), num(post[["mu"]]),
           num(post[["sigma"]] / sqrt(post[["kappa"]]))
         ),
         binary2 = sprintf("Dirichlet(alpha = (%s))",
                           num(post)),
         continuous2 = {
           f <- 1 / (post[["kappa"]] * (post[["nu"]] - 1))
           sprintf("mean ~ t(df = %s, location = (%s), scale = (%s; %s))",
                   num(post[["nu"]] - 1),
                   num(post[c("mu1", "mu2")]),
                   num(f * post[c("L11", "L12")]),
                   num(f * post[c("L12", "L22")]))
         },
         tte = sprintf("hazard ~ Gamma(shape = %s, rate = %s)",
                       num(post[["shape"]]),
                       num(post[["rate"]])))
}
