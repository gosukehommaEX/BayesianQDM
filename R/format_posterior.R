#' Format the parameters of a posterior distribution
#'
#' @param design An object of class \code{"qdm_design"}.
#' @param post Named numeric vector of posterior parameters.
#' @param digits Number of significant digits.
#'
#' @return A character string such as \code{"Beta(a = 7.5, b = 5.5)"} for a
#'   binary endpoint, or the t distribution of the mean, such as
#'   \code{"mean ~ t(df = 11, location = 18, scale = 7.217)"}, for a
#'   continuous endpoint, or \code{"Dirichlet(alpha = (2.25, 1.25, 3.25,
#'   4.25))"} for two binary endpoints.
#' @noRd
format_posterior <- function(design, post, digits = 4L) {
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
         })
}
