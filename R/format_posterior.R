#' Format the parameters of a posterior distribution
#'
#' @param design An object of class \code{"qdm_design"}.
#' @param post Named numeric vector of posterior parameters.
#' @param digits Number of significant digits.
#'
#' @return A character string such as \code{"Beta(a = 7.5, b = 5.5)"}.
#' @noRd
format_posterior <- function(design, post, digits = 4L) {
  switch(design$endpoint,
         binary = sprintf("Beta(a = %s, b = %s)",
                          format_num(post[["a"]], digits),
                          format_num(post[["b"]], digits)))
}
