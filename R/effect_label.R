#' Label of the treatment effect of one endpoint
#'
#' @param design A design for one endpoint.
#'
#' @return A plotmath expression for the axis of the effect.
#' @noRd
effect_label <- function(design) {
  if (identical(design$effect, "hazard_ratio")) {
    return(expression(theta == lambda[t] / lambda[c]))
  }
  switch(design$endpoint,
         binary = expression(theta == pi[t] - pi[c]),
         continuous = expression(theta == mu[t] - mu[c]),
         tte = expression(theta == S[t](t[0]) - S[c](t[0])))
}
