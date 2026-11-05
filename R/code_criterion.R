#' Decision criterion as R code
#'
#' @param x An object of class \code{"qdm_criterion"}.
#'
#' @return A call from \code{code_node()} of \code{criterion_posterior()} or
#'   \code{criterion_predictive()}, without the arguments equal to their
#'   defaults.
#' @noRd
code_criterion <- function(x) {
  two <- x$n_endpoints == 2L
  # A threshold given once is used for both endpoints
  gamma <- function(g) {
    if (length(g) == 2L && identical(g[1L], g[2L])) g[1L] else g
  }
  regions <- function(r) {
    if (length(r) > 1L && all(diff(r) == 1L)) {
      sprintf("%d:%d", r[1L], r[length(r)])
    } else {
      code_num(r, keep_type = FALSE)
    }
  }
  if (x$type == "posterior") {
    fun <- "criterion_posterior"
    args <- list(theta_tv = code_num(x$theta_tv),
                 theta_mav = code_num(x$theta_mav))
    nogo_default <- 9L
  } else {
    fun <- "criterion_predictive"
    args <- list(theta_null = code_num(x$theta_null),
                 m_t = code_num(x$m_t, keep_type = FALSE))
    if (!identical(x$m_c, x$m_t)) {
      args$m_c <- code_num(x$m_c, keep_type = FALSE)
    }
    nogo_default <- 4L
  }
  args$gamma_go <- code_num(gamma(x$gamma_go))
  args$gamma_nogo <- code_num(gamma(x$gamma_nogo))
  if (two && x$rule == "joint") {
    if (!identical(x$go_regions, 1L)) {
      args$go_regions <- regions(x$go_regions)
    }
    if (!identical(x$nogo_regions, nogo_default)) {
      args$nogo_regions <- regions(x$nogo_regions)
    }
  }
  if (two && x$rule != "joint") {
    args$rule <- code_value(x$rule)
  }
  if (any(x$direction != "greater")) {
    direction <- if (length(unique(x$direction)) == 1L) {
      x$direction[1L]
    } else {
      x$direction
    }
    args$direction <- code_value(direction)
  }
  code_node(fun, args)
}
