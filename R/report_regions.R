#' Table of the regions of two endpoints in a report
#'
#' The regions are numbered as by \code{qdm_regions()}: the columns are the
#' conditions on endpoint 1 and the rows those on endpoint 2, from the most
#' favorable to the least favorable.
#'
#' @param criterion A decision criterion for two endpoints.
#'
#' @return A block of type \code{"table"}.
#' @noRd
report_regions <- function(criterion) {
  pred <- criterion$type == "predictive"
  th <- if (pred) "\\tilde\\theta" else "\\theta"
  labels <- function(k) {
    s <- paste0(th, "_", k)
    less <- criterion$direction[k] == "less"
    if (pred) {
      v <- report_num(criterion$theta_null[k])
      return(if (less) {
        sprintf(c("$%s < %s$", "$%s \\ge %s$"), s, v)
      } else {
        sprintf(c("$%s > %s$", "$%s \\le %s$"), s, v)
      })
    }
    tv <- report_num(criterion$theta_tv[k])
    mav <- report_num(criterion$theta_mav[k])
    if (less) {
      c(sprintf("$%s < %s$", s, tv), sprintf("$%s \\le %s < %s$", tv, s, mav),
        sprintf("$%s \\ge %s$", s, mav))
    } else {
      c(sprintf("$%s > %s$", s, tv), sprintf("$%s < %s \\le %s$", mav, s, tv),
        sprintf("$%s \\le %s$", s, mav))
    }
  }
  n <- if (pred) 2L else 3L
  cells <- matrix(sprintf("$R_{%d}$", seq_len(n * n)), nrow = n)
  report_block("table", header = c("", labels(1L)),
               rows = unname(cbind(labels(2L), cells)))
}
