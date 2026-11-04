#' Probabilities of the criteria in a report
#'
#' LaTeX expressions of the probabilities that the decision criteria compare
#' with the thresholds, with the values of the criterion.
#'
#' @param criterion A decision criterion.
#'
#' @return A list with the elements \code{go} and \code{nogo}: one expression
#'   each, or one per endpoint for two endpoints with a marginal rule.
#' @noRd
report_events <- function(criterion) {
  pred <- criterion$type == "predictive"
  th <- if (pred) "\\tilde\\theta" else "\\theta"
  if (criterion$n_endpoints == 2L && criterion$rule == "joint") {
    regions <- function(r) paste0("R_{", r, "}", collapse = " \\cup ")
    pair <- sprintf("(%s_1, %s_2)", th, th)
    return(list(
      go = sprintf("\\Pr\\{%s \\in %s \\mid D\\}", pair,
                   regions(criterion$go_regions)),
      nogo = sprintf("\\Pr\\{%s \\in %s \\mid D\\}", pair,
                     regions(criterion$nogo_regions))
    ))
  }
  sym <- if (criterion$n_endpoints == 2L) paste0(th, "_", 1:2) else th
  less <- criterion$direction == "less"
  if (pred) {
    v_go <- v_nogo <- vapply(criterion$theta_null, report_num, character(1))
  } else {
    v_go <- vapply(criterion$theta_tv, report_num, character(1))
    v_nogo <- vapply(criterion$theta_mav, report_num, character(1))
  }
  list(go = sprintf("\\Pr(%s %s %s \\mid D)", sym, ifelse(less, "<", ">"),
                    v_go),
       nogo = sprintf("\\Pr(%s %s %s \\mid D)", sym,
                      ifelse(less, "\\ge", "\\le"), v_nogo))
}
