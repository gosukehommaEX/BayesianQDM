#' Sections of a report that describe a design
#'
#' @param design A design.
#'
#' @return A list of blocks: the study design, the statistical model, the
#'   decision criteria and the computation of the probabilities.
#' @noRd
report_design <- function(design) {
  c(report_study(design), report_model(design), report_criteria(design),
    report_computation(design))
}
