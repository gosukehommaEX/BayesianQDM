#' Study design section of a report
#'
#' @param design A design.
#'
#' @return A list of blocks.
#' @noRd
report_study <- function(design) {
  endpoints <- switch(design$endpoint,
                      binary = "one binary endpoint",
                      continuous = "one continuous endpoint",
                      binary2 = "two binary endpoints",
                      continuous2 = "two continuous endpoints",
                      tte = "one time-to-event endpoint")
  groups <- if (design$control == "concurrent") {
    sprintf(paste("%d subjects in the treatment group and %d subjects in",
                  "the control group"), design$n_t, design$n_c)
  } else {
    sprintf(paste("%d subjects in the treatment group, without a concurrent",
                  "control group"), design$n_t)
  }
  list(report_block("heading", text = "Study design"),
       report_block("par", text = sprintf(
         "The proof-of-concept study has %s and will enroll %s.", endpoints,
         groups
       )))
}
