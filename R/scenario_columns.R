#' Columns of the scenarios of a design
#'
#' @param design An object of class \code{"qdm_design"}.
#'
#' @return A list with elements \code{names} (the required columns),
#'   \code{lower} and \code{upper} (bounds of their values) and
#'   \code{lower_open} (whether the lower bound is excluded), each of the
#'   length of \code{names}.
#' @noRd
scenario_columns <- function(design) {
  concurrent <- design$control == "concurrent"
  cols <- switch(
    design$endpoint,
    binary = list(names = c("pi_t", "pi_c"), lower = c(0, 0),
                  upper = c(1, 1), lower_open = c(FALSE, FALSE)),
    continuous = list(names = c("mu_t", "sd_t", "mu_c", "sd_c"),
                      lower = c(-Inf, 0, -Inf, 0),
                      upper = c(Inf, Inf, Inf, Inf),
                      lower_open = c(FALSE, TRUE, FALSE, TRUE))
  )
  # Columns of the treatment group come first; those of the control group
  # are required only with a concurrent control group
  keep <- if (concurrent) {
    seq_along(cols$names)
  } else {
    seq_len(length(cols$names) / 2L)
  }
  lapply(cols, `[`, keep)
}
