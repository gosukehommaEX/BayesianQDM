#' Columns of the scenarios of a design
#'
#' @param design An object of class \code{"qdm_design"}.
#'
#' @return A list with elements \code{names} (the required columns),
#'   \code{lower} and \code{upper} (bounds of their values),
#'   \code{lower_open} and \code{upper_open} (whether the bounds are
#'   excluded), each of the length of \code{names}.
#' @noRd
scenario_columns <- function(design) {
  concurrent <- design$control == "concurrent"
  cols <- switch(
    design$endpoint,
    binary = list(names = c("pi_t", "pi_c"), lower = c(0, 0),
                  upper = c(1, 1), lower_open = c(FALSE, FALSE),
                  upper_open = c(FALSE, FALSE)),
    continuous = list(names = c("mu_t", "sd_t", "mu_c", "sd_c"),
                      lower = c(-Inf, 0, -Inf, 0),
                      upper = c(Inf, Inf, Inf, Inf),
                      lower_open = c(FALSE, TRUE, FALSE, TRUE),
                      upper_open = rep(FALSE, 4L)),
    binary2 = list(names = c("pi_t1", "pi_t2", "rho_t", "pi_c1", "pi_c2",
                             "rho_c"),
                   lower = rep(c(0, 0, -1), 2L), upper = rep(1, 6L),
                   lower_open = rep(c(TRUE, TRUE, FALSE), 2L),
                   upper_open = rep(c(TRUE, TRUE, FALSE), 2L)),
    continuous2 = list(names = c("mu_t1", "mu_t2", "sd_t1", "sd_t2", "rho_t",
                                 "mu_c1", "mu_c2", "sd_c1", "sd_c2",
                                 "rho_c"),
                       lower = rep(c(-Inf, -Inf, 0, 0, -1), 2L),
                       upper = rep(c(Inf, Inf, Inf, Inf, 1), 2L),
                       lower_open = rep(c(FALSE, FALSE, TRUE, TRUE, TRUE),
                                        2L),
                       upper_open = rep(c(FALSE, FALSE, FALSE, FALSE, TRUE),
                                        2L))
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
