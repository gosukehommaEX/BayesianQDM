#' Names of the scenarios that set a column of both groups
#'
#' @param endpoint The endpoint of a design.
#' @param cols The columns of the scenarios of the design (from
#'   \code{scenario_columns()}).
#'
#' @return A named list: for each name without the letter of the group, the
#'   columns of the scenarios that it sets (those of the control group only
#'   with a concurrent control group).
#' @noRd
scenario_shared <- function(endpoint, cols) {
  shared <- switch(endpoint,
                   continuous = list(sd = c("sd_t", "sd_c")),
                   binary2 = list(rho = c("rho_t", "rho_c")),
                   continuous2 = list(sd1 = c("sd_t1", "sd_c1"),
                                      sd2 = c("sd_t2", "sd_c2"),
                                      rho = c("rho_t", "rho_c")),
                   list())
  lapply(shared, intersect, cols)
}
