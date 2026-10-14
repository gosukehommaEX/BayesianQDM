#' Check the scenarios of operating characteristics
#'
#' Checks that the scenarios are a data frame with the columns required by
#' the design, with values in the allowed range.
#'
#' @param design An object of class \code{"qdm_design"}.
#' @param scenarios Object to check.
#' @param name Name of the argument, used in error messages.
#' @param one_row Logical; if \code{TRUE}, exactly one row is required.
#' @param feasible Logical; if \code{TRUE}, scenarios of two binary endpoints
#'   whose correlation is infeasible for their response rates are rejected
#'   (otherwise they are reported as not evaluable by \code{qdm_oc()}).
#'
#' @return \code{scenarios} as a plain data frame.
#' @noRd
check_scenarios <- function(design, scenarios, name = "scenarios",
                            one_row = FALSE, feasible = FALSE) {
  if (!is.data.frame(scenarios) || nrow(scenarios) == 0L) {
    stop(sprintf("'%s' must be a data frame with at least one row.", name),
         call. = FALSE)
  }
  if (one_row && nrow(scenarios) != 1L) {
    stop(sprintf("'%s' must have exactly one row.", name), call. = FALSE)
  }
  scenarios <- as.data.frame(scenarios)
  cols <- scenario_columns(design)
  if (!all(cols$names %in% names(scenarios))) {
    stop(sprintf("'%s' must have the column(s) %s.", name,
                 paste(cols$names, collapse = ", ")), call. = FALSE)
  }
  if (any(c("Go", "NoGo", "Gray", "Miss") %in% names(scenarios))) {
    stop(sprintf("'%s' must not have columns named Go, NoGo, Gray or Miss.",
                 name), call. = FALSE)
  }
  for (k in seq_along(cols$names)) {
    col <- cols$names[k]
    check_numeric(scenarios[[col]], sprintf("%s$%s", name, col), len = NULL,
                  lower = cols$lower[k], upper = cols$upper[k],
                  lower_open = cols$lower_open[k],
                  upper_open = cols$upper_open[k])
  }
  if (feasible && design$endpoint == "binary2") {
    groups <- if (design$control == "concurrent") c("t", "c") else "t"
    for (g in groups) {
      cells <- cellprob_scenarios(scenarios[[paste0("pi_", g, "1")]],
                                  scenarios[[paste0("pi_", g, "2")]],
                                  scenarios[[paste0("rho_", g)]])
      if (anyNA(cells)) {
        stop(sprintf(paste("'%s$rho_%s' is outside the feasible range for",
                           "the response rates (see cellprob_binary2())."),
                     name, g), call. = FALSE)
      }
    }
  }
  scenarios
}
