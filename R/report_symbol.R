#' Symbols of the columns of the scenarios in a report
#'
#' @param name Names of columns of the scenarios of \code{qdm_oc()} and
#'   \code{qdm_calibrate()}.
#'
#' @return The LaTeX symbol of each column, or \code{NA} for a column without
#'   a symbol.
#' @noRd
report_symbol <- function(name) {
  sym <- c(pi_t = "\\pi_t", pi_c = "\\pi_c", mu_t = "\\mu_t",
           sd_t = "\\sigma_t", mu_c = "\\mu_c", sd_c = "\\sigma_c",
           pi_t1 = "\\pi_{t1}", pi_t2 = "\\pi_{t2}", rho_t = "\\rho_t",
           pi_c1 = "\\pi_{c1}", pi_c2 = "\\pi_{c2}", rho_c = "\\rho_c",
           mu_t1 = "\\mu_{t1}", mu_t2 = "\\mu_{t2}", sd_t1 = "\\sigma_{t1}",
           sd_t2 = "\\sigma_{t2}", mu_c1 = "\\mu_{c1}", mu_c2 = "\\mu_{c2}",
           sd_c1 = "\\sigma_{c1}", sd_c2 = "\\sigma_{c2}", median_t = "M_t",
           median_c = "M_c")
  unname(sym[name])
}
