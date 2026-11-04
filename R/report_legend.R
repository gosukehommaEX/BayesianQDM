#' Meaning of the columns of the scenarios in a report
#'
#' @param design A design.
#'
#' @return A character string.
#' @noRd
report_legend <- function(design) {
  conc <- design$control == "concurrent"
  switch(
    design$endpoint,
    binary = if (conc) {
      paste("$\\pi_t$ and $\\pi_c$ are the true response rates of the",
            "treatment and control groups.")
    } else {
      "$\\pi_t$ is the true response rate of the treatment group."
    },
    continuous = paste("$\\mu_k$ and $\\sigma_k$ are the true mean and",
                       "standard deviation of group $k$."),
    binary2 = paste("$\\pi_{kj}$ is the true response rate of endpoint $j$",
                    "in group $k$, and $\\rho_k$ is the correlation between",
                    "the responses on the two endpoints in group $k$."),
    continuous2 = paste("$\\mu_{kj}$ and $\\sigma_{kj}$ are the true mean",
                        "and standard deviation of endpoint $j$ in group",
                        "$k$, and $\\rho_k$ is the correlation between the",
                        "two endpoints in group $k$."),
    tte = if (conc) {
      paste("$M_t$ and $M_c$ are the true median survival times of the",
            "treatment and control groups.")
    } else {
      "$M_t$ is the true median survival time of the treatment group."
    }
  )
}
