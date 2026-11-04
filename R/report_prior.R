#' Prior distribution in a report
#'
#' @param prior A prior distribution created by a \code{prior_*()} function.
#' @param k Group: \code{"t"} (treatment) or \code{"c"} (control).
#' @param endpoint Type of endpoint of the design, used for the vague prior.
#'
#' @return The prior distribution with LaTeX formulas, to be used in a
#'   sentence.
#' @noRd
report_prior <- function(prior, k, endpoint) {
  num <- report_num
  switch(
    class(prior)[1L],
    qdm_prior_beta = sprintf("$\\pi_%s \\sim \\textrm{Beta}(%s, %s)$", k,
                             num(prior$a), num(prior$b)),
    qdm_prior_dirichlet = sprintf("$p_%s \\sim \\textrm{Dirichlet}(%s)$", k,
                                  num(prior$alpha)),
    qdm_prior_nix = sprintf(
      paste("$\\mu_%s \\mid \\sigma_%s^2 \\sim N(%s, \\sigma_%s^2 / %s)$ and",
            "$\\sigma_%s^2 \\sim \\textrm{Inv-}\\chi^2(%s, %s^2)$"),
      k, k, num(prior$mu0), k, num(prior$kappa0), k, num(prior$nu0),
      num(prior$sigma0)
    ),
    qdm_prior_niw = sprintf(
      paste("$\\mu_%s \\mid \\Sigma_%s \\sim N_2((%s), \\Sigma_%s / %s)$ and",
            "$\\Sigma_%s \\sim \\textrm{Inv-Wishart}(%s, \\Lambda_{0%s})$",
            "with $\\Lambda_{0%s} = \\begin{pmatrix} %s & %s \\\\ %s & %s",
            "\\end{pmatrix}$"),
      k, k, num(prior$mu0), k, num(prior$kappa0), k, num(prior$nu0), k, k,
      num(prior$Lambda0[1L, 1L]), num(prior$Lambda0[1L, 2L]),
      num(prior$Lambda0[2L, 1L]), num(prior$Lambda0[2L, 2L])
    ),
    qdm_prior_vague = if (endpoint == "continuous2") {
      sprintf(paste("the vague prior $p(\\mu_%s, \\Sigma_%s) \\propto",
                    "|\\Sigma_%s|^{-3/2}$"), k, k, k)
    } else {
      sprintf(paste("the vague prior $p(\\mu_%s, \\sigma_%s^2) \\propto",
                    "\\sigma_%s^{-2}$"), k, k, k)
    },
    qdm_prior_gamma = sprintf(
      "$\\lambda_%s \\sim \\textrm{Gamma}(%s, %s)$ (shape and rate)", k,
      num(prior$shape), num(prior$rate)
    )
  )
}
