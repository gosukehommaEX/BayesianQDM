#' Posterior distribution of a group in a report
#'
#' @param design An object of class \code{"qdm_design"}.
#' @param post Named numeric vector of the parameters of the posterior
#'   distribution (\code{posterior_t} or \code{posterior_c} of a result of
#'   \code{qdm_decide()}), not a fixed control value.
#' @param k Group: \code{"t"} (treatment) or \code{"c"} (control).
#'
#' @return A character vector with the parameter (in LaTeX) and its posterior
#'   distribution, to be used in a sentence, for example \code{"$\\pi_t$"}
#'   and \code{"$\\textrm{Beta}(10, 4)$"}. Each number is formatted
#'   separately (see \code{report_num()}).
#' @noRd
report_posterior <- function(design, post, k) {
  num <- report_num
  dof <- function(v) {
    sprintf("%s degree%s of freedom", num(v), if (v == 1) "" else "s")
  }
  switch(
    design$endpoint,
    binary = c(sprintf("$\\pi_%s$", k),
               sprintf("$\\textrm{Beta}(%s, %s)$", num(post[["a"]]),
                       num(post[["b"]]))),
    continuous = c(sprintf("$\\mu_%s$", k), sprintf(
      "the t distribution with %s, location %s and scale %s",
      dof(post[["nu"]]), num(post[["mu"]]),
      num(post[["sigma"]] / sqrt(post[["kappa"]]))
    )),
    binary2 = c(sprintf("$p_%s$", k),
                sprintf("$\\textrm{Dirichlet}(%s)$", num(post))),
    continuous2 = {
      f <- 1 / (post[["kappa"]] * (post[["nu"]] - 1))
      c(sprintf("$\\mu_%s$", k), sprintf(
        paste("the bivariate t distribution with %s, location (%s) and scale",
              "matrix $\\begin{pmatrix} %s & %s \\\\ %s & %s",
              "\\end{pmatrix}$"),
        dof(post[["nu"]] - 1), num(post[c("mu1", "mu2")]),
        num(f * post[["L11"]]), num(f * post[["L12"]]),
        num(f * post[["L12"]]), num(f * post[["L22"]])
      ))
    },
    tte = c(sprintf("$\\lambda_%s$", k),
            sprintf("$\\textrm{Gamma}(%s, %s)$ (shape and rate)",
                    num(post[["shape"]]), num(post[["rate"]])))
  )
}
