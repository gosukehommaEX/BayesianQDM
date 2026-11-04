test_that("report_prior() writes each prior distribution with LaTeX formulas", {
  expect_identical(report_prior(prior_beta(0.5, 0.5), "t", "binary"),
                   "$\\pi_t \\sim \\textrm{Beta}(0.5, 0.5)$")
  expect_identical(report_prior(prior_dirichlet(rep(0.25, 4)), "c", "binary2"),
                   "$p_c \\sim \\textrm{Dirichlet}(0.25, 0.25, 0.25, 0.25)$")
  expect_identical(
    report_prior(prior_nix(5, 2, 4, 20), "t", "continuous"),
    paste("$\\mu_t \\mid \\sigma_t^2 \\sim N(5, \\sigma_t^2 / 2)$ and",
          "$\\sigma_t^2 \\sim \\textrm{Inv-}\\chi^2(4, 20^2)$")
  )
  expect_identical(
    report_prior(prior_niw(c(0, 1), 0.1, 4, matrix(c(25, 1, 1, 9), 2L)), "c",
                 "continuous2"),
    paste("$\\mu_c \\mid \\Sigma_c \\sim N_2((0, 1), \\Sigma_c / 0.1)$ and",
          "$\\Sigma_c \\sim \\textrm{Inv-Wishart}(4, \\Lambda_{0c})$ with",
          "$\\Lambda_{0c} = \\begin{pmatrix} 25 & 1 \\\\ 1 & 9",
          "\\end{pmatrix}$")
  )
  expect_identical(report_prior(prior_gamma(0.1, 1), "t", "tte"),
                   "$\\lambda_t \\sim \\textrm{Gamma}(0.1, 1)$ (shape and rate)")
})

test_that("report_prior() writes the vague prior of one or two endpoints", {
  expect_identical(
    report_prior(prior_vague(), "t", "continuous"),
    "the vague prior $p(\\mu_t, \\sigma_t^2) \\propto \\sigma_t^{-2}$"
  )
  expect_identical(
    report_prior(prior_vague(), "c", "continuous2"),
    "the vague prior $p(\\mu_c, \\Sigma_c) \\propto |\\Sigma_c|^{-3/2}$"
  )
})
