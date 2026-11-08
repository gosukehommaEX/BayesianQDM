# The parameters are those of test-format_posterior.R
test_that("report_posterior() gives the Beta and Dirichlet distributions", {
  expect_identical(report_posterior(ref_design("D1"), c(a = 7.5, b = 5.5), "t"),
                   c("$\\pi_t$", "$\\textrm{Beta}(7.5, 5.5)$"))
  expect_identical(report_posterior(ref_design("D1"), c(a = 5, b = 9), "c"),
                   c("$\\pi_c$", "$\\textrm{Beta}(5, 9)$"))
  post <- c("00" = 2.25, "01" = 1.25, "10" = 3.25, "11" = 4.25)
  expect_identical(report_posterior(ref_design_bin2("B1"), post, "t"),
                   c("$p_t$", "$\\textrm{Dirichlet}(2.25, 1.25, 3.25, 4.25)$"))
})

test_that("report_posterior() gives the t distribution of a continuous mean", {
  des <- ref_design_cont("C1")
  expect_identical(
    report_posterior(des, c(mu = 18, kappa = 12, nu = 11, sigma = 25), "t"),
    c("$\\mu_t$", paste("the t distribution with 11 degrees of freedom,",
                        "location 18 and scale 7.217"))
  )
  # One degree of freedom; the scale is 2 / sqrt(4)
  expect_identical(
    report_posterior(des, c(mu = -1, kappa = 4, nu = 1, sigma = 2), "c"),
    c("$\\mu_c$", paste("the t distribution with 1 degree of freedom,",
                        "location -1 and scale 1"))
  )
})

test_that("report_posterior() gives the t distribution of two means", {
  # kappa = 10, nu = 11: the scale matrix is Lambda / 100
  post <- c(mu1 = 1, mu2 = -2, kappa = 10, nu = 11, L11 = 200, L12 = 50,
            L22 = 100)
  expect_identical(
    report_posterior(ref_design_cont2("K1"), post, "t"),
    c("$\\mu_t$",
      paste("the bivariate t distribution with 10 degrees of freedom,",
            "location (1, -2) and scale matrix $\\begin{pmatrix} 2 & 0.5",
            "\\\\ 0.5 & 1 \\end{pmatrix}$"))
  )
})

test_that("report_posterior() gives the Gamma distribution of a hazard", {
  expect_identical(
    report_posterior(ref_design_tte("T2"), c(shape = 5.1, rate = 41), "c"),
    c("$\\lambda_c$", "$\\textrm{Gamma}(5.1, 41)$ (shape and rate)")
  )
})
