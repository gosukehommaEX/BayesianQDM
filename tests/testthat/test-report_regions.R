test_that("report_regions() gives the nine regions of the posterior criterion", {
  crit <- criterion_posterior(c(0.15, 0.25), c(0.1, 0.15), gamma_go = 0.5,
                              gamma_nogo = 0.5)
  expect_identical(
    report_regions(crit),
    list(type = "table",
         header = c("", "$\\theta_1 > 0.15$", "$0.1 < \\theta_1 \\le 0.15$",
                    "$\\theta_1 \\le 0.1$"),
         rows = unname(cbind(c("$\\theta_2 > 0.25$",
                               "$0.15 < \\theta_2 \\le 0.25$",
                               "$\\theta_2 \\le 0.15$"),
                             matrix(sprintf("$R_{%d}$", 1:9), nrow = 3L))))
  )
  crit2 <- criterion_posterior(c(10, 1), c(5, 2), gamma_go = 0.5,
                               gamma_nogo = 0.3,
                               direction = c("greater", "less"))
  expect_identical(report_regions(crit2)$rows[, 1L],
                   c("$\\theta_2 < 1$", "$1 \\le \\theta_2 < 2$",
                     "$\\theta_2 \\ge 2$"))
})

test_that("report_regions() gives the four regions of the predictive criterion", {
  crit <- criterion_predictive(c(0.1, 0), m_t = 4, m_c = 3, gamma_go = 0.4,
                               gamma_nogo = 0.3)
  expect_identical(
    report_regions(crit),
    list(type = "table",
         header = c("", "$\\tilde\\theta_1 > 0.1$",
                    "$\\tilde\\theta_1 \\le 0.1$"),
         rows = unname(cbind(c("$\\tilde\\theta_2 > 0$",
                               "$\\tilde\\theta_2 \\le 0$"),
                             matrix(sprintf("$R_{%d}$", 1:4), nrow = 2L))))
  )
  crit2 <- criterion_predictive(c(0.1, 0), m_t = 4, m_c = 3, gamma_go = 0.4,
                                gamma_nogo = 0.3, direction = "less")
  expect_identical(report_regions(crit2)$header[-1L],
                   c("$\\tilde\\theta_1 < 0.1$", "$\\tilde\\theta_1 \\ge 0.1$"))
})

test_that("report_regions() numbers the regions as qdm_regions()", {
  for (prob in c("posterior", "predictive")) {
    crit <- if (prob == "posterior") {
      criterion_posterior(c(0.15, 0.25), c(0.1, 0.15), gamma_go = 0.5,
                          gamma_nogo = 0.5)
    } else {
      criterion_predictive(c(0.1, 0), m_t = 4, m_c = 3, gamma_go = 0.4,
                           gamma_nogo = 0.3)
    }
    rows <- report_regions(crit)$rows
    expect_identical(unname(gsub("[${}_]", "", rows[, -1L, drop = FALSE])),
                     unname(unclass(qdm_regions(prob))), info = prob)
  }
})
