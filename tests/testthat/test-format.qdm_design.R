test_that("format.qdm_design() describes a design with concurrent control", {
  expect_identical(
    format(ref_design("D2")),
    c("Design: one binary endpoint",
      "  Treatment: n = 10, prior Beta(a = 1, b = 1)",
      "    External data: binary, n = 10, y = 6, a0 = 0.5",
      "  Control: concurrent, n = 8, prior Beta(a = 2, b = 3)",
      "    External data: binary, n = 20, y = 8, a0 = 0.5",
      "Posterior probability criterion",
      "  Go:   P(theta > 0.2 | data) >= 0.5",
      "  NoGo: P(theta <= 0.1 | data) >= 0.8")
  )
})

test_that("format.qdm_design() describes the other control forms", {
  f3 <- format(ref_design("D3"))
  expect_identical(f3[3:4],
                   c("  Control: hypothetical, prior Beta(a = 1, b = 1)",
                     "    Hypothetical control: data (binary, n = 12, y = 4)"))
  f4 <- format(ref_design("D4"))
  expect_identical(f4[3:4],
                   c("  Control: external data only, prior Beta(a = 1, b = 1)",
                     "    External data: binary, n = 20, y = 6, a0 = 0.5"))
  f5 <- format(ref_design("D5"))
  expect_identical(f5[4], paste("Posterior predictive probability criterion",
                                "(m_t = 30, m_c = 20)"))
})

test_that("format.qdm_design() rejects unknown endpoints", {
  des <- ref_design("D1")
  des$endpoint <- "unknown"
  expect_error(format(des), "Unknown design.")
})
