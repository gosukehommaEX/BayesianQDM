test_that("report_scenarios() gives the scenarios and the probabilities", {
  x <- data.frame(pi_t = c(0.2, 0.4), pi_c = 0.4, Go = c(0.01234, 0.5),
                  NoGo = c(0.9, NA), Gray = c(0.08766, 0.2), Miss = 0)
  expect_identical(
    report_scenarios(x, 3L),
    list(type = "table",
         header = c("$\\pi_t$", "$\\pi_c$", "Go", "NoGo", "Gray", "Miss"),
         rows = matrix(c("0.2", "0.4", "0.4", "0.4", "0.012", "0.500",
                         "0.900", "NA", "0.088", "0.200", "0.000", "0.000"),
                       nrow = 2L))
  )
  expect_identical(report_scenarios(x, 1L)$rows[, 3L], c("0.0", "0.5"))
})

test_that("report_scenarios() keeps a column without a symbol and one row", {
  x <- data.frame(Scenario = "A", pi_t = 0.3, Go = 0.1, NoGo = 0.2,
                  Gray = 0.7, Miss = 0)
  res <- report_scenarios(x, 2L)
  expect_identical(res$header[1:2], c("Scenario", "$\\pi_t$"))
  expect_identical(res$rows,
                   matrix(c("A", "0.3", "0.10", "0.20", "0.70", "0.00"),
                          nrow = 1L))
})

test_that("report_scenarios() gives the probabilities of qdm_oc()", {
  oc <- qdm_oc(ref_design("D1"), scen_c)
  res <- report_scenarios(oc, 3L)
  expect_identical(res$rows[, 3L], formatC(oc$Go, format = "f", digits = 3))
  expect_identical(res$rows[, 6L], formatC(oc$Miss, format = "f", digits = 3))
})
