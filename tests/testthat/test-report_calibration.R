test_that("report_calibration() describes how the thresholds were chosen", {
  cal <- run_cal("D1")
  prob <- function(p) formatC(p, format = "f", digits = 3)
  res <- report_calibration(cal, 3L)
  expect_length(res, 5L)
  expect_identical(res[[1L]],
                   list(type = "heading",
                        text = "Calibration of the thresholds"))
  expect_identical(
    res[[2L]][["text"]],
    paste(sprintf(paste("The thresholds $\\gamma_{go} = %s$ and",
                        "$\\gamma_{nogo} = %s$ were chosen from 99 candidate",
                        "values from 0.01 to 0.99 (%s)."),
                  report_num(cal$gamma_go), report_num(cal$gamma_nogo),
                  cal$method),
          sprintf(paste("$\\gamma_{go}$ is the smallest candidate for which",
                        "the probability that the Go criterion is met is at",
                        "most 0.05 under the scenario $\\pi_t = 0.4$,",
                        "$\\pi_c = 0.4$ (the probability is %s), and",
                        "$\\gamma_{nogo}$ is the smallest candidate for which",
                        "the probability that the NoGo criterion is met is",
                        "at most 0.2 under the scenario $\\pi_t = 0.7$,",
                        "$\\pi_c = 0.4$ (the probability is %s)."),
                  prob(cal$prob_go), prob(cal$prob_nogo)))
  )
  table <- res[[4L]]
  expect_identical(table$header,
                   c("Scenario", "$\\pi_t$", "$\\pi_c$", "Go", "NoGo",
                     "Gray", "Miss"))
  expect_identical(table$rows[, 1L],
                   c("Go would be wrong", "NoGo would be wrong"))
  expect_identical(table$rows[, 4L], prob(cal$oc$Go))
  expect_identical(res[[5L]][["text"]], report_legend(cal$design))
})
