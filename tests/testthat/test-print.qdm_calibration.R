test_that("print.qdm_calibration() prints the thresholds and the table", {
  cal <- run_cal("D3")
  out <- capture.output(res <- withVisible(print(cal)))
  expect_false(res$visible)
  expect_identical(res$value, cal)
  expect_match(out[1], "^Calibration of the decision thresholds")
  expect_identical(out[4], "Decision probabilities with the chosen thresholds")
  expect_match(out[5], "scenario +pi_t +Go +NoGo +Gray +Miss")
  expect_match(out[6], "^ +scenario_go ")
})

test_that("print.qdm_calibration() omits the table without thresholds", {
  a <- cal_args$D1
  cal <- suppressWarnings(
    qdm_calibrate(ref_design("D1"), a$scenario_go, a$scenario_nogo,
                  target_go = 1e-6, target_nogo = 0.2)
  )
  out <- capture.output(print(cal))
  expect_length(out, 3L)
})
