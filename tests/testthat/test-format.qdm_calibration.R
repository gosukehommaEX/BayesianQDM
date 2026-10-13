test_that("format.qdm_calibration() reports the chosen thresholds", {
  cal <- run_cal("D1")
  expect_identical(
    format(cal),
    c(paste("Calibration of the decision thresholds (exact enumeration of",
            "the outcomes)"),
      paste("  gamma_go   = 0.7: P(Go criterion met | scenario_go) = 0.0454",
            "(target 0.05)"),
      paste("  gamma_nogo = 0.41: P(NoGo criterion met | scenario_nogo) =",
            "0.182 (target 0.2)"))
  )
})

test_that("format.qdm_calibration() reports thresholds that were not found", {
  a <- cal_args$D1
  cal <- suppressWarnings(
    qdm_calibrate(ref_design("D1"), a$scenario_go, a$scenario_nogo,
                  target_go = 1e-6, target_nogo = 0.2)
  )
  expect_identical(format(cal)[2],
                   "  gamma_go  : no value in the grid meets the target 1e-06")
})
