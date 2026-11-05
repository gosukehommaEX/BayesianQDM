test_that("as.data.frame() gives the chosen thresholds of a calibration", {
  cal <- run_cal("D1")
  df <- as.data.frame(cal)
  expect_identical(df, data.frame(gamma_go = cal$gamma_go,
                                  prob_go = cal$prob_go, target_go = 0.05,
                                  gamma_nogo = cal$gamma_nogo,
                                  prob_nogo = cal$prob_nogo,
                                  target_nogo = 0.2))
  # Exact references (helper-calibration.R)
  expect_equal(c(df$gamma_go, df$gamma_nogo), c(0.7, 0.41))
  expect_equal(c(df$prob_go, df$prob_nogo),
               c(ref_cal$D1$prob_go, ref_cal$D1$prob_nogo),
               tolerance = 1e-12)
  expect_identical(rownames(as.data.frame(cal, row.names = "D1")), "D1")
})

test_that("as.data.frame() gives NA for thresholds that were not found", {
  a <- cal_args$D1
  cal <- suppressWarnings(
    qdm_calibrate(ref_design("D1"), a$scenario_go, a$scenario_nogo,
                  target_go = 1e-6, target_nogo = 0.2)
  )
  df <- as.data.frame(cal)
  expect_true(is.na(df$gamma_go))
  expect_true(is.na(df$prob_go))
  expect_equal(df$gamma_nogo, 0.41)
})
