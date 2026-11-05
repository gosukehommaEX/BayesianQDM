test_that("report_title() depends on the object of the report", {
  expect_identical(report_title(ref_design("D1")), "Decision criteria")
  expect_identical(report_title(spec_d1()), "Decision criteria")
  expect_identical(report_title(structure(list(), class = "qdm_oc")),
                   "Decision criteria and operating characteristics")
  expect_identical(report_title(structure(list(), class = "qdm_calibration")),
                   "Decision criteria and calibration of the thresholds")
  expect_identical(report_title(structure(list(), class = "qdm_n_search")),
                   "Decision criteria and sample size")
  expect_identical(report_title(structure(list(), class = "qdm_boundary")),
                   "Decision criteria and decision boundaries")
  expect_identical(report_title(structure(list(), class = "qdm_decision")),
                   "Decision criteria and decision for the observed data")
})
