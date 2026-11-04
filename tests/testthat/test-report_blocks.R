test_that("report_blocks() describes a design or a specification", {
  des <- ref_design("D1")
  expect_identical(report_blocks(des, 3L), report_design(des))
  expect_identical(report_blocks(spec_d1(), 3L), report_design(des))
})

test_that("report_blocks() adds the result to the design", {
  des <- ref_design("D1")
  oc <- qdm_oc(des, scen_c)
  expect_identical(report_blocks(oc, 2L),
                   c(report_design(des), report_oc(oc, 2L)))
  cal <- run_cal("D1")
  expect_identical(report_blocks(cal, 3L),
                   c(report_design(cal$design), report_calibration(cal, 3L)))
  b <- qdm_boundary(des)
  expect_identical(report_blocks(b, 3L),
                   c(report_design(des), report_boundary(b)))
  dec <- qdm_decide(des, data_binary(12, 9), data_binary(12, 4))
  expect_identical(report_blocks(dec, 3L),
                   c(report_design(des), report_decision(dec, 3L)))
})

test_that("report_blocks() rejects other objects", {
  cal <- run_cal("D1")
  cal$design <- NULL
  expect_error(report_blocks(cal, 3L),
               "The calibration did not find both thresholds", fixed = TRUE)
  expect_error(report_blocks(list(), 3L),
               "'x' must be a design, a design specification, or a result of",
               fixed = TRUE)
})
