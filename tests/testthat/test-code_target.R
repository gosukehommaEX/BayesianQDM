test_that("code_target() gives the design or the specification", {
  des <- ref_design("D1")
  expect_identical(code_target(des), des)
  expect_identical(code_target(spec_d1()), spec_d1())
  expect_identical(code_target(qdm_oc(des, scen_c)), des)
  expect_identical(code_target(qdm_boundary(des)), des)
  dec <- qdm_decide(des, data_binary(12, 9), data_binary(12, 4))
  expect_identical(code_target(dec), des)
  cal <- run_cal("D1")
  expect_identical(code_target(cal), cal$design)
})

test_that("code_target() rejects objects without a design", {
  expect_error(code_target(list()),
               "'x' must be a design, a design specification, or a result",
               fixed = TRUE)
  none <- suppressWarnings(find_n_d1(n_t = c(6, 8)))
  expect_error(code_target(none), "'x' has no design", fixed = TRUE)
  a <- cal_args$D1
  cal <- suppressWarnings(
    qdm_calibrate(ref_design("D1"), a$scenario_go, a$scenario_nogo,
                  target_go = 1e-6, target_nogo = 0.2)
  )
  expect_error(code_target(cal), "'x' has no design", fixed = TRUE)
})
