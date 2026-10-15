test_that("cellprob_scenarios() gives NA for infeasible correlations", {
  out <- cellprob_scenarios(c(0.3, 0.35), c(0.4, 0.15), c(0.3, 0.6))
  expect_equal(out[1L, ], cellprob_binary2(0.3, 0.4, 0.3))
  expect_true(all(is.na(out[2L, ])))
  expect_identical(colnames(out), c("00", "01", "10", "11"))
})
