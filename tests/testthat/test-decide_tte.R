test_that("decide_tte() matches the references", {
  for (name in c("T1", "T2", "T5", "T6")) {
    d <- ref_data_tte[[name]]
    res <- decide_tte(ref_design_tte(name), d$t, d$c)
    expect_equal(c(res$prob_go, res$prob_nogo), unname(ref_prob_tte[[name]]),
                 tolerance = 1e-10, info = name)
  }
  res <- decide_tte(ref_design_tte("T6"), ref_data_tte$T6$t, NULL)
  expect_equal(res$posterior_t, c(shape = 8.1, rate = 166))
  expect_equal(res$posterior_c, c(shape = 10.1, rate = 51))
  expect_identical(decide_tte(ref_design_tte("T1"), ref_data_tte$T1$t,
                              NULL)$posterior_c, c(value = 0.72))
  res2 <- decide_tte(ref_design_tte("T2"), ref_data_tte$T2$t,
                     ref_data_tte$T2$c)
  expect_equal(res2$posterior_c, c(shape = 9.1, rate = 36))
})

test_that("decide_tte() validates the data", {
  expect_error(decide_tte(ref_design_tte("T2"), data_binary(10, 2), NULL),
               "'data_t' must be created by data_tte\\(\\).")
  expect_error(decide_tte(ref_design_tte("T2"), data_tte(5, 40), NULL),
               "'data_c' must be created by data_tte\\(\\) for a design")
  expect_error(decide_tte(ref_design_tte("T1"), data_tte(5, 40),
                          data_tte(5, 40)),
               "'data_c' must be NULL for a design without a concurrent")
})
