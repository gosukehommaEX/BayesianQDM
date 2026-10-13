test_that("decide_binary() returns probabilities and posterior parameters", {
  res <- decide_binary(ref_design("D4"), data_binary(12, 7), NULL)
  expect_named(res, c("prob_go", "prob_nogo", "posterior_t", "posterior_c"))
  expect_equal(res$posterior_t, c(a = 8, b = 6))
  # Control: Beta(1, 1) with external data 6/20 at a0 = 0.5
  expect_equal(res$posterior_c, c(a = 4, b = 8))
  expect_equal(res$prob_go, pbetadiff(0.2, 8, 6, 4, 8, lower.tail = FALSE))
  expect_equal(res$prob_nogo, pbetadiff(0.1, 8, 6, 4, 8))
})
