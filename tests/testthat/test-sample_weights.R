test_that("sample_weights() returns the weights of one scenario", {
  expect_identical(sample_weights(list(weights = cbind(1:2, 3:4)), 2L), 3:4)
  s <- list(weights_t = cbind(1:2, 3:4), weights_c = cbind(5:6, 7:8))
  expect_identical(sample_weights(s, 1L), list(t = 1:2, c = 5:6))
})
