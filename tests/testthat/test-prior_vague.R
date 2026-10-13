test_that("prior_vague() creates a vague prior object", {
  p <- prior_vague()
  expect_s3_class(p, c("qdm_prior_vague", "qdm_prior"), exact = TRUE)
  expect_length(p, 0L)
})
