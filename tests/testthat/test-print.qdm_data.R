test_that("print.qdm_data() prints the description and returns invisibly", {
  d <- data_binary(12, 7)
  expect_output(res <- withVisible(print(d)), "Data: binary, n = 12, y = 7")
  expect_false(res$visible)
  expect_identical(res$value, d)
})
