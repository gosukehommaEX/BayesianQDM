test_that("print.qdm_data() prints the description and returns invisibly", {
  d <- data_binary(12, 7)
  expect_output(print(d), "Data: binary, n = 12, y = 7")
  expect_invisible(print(d))
})
