test_that("check_external_binary() accepts NULL and binary external data", {
  ext <- borrow(data_binary(20, 6), 0.5)
  expect_null(check_external_binary(NULL, "e"))
  expect_identical(check_external_binary(ext, "e"), ext)
  expect_error(check_external_binary(data_binary(20, 6), "ext"),
               "'ext' must be NULL or created by borrow\\(\\) with data_binary")
  expect_error(check_external_binary(borrow(data_binary2(c(1, 2, 3, 4)),
                                            0.5), "ext"),
               "'ext' must be NULL")
})
