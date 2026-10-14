test_that("check_external_binary2() accepts NULL and data_binary2()", {
  expect_null(check_external_binary2(NULL, "external_t"))
  ext <- borrow(data_binary2(c(1, 2, 3, 4)), 0.5)
  expect_identical(check_external_binary2(ext, "external_t"), ext)
  expect_error(check_external_binary2(borrow(data_binary(10, 3), 0.5),
                                      "external_c"),
               paste("'external_c' must be NULL or created by borrow\\(\\)",
                     "with data_binary2\\(\\)."))
})
