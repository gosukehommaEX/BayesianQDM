test_that("check_external_continuous() accepts NULL and continuous data", {
  expect_null(check_external_continuous(NULL, "external_t"))
  ext <- borrow(data_continuous(20, 1, 2), 0.5)
  expect_identical(check_external_continuous(ext, "external_t"), ext)
})

test_that("check_external_continuous() rejects other objects", {
  msg <- paste("'external_c' must be NULL or created by borrow\\(\\) with",
               "data_continuous\\(\\).")
  expect_error(check_external_continuous(borrow(data_binary(20, 8), 0.5),
                                         "external_c"), msg)
  expect_error(check_external_continuous(data_continuous(20, 1, 2),
                                         "external_c"), msg)
})
