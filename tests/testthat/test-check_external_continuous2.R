test_that("check_external_continuous2() accepts two continuous endpoints", {
  expect_null(check_external_continuous2(NULL, "e"))
  ext <- borrow(data_continuous2(10, c(1, 2), diag(2)), 0.5)
  expect_identical(check_external_continuous2(ext, "e"), ext)
  expect_error(check_external_continuous2(borrow(data_binary2(c(1, 2, 3, 4)),
                                                 0.5), "external_c"),
               paste("'external_c' must be NULL or created by borrow\\(\\)",
                     "with data_continuous2\\(\\)."))
  expect_error(check_external_continuous2(data_continuous2(10, c(1, 2),
                                                           diag(2)),
                                          "external_t"),
               "'external_t' must be NULL")
})
