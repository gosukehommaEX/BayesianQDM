test_that("check_external_tte() accepts NULL and time-to-event data", {
  ext <- borrow(data_tte(4, 30), 0.5)
  expect_null(check_external_tte(NULL, "external_t"))
  expect_identical(check_external_tte(ext, "external_t"), ext)
  expect_error(check_external_tte(borrow(data_binary(10, 2), 0.5),
                                  "external_c"),
               "'external_c' must be NULL or created by borrow\\(\\) with data_tte")
  expect_error(check_external_tte(data_tte(4, 30), "external_t"),
               "'external_t' must be NULL or created by borrow")
})
