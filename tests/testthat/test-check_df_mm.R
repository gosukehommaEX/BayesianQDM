test_that("check_df_mm() requires more than 4 degrees of freedom for MM", {
  expect_null(check_df_mm(4.5, 30, "MM"))
  expect_null(check_df_mm(2, 1, "NI"))
  expect_null(check_df_mm(2, 1, "MC"))
  expect_error(check_df_mm(4, 30, "MM"),
               "more than 4 degrees of freedom \\(here 4 and 30\\)")
  expect_error(check_df_mm(11, 3.5, "MM"), "\\(here 11 and 3.5\\)")
})
