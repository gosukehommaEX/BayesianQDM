test_that("check_color() accepts color names and hexadecimal codes", {
  expect_invisible(check_color("red", "col_go"))
  expect_identical(check_color("#658D1B", "col_go"), "#658D1B")
  expect_identical(check_color("#658D1B80", "col_go"), "#658D1B80")
  expect_identical(check_color("gray80", "col_gray"), "gray80")
})

test_that("check_color() rejects anything else", {
  msg <- "'col_go' must be a single color name or hexadecimal color code."
  expect_error(check_color("nocolor", "col_go"), msg, fixed = TRUE)
  expect_error(check_color("#12345", "col_go"), msg, fixed = TRUE)
  expect_error(check_color(c("red", "blue"), "col_go"), msg, fixed = TRUE)
  expect_error(check_color(NA_character_, "col_go"), msg, fixed = TRUE)
  expect_error(check_color(2, "col_go"), msg, fixed = TRUE)
  expect_error(check_color(character(0), "col_go"), msg, fixed = TRUE)
})
