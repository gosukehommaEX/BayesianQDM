test_that("check_direction() recycles a single direction", {
  expect_identical(check_direction("greater", 1L), "greater")
  expect_identical(check_direction("less", 2L), c("less", "less"))
  expect_identical(check_direction(c("greater", "less"), 2L),
                   c("greater", "less"))
})

test_that("check_direction() rejects invalid directions", {
  msg1 <- "'direction' must be \"greater\" or \"less\"."
  msg2 <- paste0("'direction' must be \"greater\" or \"less\", or one of",
                 " them for each of the two endpoints.")
  expect_error(check_direction("up", 1L), msg1, fixed = TRUE)
  expect_error(check_direction(c("greater", "less"), 1L), msg1, fixed = TRUE)
  expect_error(check_direction(NA_character_, 1L), msg1, fixed = TRUE)
  expect_error(check_direction(1, 2L), msg2, fixed = TRUE)
  expect_error(check_direction(c("less", "less", "less"), 2L), msg2,
               fixed = TRUE)
  expect_error(check_direction(character(0), 2L), msg2, fixed = TRUE)
})
