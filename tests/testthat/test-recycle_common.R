test_that("recycle_common() recycles length-1 arguments", {
  out <- recycle_common(list(a = 1, b = c(2, 3, 4)))
  expect_identical(out, list(a = c(1, 1, 1), b = c(2, 3, 4)))
})

test_that("recycle_common() rejects incompatible lengths", {
  expect_error(recycle_common(list(a = 1:2, b = 1:3)),
               "'a' must have length 1 or 3.")
})

test_that("recycle_common() keeps equal lengths and lists every bad argument", {
  out <- recycle_common(list(a = 1, b = 2))
  expect_identical(out, list(a = 1, b = 2))
  expect_error(recycle_common(list(a = 1:2, b = 1:4, c = 1:3)),
               "'a', 'c' must have length 1 or 4.")
})
