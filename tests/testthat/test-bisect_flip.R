test_that("bisect_flip() locates the change of a logical function", {
  pred <- function(i, x) x >= c(0.3, -2)[i]
  expect_equal(bisect_flip(pred, c(1L, 2L), c(0, -5), c(1, 5)), c(0.3, -2),
               tolerance = 1e-10)
  # A function that changes from TRUE to FALSE
  expect_equal(bisect_flip(function(i, x) x <= 1.5, 1L, 0, 4), 1.5,
               tolerance = 1e-10)
  expect_identical(bisect_flip(pred, integer(0), numeric(0), numeric(0)),
                   numeric(0))
  # The midpoint of the final interval: one step from [0, 1] gives [0, 0.5]
  expect_identical(bisect_flip(function(i, x) x >= 0.3, 1L, 0, 1,
                               n_iter = 1L), 0.25)
})
