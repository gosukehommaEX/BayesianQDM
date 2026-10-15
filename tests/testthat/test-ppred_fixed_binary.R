# Reference: Python_reference/ref_fixed_tte_qdm.py (rational arithmetic)
test_that("ppred_fixed_binary() does not count a future effect equal to q", {
  # P(Y / 10 - 0.3 > 0) for Y ~ Beta-binomial(10, 4, 8) is P(Y >= 4) = 283/646;
  # counting Y = 3 would give 0.632684652808
  p <- ppred_fixed_binary(0, 10, 4, 8, 0.3, lower.tail = FALSE)
  expect_equal(p, 283 / 646, tolerance = 1e-12)
  expect_equal(ppred_fixed_binary(0, 10, 4, 8, 0.3), 1 - 283 / 646,
               tolerance = 1e-12)
  # A tie produced with rounding errors (0.1 + 0.2 is not exactly 0.3)
  expect_equal(ppred_fixed_binary(0.1, 10, 4, 8, 0.2, lower.tail = FALSE),
               283 / 646, tolerance = 1e-12)
})

test_that("ppred_fixed_binary() is vectorized over the Beta parameters", {
  ref <- c(0.0350877192982, 0.125902992776, 0.267801857585)
  y <- 0:2
  expect_equal(ppred_fixed_binary(0, 10, 1 + y, 11 - y, 0.3,
                                  lower.tail = FALSE), ref, tolerance = 1e-10)
  expect_equal(ppred_fixed_binary(-2, 10, c(1, 5), 3, 0.5,
                                  lower.tail = FALSE), c(1, 1))
  expect_equal(ppred_fixed_binary(2, 10, c(1, 5), 3, 0.5), c(1, 1))
})
