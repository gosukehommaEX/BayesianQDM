test_that("dmultinom_counts() matches dmultinom()", {
  x <- outcomes_binary2(5)
  p <- c(0.1, 0.2, 0.3, 0.4)
  ref <- apply(x, 1L, stats::dmultinom, prob = p)
  expect_equal(dmultinom_counts(x, p), ref, tolerance = 1e-12)
  expect_equal(sum(dmultinom_counts(x, p)), 1, tolerance = 1e-12)
})

test_that("dmultinom_counts() handles cells with probability 0 and NA", {
  x <- outcomes_binary2(3)
  p <- c(0.5, 0, 0.5, 0)
  w <- dmultinom_counts(x, p)
  expect_equal(w, apply(x, 1L, stats::dmultinom, prob = p),
               tolerance = 1e-12)
  expect_true(all(w[x[, "01"] > 0 | x[, "11"] > 0] == 0))
  expect_identical(dmultinom_counts(x, rep(NA_real_, 4)),
                   rep(NA_real_, nrow(x)))
})
