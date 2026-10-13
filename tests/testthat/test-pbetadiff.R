# Reference values: Python_reference/ref_distributions_mpmath.py (mpmath, 40
# digits, two integral representations), cross-checked with scipy.
ref_pbetadiff <- data.frame(
  q   = c(0, 0.1, -0.2, 0.05, 0.2, 0.001, 0.95, -0.95),
  a_t = c(5, 6, 2.5, 0.5, 7.5, 30.5, 20, 1),
  b_t = c(5, 3, 10.5, 12.5, 5.5, 9970.5, 1, 20),
  a_c = c(5, 5, 4.5, 0.5, 4.5, 10.5, 1, 20),
  b_c = c(5, 8, 8.5, 12.5, 8.5, 9990.5, 20, 1),
  upper = c(0.5, 0.81832837711327888, 0.60911710328926547,
            0.16845257184576472, 0.57601763350543618, 0.94619967208045943,
            0.27051374743310446, 0.72948625256689554)
)

test_that("pbetadiff() matches high-precision upper tail references", {
  r <- ref_pbetadiff
  for (i in seq_len(nrow(r))) {
    expect_equal(pbetadiff(r$q[i], r$a_t[i], r$b_t[i], r$a_c[i], r$b_c[i],
                           lower.tail = FALSE),
                 r$upper[i], tolerance = 1e-9, info = paste("case", i))
  }
})

test_that("pbetadiff() lower tail is the complement of the upper tail", {
  r <- ref_pbetadiff
  lower <- pbetadiff(r$q, r$a_t, r$b_t, r$a_c, r$b_c)
  expect_equal(lower, 1 - r$upper, tolerance = 1e-9)
})

test_that("pbetadiff() is vectorised with recycling", {
  q <- c(-0.1, 0, 0.1)
  vec <- pbetadiff(q, 6, 3, 5, 8, lower.tail = FALSE)
  one <- vapply(q, pbetadiff, numeric(1), a_t = 6, b_t = 3, a_c = 5, b_c = 8,
                lower.tail = FALSE)
  expect_identical(vec, one)
  expect_true(all(diff(vec) < 0))
})

test_that("pbetadiff() returns 0 and 1 outside the support", {
  expect_identical(pbetadiff(1, 2, 3, 4, 5, lower.tail = FALSE), 0)
  expect_identical(pbetadiff(-1, 2, 3, 4, 5, lower.tail = FALSE), 1)
})

test_that("pbetadiff() validates its arguments", {
  expect_error(pbetadiff(0, 0, 1, 1, 1), "'a_t' must be greater than 0.")
  expect_error(pbetadiff(c(0, 0.1), c(1, 2, 3), 1, 1, 1),
               "'q' must have length 1 or 3.")
  expect_error(pbetadiff(0, 1, 1, 1, 1, lower.tail = NA),
               "'lower.tail' must be TRUE or FALSE.")
})
