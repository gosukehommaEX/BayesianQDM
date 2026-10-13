# Reference values: Python_reference/ref_distributions_mpmath.py (exact
# rational summation with mpmath, 40 digits), cross-checked with scipy.
ref_pbetabinomdiff <- data.frame(
  q   = c(0, 0.1, 0.05, -0.1, 0),
  m_t = c(20, 30, 20, 15, 1000),
  a_t = c(7.5, 3, 7.5, 1.5, 50.5),
  b_t = c(5.5, 2, 5.5, 3.5, 50.5),
  m_c = c(20, 20, 20, 25, 1000),
  a_c = c(4.5, 2, 4.5, 2.5, 40.5),
  b_c = c(8.5, 3, 8.5, 2.5, 60.5),
  upper = c(0.80686057751115994, 0.62263529749757634, 0.74700508465629342,
            0.36688900774350382, 0.91128008512545626)
)

test_that("pbetabinomdiff() matches exact upper tail references", {
  r <- ref_pbetabinomdiff
  for (i in seq_len(nrow(r))) {
    expect_equal(pbetabinomdiff(r$q[i], r$m_t[i], r$a_t[i], r$b_t[i],
                                r$m_c[i], r$a_c[i], r$b_c[i],
                                lower.tail = FALSE),
                 r$upper[i], tolerance = 1e-9, info = paste("case", i))
  }
})

test_that("pbetabinomdiff() lower tail is the complement of the upper tail", {
  r <- ref_pbetabinomdiff
  for (i in seq_len(nrow(r))) {
    expect_equal(pbetabinomdiff(r$q[i], r$m_t[i], r$a_t[i], r$b_t[i],
                                r$m_c[i], r$a_c[i], r$b_c[i]),
                 1 - r$upper[i], tolerance = 1e-9, info = paste("case", i))
  }
})

test_that("pbetabinomdiff() puts exact ties in the lower tail", {
  # With m_t = m_c = 20 and q = 0.05, the difference equals q exactly when
  # Y_t - Y_c = 1. Brute-force enumeration with integer arithmetic.
  k <- 0:20
  p_t <- dbetabinom(k, 20, 7.5, 5.5)
  p_c <- dbetabinom(k, 20, 4.5, 8.5)
  above <- outer(k, k, "-") > 1
  expected <- sum(outer(p_t, p_c) * above)
  expect_equal(pbetabinomdiff(0.05, 20, 7.5, 5.5, 20, 4.5, 8.5,
                              lower.tail = FALSE),
               expected, tolerance = 1e-12)
})

test_that("pbetabinomdiff() is vectorised over q and the shapes", {
  vec <- pbetabinomdiff(c(0, 0.1), 20, c(7.5, 6.5), 5.5, 20, 4.5, 8.5)
  one <- c(pbetabinomdiff(0, 20, 7.5, 5.5, 20, 4.5, 8.5),
           pbetabinomdiff(0.1, 20, 6.5, 5.5, 20, 4.5, 8.5))
  expect_identical(vec, one)
})

test_that("pbetabinomdiff() validates its arguments", {
  expect_error(pbetabinomdiff(0, 0, 1, 1, 10, 1, 1), "'m_t' must be at least")
  expect_error(pbetabinomdiff(0, 10, 1, 1, c(10, 20), 1, 1),
               "'m_c' must have length 1.")
  expect_error(pbetabinomdiff(0, 10, 1, -1, 10, 1, 1),
               "'b_t' must be greater than 0.")
})
