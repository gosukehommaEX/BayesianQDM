# References: Python_reference/ref_tte_psurvdiff.py (integration over the
# treatment hazard, independent of the integration over the control hazard
# in the package, and a Monte Carlo check of the reference itself).
ref_psurv <- data.frame(
  q = c(0.1, -0.2, 0, 0.15, -0.05, 0.3, -0.6, 0.95, -0.95, 0, 0.05),
  t0 = c(3, 3, 3, 12, 12, 6, 1, 1, 1, 2, 5),
  a_t = c(5.1, 5.1, 0.5, 20.1, 20.1, 2, 0.3, 0.1, 50, 200, 1),
  b_t = c(41, 41, 10, 300, 300, 30, 0.5, 1, 10, 400, 1),
  a_c = c(3.1, 0.6, 0.8, 15.1, 15.1, 40, 0.2, 0.1, 0.1, 180, 1),
  b_c = c(20, 8, 12, 160, 160, 200, 0.4, 1, 1, 400, 1),
  p = c(0.385583382387744, 0.587302972430312, 0.608102285596848,
        0.407897841150297, 0.913904306893679, 0.703965564247968,
        0.869776929562393, 0.001098668828026, 0.237991538782146,
        0.152129549407913, 0.326187822705525)
)

test_that("psurvdiff() matches the independent references", {
  r <- ref_psurv
  up <- psurvdiff(r$q, r$t0, r$a_t, r$b_t, r$a_c, r$b_c, lower.tail = FALSE)
  expect_equal(up, r$p, tolerance = 1e-10)
  low <- psurvdiff(r$q, r$t0, r$a_t, r$b_t, r$a_c, r$b_c)
  expect_equal(low, 1 - r$p, tolerance = 1e-10)
})

test_that("psurvdiff() handles values outside (-1, 1) and recycles", {
  expect_identical(psurvdiff(c(-1.5, -1, 1, 2), 3, 2, 10, 3, 12,
                             lower.tail = FALSE), c(1, 1, 0, 0))
  expect_identical(psurvdiff(c(-1.5, 1, 2), 3, 2, 10, 3, 12), c(0, 1, 1))
  expect_length(psurvdiff(0.1, c(1, 2, 3), 2, 10, 3, 12), 3L)
  # Identical distributions: P(S_t - S_c > 0) = 0.5
  expect_equal(psurvdiff(0, 3, 4, 20, 4, 20), 0.5, tolerance = 1e-10)
})

test_that("psurvdiff() validates its arguments", {
  expect_error(psurvdiff(0, 0, 1, 1, 1, 1), "'milestone' must be greater")
  expect_error(psurvdiff(0, 1, -1, 1, 1, 1), "'shape_t' must be greater")
  expect_error(psurvdiff(0, 1, 1, 1, 1, 1, lower.tail = NA),
               "'lower.tail' must be TRUE or FALSE.")
  expect_error(psurvdiff(c(0, 1), 1, c(1, 2, 3), 1, 1, 1),
               "must have length 1 or 3.")
})

test_that("psurvdiff() handles long vectors", {
  # 1000 values reach the check for user interrupts in the C++ loop
  p1 <- psurvdiff(0.1, 3, 5.1, 41, 9.1, 36, lower.tail = FALSE)
  p <- psurvdiff(rep(0.1, 1000L), 3, 5.1, 41, 9.1, 36, lower.tail = FALSE)
  expect_length(p, 1000L)
  expect_identical(p, rep(p1, 1000L))
})
