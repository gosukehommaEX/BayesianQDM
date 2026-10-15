fake_spec <- function(values, discrete, go, nogo) {
  list(values = values, discrete = discrete,
       prob = function(idx, x) list(go = go(idx, x), nogo = nogo(idx, x)))
}
# Go if P(Go region) >= 0.5, NoGo if P(NoGo region) >= 0.8
crit_fake <- criterion_posterior(0.2, 0.1, 0.5, 0.8)
far <- function(x) sort(c(min(x) - c(1e4, 1e2), x, max(x) + c(1e2, 1e4)))

test_that("boundary_intervals() merges the values of a discrete grid", {
  s <- fake_spec(list(0:5, 0:3), TRUE,
                 go = function(i, x) pmin(1, x / c(5, 2)[i]),
                 nogo = function(i, x) pmax(0, 1 - x / 3))
  iv <- boundary_intervals(s, crit_fake)
  expect_identical(names(iv), c("row", "decision", "from", "to"))
  expect_identical(iv$row, c(1L, 1L, 1L, 2L, 2L))
  expect_identical(iv$decision, c("NoGo", "Gray", "Go", "NoGo", "Go"))
  expect_identical(iv$from, c(0L, 1L, 3L, 0L, 1L))
  expect_identical(iv$to, c(0L, 2L, 5L, 0L, 3L))
})

test_that("boundary_intervals() locates the boundaries of a continuous grid", {
  # Go if x >= 1; NoGo if x <= 3 - log(4), so both are met (Miss) between
  s <- fake_spec(list(far(seq(-5, 5, by = 1))), FALSE,
                 go = function(i, x) stats::plogis(x - 1),
                 nogo = function(i, x) stats::plogis(3 - x))
  iv <- boundary_intervals(s, crit_fake)
  expect_identical(iv$decision, c("NoGo", "Miss", "Go"))
  expect_equal(iv$from, c(-Inf, 1, 3 - log(4)), tolerance = 1e-9)
  expect_equal(iv$to, c(1, 3 - log(4), Inf), tolerance = 1e-9)
  # A Go region between two boundaries (not monotone)
  s2 <- fake_spec(list(far(seq(-3, 3, by = 0.5)), far(seq(-3, 3, by = 2))),
                  FALSE,
                  go = function(i, x) exp(-x ^ 2 / 2),
                  nogo = function(i, x) 0 * x)
  iv2 <- boundary_intervals(s2, crit_fake)
  b <- sqrt(2 * log(2))
  expect_identical(iv2$decision, rep(c("Gray", "Go", "Gray"), 2))
  expect_equal(iv2$from, rep(c(-Inf, -b, b), 2), tolerance = 1e-9)
  # A condition without any change of the decision
  s3 <- fake_spec(list(far(c(-1, 1))), FALSE,
                  go = function(i, x) 0 * x, nogo = function(i, x) 0 * x)
  iv3 <- boundary_intervals(s3, crit_fake)
  expect_identical(iv3$decision, "Gray")
  expect_identical(c(iv3$from, iv3$to), c(-Inf, Inf))
})

test_that("boundary_intervals() drops intervals of zero width", {
  # The Go and NoGo criteria change at the same point 0: the Miss decision
  # at that single point is not reported
  s <- fake_spec(list(c(-0.5, 0.5)), FALSE,
                 go = function(i, x) as.numeric(x >= 0),
                 nogo = function(i, x) as.numeric(x <= 0))
  iv <- boundary_intervals(s, crit_fake)
  expect_identical(iv$decision, c("NoGo", "Go"))
  expect_lt(abs(iv$to[1L]), 1e-11)
  expect_lt(abs(iv$from[2L]), 1e-11)
})
