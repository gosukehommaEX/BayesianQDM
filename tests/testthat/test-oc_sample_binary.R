test_that("oc_sample_binary() enumerates the outcomes with y_t fastest", {
  des <- ref_design("D2")
  s <- oc_sample_binary(des, data.frame(pi_t = c(0.3, 0.6), pi_c = 0.2))
  n_out <- (des$n_t + 1) * (des$n_c + 1)
  expect_length(s$p_go, n_out)
  expect_length(s$p_nogo, n_out)
  expect_identical(dim(s$weights), c(as.integer(n_out), 2L))
  expect_equal(colSums(s$weights), c(1, 1), tolerance = 1e-12)
  # Outcome (y_t, y_c) = (3, 2) is in row 3 + 2 * (n_t + 1) + 1
  k <- 3 + 2 * (des$n_t + 1) + 1
  expect_equal(s$weights[k, 2], stats::dbinom(3, 10, 0.6) *
                 stats::dbinom(2, 8, 0.2))
  post_t <- posterior_binary(des$prior_t, 10, 3, des$external_t)
  post_c <- posterior_binary(des$prior_c, 8, 2, des$external_c)
  expect_equal(s$p_go[k], pbetadiff(0.2, post_t$a, post_t$b, post_c$a,
                                    post_c$b, lower.tail = FALSE))
  expect_identical(s$method, "exact enumeration of the outcomes")
})

test_that("oc_sample_binary() uses a fixed control distribution without a concurrent control", {
  des <- ref_design("D3")
  s <- oc_sample_binary(des, data.frame(pi_t = 0.5))
  expect_length(s$p_go, des$n_t + 1)
  expect_identical(dim(s$weights), c(13L, 1L))
  expect_equal(s$p_go[10], pbetadiff(0.2, 10, 4, 5, 9, lower.tail = FALSE))
})

test_that("oc_sample_binary() handles a fixed control", {
  des <- ref_design("D8")
  s <- oc_sample_binary(des, data.frame(pi_t = c(0.2, 0.6)))
  expect_length(s$p_go, 21)
  expect_identical(dim(s$weights), c(21L, 2L))
  y <- 0:20
  expect_equal(s$p_go, stats::pbeta(0.4, 1 + y, 21 - y, lower.tail = FALSE))
  expect_equal(s$p_nogo, stats::pbeta(0.25, 1 + y, 21 - y))
  expect_equal(s$weights[, 2], stats::dbinom(y, 20, 0.6))
})
