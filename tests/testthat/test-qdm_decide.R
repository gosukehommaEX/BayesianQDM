# Exact references: Python_reference/ref_binary_qdm.py (rational arithmetic).
# In design D7 the probability of the Go region equals the threshold 0.5
# exactly (a tie), which meets the Go criterion.
ref_decide <- data.frame(
  design = c("D1", "D1", "D1", "D2", "D3",
             "D4", "D5", "D5", "D6", "D7"),
  n_t = c(12, 12, 12, 10, 12, 12, 10, 10, 8, 6),
  y_t = c(9, 6, 3, 7, 9, 7, 8, 4, 1, 3),
  n_c = c(12, 12, 12, 8, NA, NA, 10, 10, 8, 6),
  y_c = c(4, 5, 6, 2, NA, NA, 3, 4, 0, 3),
  a_t = c(10, 7, 4, 11, 10, 8, 9, 5, 2, 4),
  b_t = c(4, 7, 10, 6, 4, 6, 3, 7, 8, 4),
  a_c = c(5, 6, 7, 8, 5, 4, 4, 5, 1, 4),
  b_c = c(9, 8, 7, 15, 9, 8, 8, 7, 9, 4),
  prob_go = c(0.8207320756665553, 0.24509783589486162, 0.010407141009605608,
              0.7502162840687308, 0.8207320756665553, 0.5924430094360085,
              0.9402635338185402, 0.40597304686911645, 0.23621980701275141,
              0.5),
  prob_nogo = c(0.07292664584103187, 0.5578589165230211, 0.9606462168290966,
                0.09572468725900227, 0.07292664584103187, 0.22499768661993003,
                0.059736466181459745, 0.5940269531308836, 0.5205968501563143,
                0.34098075320453497),
  decision = c("Go", "Gray", "NoGo", "Go", "Go",
               "Go", "Go", "Gray", "Miss", "Go")
)

test_that("qdm_decide() matches the exact references", {
  r <- ref_decide
  for (i in seq_len(nrow(r))) {
    des <- ref_design(r$design[i])
    data_c <- if (is.na(r$n_c[i])) NULL else data_binary(r$n_c[i], r$y_c[i])
    res <- qdm_decide(des, data_binary(r$n_t[i], r$y_t[i]), data_c)
    info <- paste("case", i)
    expect_equal(res$posterior_t, c(a = r$a_t[i], b = r$b_t[i]), info = info)
    expect_equal(res$posterior_c, c(a = r$a_c[i], b = r$b_c[i]), info = info)
    expect_equal(res$prob_go, r$prob_go[i], tolerance = 1e-10, info = info)
    expect_equal(res$prob_nogo, r$prob_nogo[i], tolerance = 1e-10,
                 info = info)
    expect_identical(res$decision, r$decision[i], info = info)
  }
})

test_that("qdm_decide() matches high-precision references with Jeffreys priors", {
  # Python_reference/ref_distributions_mpmath.py: Beta(7.5, 5.5) versus
  # Beta(4.5, 8.5), the posteriors for 7/12 versus 4/12 responders
  jeff <- prior_beta(0.5, 0.5)
  d_post <- design_binary(12, 12, jeff, jeff, crit_post(tv = 0.2, mav = 0))
  res <- qdm_decide(d_post, data_binary(12, 7), data_binary(12, 4))
  expect_equal(res$prob_go, 0.57601763350543618, tolerance = 1e-10)
  expect_identical(res$decision, "Go")
  d_pred <- design_binary(12, 12, jeff, jeff,
                          criterion_predictive(0.05, m_t = 20,
                                               gamma_go = 0.7,
                                               gamma_nogo = 0.5))
  res <- qdm_decide(d_pred, data_binary(12, 7), data_binary(12, 4))
  expect_equal(res$prob_go, 0.74700508465629342, tolerance = 1e-12)
  expect_equal(res$prob_nogo, 1 - 0.74700508465629342, tolerance = 1e-12)
  expect_identical(res$decision, "Go")
})

test_that("qdm_decide() uses the sample sizes of the observed data", {
  des <- ref_design("D1")
  res <- qdm_decide(des, data_binary(20, 15), data_binary(5, 1))
  expect_equal(res$posterior_t, c(a = 16, b = 6))
  expect_equal(res$posterior_c, c(a = 2, b = 5))
  expect_equal(res$prob_go, pbetadiff(0.2, 16, 6, 2, 5, lower.tail = FALSE))
})

test_that("qdm_decide() returns the inputs and the criteria met", {
  des <- ref_design("D6")
  d_t <- data_binary(8, 1)
  d_c <- data_binary(8, 0)
  res <- qdm_decide(des, d_t, d_c)
  expect_s3_class(res, "qdm_decision", exact = TRUE)
  expect_identical(res$decision, "Miss")
  expect_true(res$met_go)
  expect_true(res$met_nogo)
  expect_identical(res$design, des)
  expect_identical(res$data_t, d_t)
  expect_identical(res$data_c, d_c)
  res_u <- qdm_decide(ref_design("D3"), data_binary(12, 2))
  expect_null(res_u$data_c)
  expect_identical(res_u$decision, "NoGo")
})

test_that("qdm_decide() validates the data", {
  des <- ref_design("D1")
  expect_error(qdm_decide(des, list(n = 12, y = 3), data_binary(12, 3)),
               "'data_t' must be created by data_binary().")
  expect_error(qdm_decide(des, data_binary(12, 3)),
               "'data_c' must be created by data_binary\\(\\) for a design")
  expect_error(qdm_decide(des, data_binary(12, 3),
                          data_continuous(12, 1, 1)),
               "'data_c' must be created by data_binary\\(\\) for a design")
  expect_error(qdm_decide(ref_design("D3"), data_binary(12, 3),
                          data_binary(12, 3)),
               "'data_c' must be NULL for a design without a concurrent")
  expect_error(qdm_decide("design", data_binary(12, 3)),
               "'design' must be created by a design_\\*\\(\\) function")
})
