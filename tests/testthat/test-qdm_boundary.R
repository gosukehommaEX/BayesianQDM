test_that("qdm_boundary() returns the decision at each point of the grid", {
  des <- ref_design("D8")
  b <- qdm_boundary(des)
  expect_identical(class(b), c("qdm_boundary", "data.frame"))
  expect_identical(names(b), c("n_t", "y_t", "prob_go", "prob_nogo",
                               "decision"))
  expect_identical(b$n_t, rep(20L, 21L))
  expect_identical(levels(b$decision), c("Go", "NoGo", "Gray", "Miss"))
  expected <- decision_label(meets_threshold(b$prob_go, 0.8),
                             meets_threshold(b$prob_nogo, 0.8))
  expect_identical(as.character(b$decision), expected)
  expect_identical(as.character(b$decision[12L]),
                   qdm_decide(des, data_binary(20, 11))$decision)
  expect_identical(attr(b, "design"), des)
  expect_named(attr(b, "axes"), c("x", "y", "xlab", "ylab", "caption"))
  expect_identical(rownames(b), as.character(seq_len(21L)))
})

test_that("qdm_boundary() shows a Miss decision", {
  # Beta(2, 2) posterior for y = 1 of 2: P(p > 0.4) = 0.648 and
  # P(p < 0.25) = 0.15625, so both criteria are met
  des <- design_binary(2, 0, prior_beta(1, 1),
                       criterion = crit_post(0.6, 0.15, tv = 0.2, mav = 0.05),
                       fixed_c = 0.2)
  b <- qdm_boundary(des)
  expect_identical(as.character(b$decision), c("NoGo", "Miss", "Go"))
})

test_that("qdm_boundary() handles all supported endpoints", {
  b <- qdm_boundary(ref_design_cont("C1"), mean_t = c(0, 30), sd_t = 25,
                    data_c = ref_data_cont$C1$c)
  expect_identical(names(b), c("n_t", "mean_t", "sd_t", "diff", "prob_go",
                               "prob_nogo", "decision"))
  dec <- qdm_decide(ref_design_cont("C1"), data_continuous(12, 30, 25),
                    ref_data_cont$C1$c)
  expect_equal(b$prob_go[2L], dec$prob_go)
  expect_identical(as.character(b$decision[2L]), dec$decision)
  b2 <- qdm_boundary(ref_design_cont2("K5"), mean_t = list(3, 1),
                     cov_t = cov2(75, 5.5, 0.4))
  expect_identical(names(b2), c("n_t", "mean_t1", "mean_t2", "prob_go1",
                                "prob_go2", "prob_nogo1", "prob_nogo2",
                                "decision"))
  dec2 <- qdm_decide(ref_design_cont2("K5"), ref_data_cont2$K5$t)
  expect_equal(c(b2$prob_go1, b2$prob_go2), unname(dec2$prob_go))
  expect_identical(as.character(b2$decision), dec2$decision)
  b_t <- qdm_boundary(ref_design_tte("T1"), n_t = 10, exposure_t = 100)
  expect_equal(b_t$events_t, 0:10)
  expect_identical(names(b_t), c("events_t", "exposure_t", "prob_go",
                                 "prob_nogo", "decision"))
  b3 <- qdm_boundary(ref_design_tte("T2"), events_t = 5, exposure_t = 40,
                     data_c = data_tte(9, 35))
  expect_equal(b3$prob_go, ref_prob_tte$T2[["go"]], tolerance = 1e-10)
  expect_identical(as.character(b3$decision),
                   qdm_decide(ref_design_tte("T2"), data_tte(5, 40),
                              data_tte(9, 35))$decision)
})

test_that("qdm_boundary() restores the random number generator", {
  des <- ref_design_cont("C1", method = "MC", nMC = 2000L)
  run <- function(seed) {
    qdm_boundary(des, mean_t = c(10, 30), sd_t = 25,
                 data_c = ref_data_cont$C1$c, seed = seed)
  }
  set.seed(1)
  state <- .Random.seed
  b1 <- run(5)
  expect_identical(.Random.seed, state)
  expect_identical(run(5), b1)
  run(NULL)
  expect_false(identical(.Random.seed, state))
  expect_error(run(1.5), "'seed' must contain whole numbers only.")
})

test_that("qdm_boundary() validates its arguments", {
  expect_error(qdm_boundary("design"),
               "'design' must be created by a design_\\*\\(\\) function")
  expect_error(qdm_boundary(ref_design_bin2("B1")),
               "qdm_boundary\\(\\) does not support two binary endpoints")
  expect_error(qdm_boundary(ref_design("D8"), mean_t = 1),
               "'mean_t' is not used for this endpoint.")
  expect_error(qdm_boundary(ref_design("D8"), sd_t = 1,
                            data_c = data_binary(10, 2)),
               "'sd_t', 'data_c' are not used for this endpoint.")
  expect_error(qdm_boundary(ref_design_tte("T1"), exposure_t = 100,
                            cov_t = diag(2)),
               "'cov_t' is not used for this endpoint.")
  expect_error(qdm_boundary(ref_design_cont("C1"), mean_t = 1, sd_t = 1),
               "'data_c' must be given for a design with a concurrent")
  expect_error(qdm_boundary(ref_design_cont("C3"), mean_t = 1, sd_t = 1,
                            data_c = ref_data_cont$C1$c),
               "'data_c' must be NULL for a design without a concurrent")
})
