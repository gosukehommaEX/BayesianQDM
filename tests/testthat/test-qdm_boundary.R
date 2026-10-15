# References: Python_reference/ref_boundary_qdm.py (decisions for every
# number of responders or events from closed forms and numerical
# integration; boundary values of observed means by root finding).

test_that("qdm_boundary() gives the decision intervals of the responders", {
  des <- ref_design("D8")
  b <- qdm_boundary(des)
  expect_identical(class(b), c("qdm_boundary", "data.frame"))
  expect_identical(names(b), c("n_t", "decision", "from", "to", "rate_from",
                               "rate_to", "effect_from", "effect_to"))
  expect_identical(levels(b$decision), c("Go", "NoGo", "Gray", "Miss"))
  expect_identical(as.character(b$decision), c("NoGo", "Gray", "Go"))
  expect_equal(b$from, c(0, 4, 10))
  expect_equal(b$to, c(3, 9, 20))
  expect_equal(b$rate_from, c(0, 0.2, 0.5))
  expect_equal(b$rate_to, c(0.15, 0.45, 1))
  expect_equal(b$effect_from, c(0, 0.2, 0.5) - 0.2)
  expect_identical(b$n_t, rep(20L, 3))
  expect_identical(rownames(b), as.character(1:3))
  expect_identical(attr(b, "design"), des)
  expect_true(attr(b, "discrete"))
  expect_true(attr(b, "increasing"))
  expect_null(attr(b, "range"))
  # The decisions agree with qdm_decide()
  expect_identical(qdm_decide(des, data_binary(20, 10))$decision, "Go")
  expect_identical(qdm_decide(des, data_binary(20, 9))$decision, "Gray")
  expect_identical(qdm_decide(des, data_binary(20, 3))$decision, "NoGo")
  # Several sample sizes
  b2 <- qdm_boundary(des, n_t = c(30, 10, 20))
  expect_identical(unique(b2$n_t), c(10L, 20L, 30L))
  expect_equal(b2$from[b2$decision == "Go"], c(6, 10, 15))
  expect_equal(b2$to[b2$decision == "NoGo"], c(1, 3, 5))
})

test_that("qdm_boundary() conditions on the number of control responders", {
  b <- qdm_boundary(ref_design("D1"))
  expect_identical(names(b)[1:3], c("n_t", "n_c", "y_c"))
  nogo_to <- c(0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12)
  go_from <- c(3:12, NA, NA, NA)
  for (yc in 0:12) {
    r <- b[b$y_c == yc, ]
    expect_equal(r$to[r$decision == "NoGo"], nogo_to[yc + 1L],
                 info = paste("y_c =", yc))
    expect_equal(r$from[r$decision == "NoGo"], 0)
    if (is.na(go_from[yc + 1L])) {
      expect_false(any(r$decision == "Go"), info = paste("y_c =", yc))
    } else {
      expect_equal(r$from[r$decision == "Go"], go_from[yc + 1L],
                   info = paste("y_c =", yc))
    }
  }
  # Observed difference in response rates at the smallest Go count
  r4 <- b[b$y_c == 4 & b$decision == "Go", ]
  expect_equal(r4$effect_from, 7 / 12 - 4 / 12)
})

test_that("qdm_boundary() shows a Miss decision", {
  # Beta(2, 2) posterior for y = 1 of 2: P(p > 0.4) = 0.648 and
  # P(p < 0.25) = 0.15625, so both criteria are met
  des <- design_binary(2, 0, prior_beta(1, 1),
                       criterion = crit_post(0.6, 0.15, tv = 0.2, mav = 0.05),
                       fixed_c = 0.2)
  b <- qdm_boundary(des)
  expect_identical(as.character(b$decision), c("NoGo", "Miss", "Go"))
  expect_equal(b$from, 0:2)
})

test_that("qdm_boundary() handles smaller effects being favorable", {
  # Rate of an adverse event in a single-arm study
  des <- design_binary(20, 0, prior_beta(1, 1),
                       criterion = criterion_posterior(0.1, 0.25, 0.8, 0.8,
                                                       direction = "less"),
                       fixed_c = 0)
  b <- qdm_boundary(des)
  expect_identical(as.character(b$decision), c("Go", "Gray", "NoGo"))
  expect_equal(b$from, c(0, 1, 7))
  expect_equal(b$to, c(0, 6, 20))
  expect_false(attr(b, "increasing"))
})

test_that("qdm_boundary() finds the boundary values of observed means", {
  b <- qdm_boundary(ref_design_cont("C1"), sd_t = c(35, 20, 25),
                    data_c = ref_data_cont$C1$c)
  expect_identical(names(b), c("n_t", "sd_t", "decision", "from", "to",
                               "effect_from", "effect_to"))
  expect_identical(unique(b$sd_t), c(20, 25, 35))
  expect_identical(as.character(b$decision), rep(c("NoGo", "Gray", "Go"), 3))
  go <- b[b$decision == "Go", ]
  nogo <- b[b$decision == "NoGo", ]
  expect_equal(go$from, c(20.4769259037, 20.7390764486, 21.3221524551),
               tolerance = 1e-7)
  expect_equal(nogo$to, c(-2.14198897276, -2.68579042435, -3.89724984505),
               tolerance = 1e-7)
  expect_identical(go$to, rep(Inf, 3))
  expect_identical(nogo$from, rep(-Inf, 3))
  expect_equal(go$effect_from, go$from + 12)
  # The boundaries separate the decisions of qdm_decide()
  dc <- ref_data_cont$C1$c
  expect_identical(qdm_decide(ref_design_cont("C1"),
                              data_continuous(12, 20.75, 25), dc)$decision,
                   "Go")
  expect_identical(qdm_decide(ref_design_cont("C1"),
                              data_continuous(12, 20.73, 25), dc)$decision,
                   "Gray")
  expect_length(attr(b, "range"), 2L)
})

test_that("qdm_boundary() matches the closed form of a predictive criterion", {
  # Design C8 (vague prior, fixed control 0, 30 future subjects): the mean of
  # the future subjects is t_14(ybar, s^2 (1 / 15 + 1 / 30)), so Go iff
  # ybar >= 5 + qt(0.7, 14) s sqrt(1 / 15 + 1 / 30) and NoGo iff
  # ybar <= 5 - qt(0.6, 14) s sqrt(1 / 15 + 1 / 30) (values by scipy)
  b <- qdm_boundary(ref_design_cont("C8"), sd_t = 10)
  expect_identical(as.character(b$decision), c("NoGo", "Gray", "Go"))
  expect_equal(b$from[3L], 6.696726971719765, tolerance = 1e-9)
  expect_equal(b$to[1L], 4.18345989303384, tolerance = 1e-9)
})

test_that("qdm_boundary() finds a bounded Go interval", {
  # With a Normal-Inverse-chi-squared prior centred at 10, the posterior
  # scale grows with the distance of the observed mean from 10, and the
  # probability of the Go region decreases again for large observed means
  des <- design_continuous(5, 0, prior_nix(10, 10, 3, 5),
                           criterion = criterion_posterior(2, 0, 0.99, 0.9),
                           fixed_c = 0)
  b <- qdm_boundary(des, sd_t = 5)
  expect_identical(as.character(b$decision), c("NoGo", "Gray", "Go", "Gray"))
  expect_equal(b$from[-1L], c(-90.0786294079, 1.82272184706, 61.9182489948),
               tolerance = 1e-7)
  expect_equal(b$to[-4L], c(-90.0786294079, 1.82272184706, 61.9182489948),
               tolerance = 1e-7)
})

test_that("qdm_boundary() handles two continuous endpoints", {
  b <- qdm_boundary(ref_design_cont2("KF1"),
                    mean_t = list(c(-20, 5, 20), seq(-10, 10, by = 0.5)),
                    cov_t = matrix(c(6400, 190, 190, 36), 2L))
  expect_identical(names(b), c("n_t", "mean_t1", "effect1", "decision",
                               "from", "to", "effect_from", "effect_to"))
  r <- function(m1, d) b[b$mean_t1 == m1 & b$decision == d, ]
  expect_equal(r(-20, "NoGo")$to, 0.756294245164, tolerance = 1e-7)
  expect_false(any(b$mean_t1 == -20 & b$decision == "Go"))
  expect_equal(r(5, "Go")$from, 2.74258829889, tolerance = 1e-7)
  expect_equal(r(5, "NoGo")$to, -1.10426503827, tolerance = 1e-7)
  expect_equal(r(20, "Go")$from, 1.47617391171, tolerance = 1e-7)
  expect_false(any(b$mean_t1 == 20 & b$decision == "NoGo"))
  expect_equal(r(5, "Go")$effect_from, 2.74258829889 + 1, tolerance = 1e-7)
})

test_that("qdm_boundary() gives the intervals of the events", {
  b <- qdm_boundary(ref_design_tte("T1"), n_t = 10,
                    exposure_t = c(150, 100))
  expect_identical(names(b), c("exposure_t", "decision", "from", "to",
                               "surv_from", "surv_to", "effect_from",
                               "effect_to"))
  expect_identical(as.character(b$decision), rep(c("Go", "Gray"), 2))
  expect_equal(b$to[b$decision == "Go"], c(4, 6))
  expect_equal(b$surv_to[1L], exp(-3 * 4 / 100))
  expect_false(attr(b, "increasing"))
  # Hazard ratio with a concurrent control
  bh <- qdm_boundary(ref_design_hr("H1"), exposure_t = 200,
                     data_c = ref_data_hr$c)
  expect_identical(names(bh), c("exposure_t", "decision", "from", "to",
                                "effect_from", "effect_to"))
  expect_identical(as.character(bh$decision), c("Go", "Gray", "NoGo"))
  expect_equal(bh$from, c(0, 15, 26))
  expect_equal(bh$to, c(14, 25, 30))
  expect_equal(bh$effect_to[1L], (14 / 200) / (20 / 180))
})

test_that("qdm_boundary() restores the random number generator", {
  des <- ref_design_cont("C1", method = "MC", nMC = 2000L)
  run <- function(seed) {
    qdm_boundary(des, sd_t = 25, mean_t = seq(-10, 30, by = 10),
                 data_c = ref_data_cont$C1$c, seed = seed)
  }
  set.seed(1)
  state <- .Random.seed
  b1 <- run(5)
  expect_identical(.Random.seed, state)
  expect_identical(run(5), b1)
  # Monte Carlo error makes the boundary values approximate
  go <- b1$from[b1$decision == "Go"]
  expect_length(go, 1L)
  expect_gt(go, 18)
  expect_lt(go, 23.5)
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
  expect_error(qdm_boundary(ref_design("D8"), sd_t = 1, cov_t = 1),
               "'sd_t', 'cov_t' are not used for this endpoint.")
  expect_error(qdm_boundary(ref_design_cont("C1"), sd_t = 25),
               "'data_c' must be given for a design with a concurrent")
  expect_error(qdm_boundary(ref_design_cont("C3"), sd_t = 25,
                            data_c = data_continuous(12, 0, 25)),
               "'data_c' must be NULL for a design without a concurrent")
})
