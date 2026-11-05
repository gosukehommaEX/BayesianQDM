test_that("qdm_find_n() gives the probabilities of the candidates", {
  res <- find_n_d1(n_t = seq(6, 20, by = 2))
  expect_s3_class(res, "qdm_n_search")
  expect_identical(res$curve$n_t, seq(6L, 20L, by = 2L))
  expect_identical(res$curve$n_c, seq(6L, 20L, by = 2L))
  # Exact references (helper-find-n.R)
  expect_equal(res$curve$prob_go, ref_find_n$prob_go, tolerance = 1e-10)
  expect_equal(res$curve$prob_nogo, ref_find_n$prob_nogo, tolerance = 1e-10)
  expect_identical(res$curve$met,
                   c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE))
  # The smallest candidate that meets the targets (14) and the smallest from
  # which all larger candidates meet them (20)
  expect_identical(c(res$n_t, res$n_c, res$n_t_stable, res$n_c_stable),
                   c(14L, 14L, 20L, 20L))
  expect_identical(res$prob_go, res$curve$prob_go[5L])
  expect_identical(res$prob_nogo, res$curve$prob_nogo[5L])
  expect_identical(res$design, design_binary(14, 14, prior_beta(1, 1),
                                             prior_beta(1, 1), crit_post()))
  expect_identical(res$method, "exact enumeration of the outcomes")
  expect_identical(c(res$target_go, res$target_nogo), c(0.1, 0.1))
  expect_identical(res$scenario_go, cal_args$D1$scenario_go)
  expect_identical(res$scenario_nogo, cal_args$D1$scenario_nogo)
})

test_that("qdm_find_n() agrees with qdm_oc() and qdm_calibrate()", {
  a <- cal_args$D1
  scen <- rbind(a$scenario_go, a$scenario_nogo)
  # 14 meets the targets, so that no warning is given
  n <- c(12, 14)
  res <- find_n_d1(n_t = n)
  for (k in seq_along(n)) {
    des <- qdm_build(qdm_sample_size(ref_design("D1"), n[k], n[k]))
    oc <- qdm_oc(des, scen)
    expect_equal(res$curve$prob_go[k], oc$Go[1] + oc$Miss[1],
                 tolerance = 1e-12)
    expect_equal(res$curve$prob_nogo[k], oc$NoGo[2] + oc$Miss[2],
                 tolerance = 1e-12)
  }
  # With 12 subjects per group (design D1), the probabilities are those of
  # the calibration at the thresholds of D1 (0.5 and 0.8)
  cal <- run_cal("D1")
  expect_equal(res$curve$prob_go[1L],
               cal$curve$prob_go[which.min(abs(cal$curve$gamma - 0.5))],
               tolerance = 1e-12)
  expect_equal(res$curve$prob_nogo[1L],
               cal$curve$prob_nogo[which.min(abs(cal$curve$gamma - 0.8))],
               tolerance = 1e-12)
})

test_that("qdm_find_n() gives the decision probabilities of the choice", {
  res <- find_n_d1(n_t = seq(6, 20, by = 2))
  a <- cal_args$D1
  oc <- qdm_oc(res$design, rbind(a$scenario_go, a$scenario_nogo))
  expect_identical(res$oc$scenario, c("scenario_go", "scenario_nogo"))
  for (col in c("pi_t", "pi_c", "Go", "NoGo", "Gray", "Miss")) {
    expect_identical(res$oc[[col]], oc[[col]], label = col)
  }
})

test_that("qdm_find_n() sorts the candidates and takes the control sizes", {
  # With fewer control subjects, no candidate meets the targets
  expect_warning(res <- find_n_d1(n_t = c(20, 14, 16), n_c = c(10, 7, 8)),
                 "No candidate sample size gives", fixed = TRUE)
  expect_identical(res$curve$n_t, c(14L, 16L, 20L))
  expect_identical(res$curve$n_c, c(7L, 8L, 10L))
  expect_warning(res2 <- find_n_d1(n_t = c(14, 16), n_c = 10),
                 "No candidate sample size gives", fixed = TRUE)
  expect_identical(res2$curve$n_c, c(10L, 10L))
  # Without a concurrent control group, n_c is 0
  a <- cal_args$D3
  res3 <- suppressWarnings(
    qdm_find_n(ref_design("D3"), a$scenario_go, a$scenario_nogo,
               a$target_go, a$target_nogo, n_t = c(10, 20))
  )
  expect_identical(res3$curve$n_c, c(0L, 0L))
})

test_that("qdm_find_n() gives the stable sample size", {
  # All candidates meet the targets
  res <- find_n_d1(n_t = c(20, 24))
  expect_identical(c(res$n_t, res$n_t_stable), c(20L, 20L))
  # The largest candidate does not meet them
  res2 <- find_n_d1(n_t = c(14, 16))
  expect_identical(res2$n_t, 14L)
  expect_true(is.na(res2$n_t_stable))
  expect_true(is.na(res2$n_c_stable))
})

test_that("qdm_find_n() warns when no candidate meets the targets", {
  expect_warning(res <- find_n_d1(n_t = c(6, 8)),
                 paste("No candidate sample size gives a probability of",
                       "meeting the Go criterion of at most 0.1"),
                 fixed = TRUE)
  expect_true(is.na(res$n_t))
  expect_true(is.na(res$prob_go))
  expect_null(res$oc)
  expect_null(res$design)
})

test_that("qdm_find_n() uses the same random numbers for every candidate", {
  des <- ref_design_cont("C1", method = "MM")
  sg <- data.frame(mu_t = 0, sd_t = 25, mu_c = 0, sd_c = 25)
  sn <- data.frame(mu_t = 40, sd_t = 25, mu_c = 0, sd_c = 25)
  res <- suppressWarnings(qdm_find_n(des, sg, sn, 0.2, 0.2, n_t = c(8, 12),
                                     nsim = 200L, seed = 2))
  des12 <- qdm_build(qdm_sample_size(des, 12, 12))
  oc <- qdm_oc(des12, rbind(sg, sn), nsim = 200L, seed = 2)
  expect_equal(res$curve$prob_go[2L], oc$Go[1] + oc$Miss[1],
               tolerance = 1e-12)
  expect_equal(res$curve$prob_nogo[2L], oc$NoGo[2] + oc$Miss[2],
               tolerance = 1e-12)
  expect_identical(res$method, attr(oc, "method"))
})

test_that("qdm_find_n() rejects invalid arguments", {
  a <- cal_args$D1
  expect_error(find_n_d1(n_t = c(12, 12)),
               "'n_t' must not contain repeated values.", fixed = TRUE)
  expect_error(find_n_d1(n_t = c(12, 14), n_c = c(1, 2, 3)),
               "'n_c' must have length 1 or the length of 'n_t'.",
               fixed = TRUE)
  expect_error(find_n_d1(n_t = 0), "'n_t' must be at least 1.", fixed = TRUE)
  expect_error(find_n_d1(n_t = 12, n_c = -1), "'n_c' must be at least 0.",
               fixed = TRUE)
  expect_error(find_n_d1(n_t = 12, nsim = 0), "'nsim' must be at least 1",
               fixed = TRUE)
  expect_error(find_n_d1(n_t = 12, seed = 1.5),
               "'seed' must contain whole numbers only.", fixed = TRUE)
  expect_error(qdm_find_n(ref_design("D1"), a$scenario_go, a$scenario_nogo,
                          target_go = 1, target_nogo = 0.1, n_t = 12),
               "'target_go' must be greater than 0 and less than 1.",
               fixed = TRUE)
  expect_error(qdm_find_n(ref_design("D1"), a$scenario_go, a$scenario_nogo,
                          target_go = 0.1, target_nogo = 0, n_t = 12),
               "'target_nogo' must be greater than 0 and less than 1.",
               fixed = TRUE)
  expect_error(qdm_find_n(ref_design("D1"),
                          rbind(a$scenario_go, a$scenario_go),
                          a$scenario_nogo, 0.1, 0.1, n_t = 12),
               "'scenario_go' must have exactly one row.", fixed = TRUE)
})
