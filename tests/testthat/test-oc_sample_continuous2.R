test_that("oc_sample_continuous2() gives one column per scenario", {
  des <- ref_design_cont2("K1")
  sc <- scen_cont2("K1")
  set.seed(1)
  s <- oc_sample_continuous2(des, sc, 50L)
  expect_identical(dim(s$p_go), c(50L, 3L))
  expect_true(s$per_scenario)
  expect_equal(s$weights, matrix(1 / 50, 50, 3))
  expect_identical(s$method, paste("simulation of 50 studies; probabilities",
                                   "by moment-matching approximation"))
  # The first study of scenario 2, recomputed with decide_continuous2()
  set.seed(1)
  dt <- draws_continuous2(50L, 12)
  dc <- draws_continuous2(50L, 12)
  one <- function(d) lapply(d, `[`, 1L)
  st <- stats_continuous2(one(dt), 12, c(0, 1), c(70, 6), 0.4)
  sc2 <- stats_continuous2(one(dc), 12, c(-10, -1), c(70, 6), 0.4)
  res <- decide_continuous2(
    des,
    data_continuous2(12, st$mean[1L, ], matrix(st$cov[1L, c(1, 2, 2, 3)], 2)),
    data_continuous2(12, sc2$mean[1L, ], matrix(sc2$cov[1L, c(1, 2, 2, 3)], 2))
  )
  expect_equal(s$p_go[1L, 2L], res$prob_go, tolerance = 1e-12)
  expect_equal(s$p_nogo[1L, 2L], res$prob_nogo, tolerance = 1e-12)
  # Identical scenarios give identical probabilities (common random numbers)
  set.seed(2)
  s2 <- oc_sample_continuous2(des, sc[c(2, 2), ], 20L)
  expect_identical(s2$p_go[, 1L], s2$p_go[, 2L])
})

test_that("oc_sample_continuous2() handles marginal rules and no control", {
  set.seed(3)
  s <- oc_sample_continuous2(ref_design_cont2("K5"), scen_cont2("K5"), 20L)
  expect_length(s$p_go, 2L)
  expect_identical(dim(s$p_nogo[[2]]), c(20L, 3L))
  # Without a concurrent control group only the treatment group is simulated
  set.seed(4)
  oc_sample_continuous2(ref_design_cont2("K4"), scen_cont2("K4"), 20L)
  after <- stats::runif(1)
  set.seed(4)
  draws_continuous2(20L, 15)
  expect_identical(stats::runif(1), after)
  set.seed(5)
  s6 <- oc_sample_continuous2(ref_design_cont2("K1", "MC", 300L),
                              scen_cont2("K1")[1, ], 5L)
  expect_identical(dim(s6$p_go), c(5L, 1L))
  expect_match(s6$method, "Monte Carlo simulation \\(nMC = 300\\)")
})

test_that("oc_sample_continuous2() handles a fixed control", {
  des <- ref_design_cont2("KF1")
  sc <- scen_cont2("K4")
  set.seed(7)
  s <- oc_sample_continuous2(des, sc, 10L)
  expect_identical(s$method, paste("simulation of 10 studies; probabilities",
                                   "by exact computation (fixed control)"))
  set.seed(7)
  dt <- draws_continuous2(10L, 10)
  one <- function(d) lapply(d, `[`, 3L)
  st <- stats_continuous2(one(dt), 10, c(sc$mu_t1[2], sc$mu_t2[2]),
                          c(sc$sd_t1[2], sc$sd_t2[2]), sc$rho_t[2])
  res <- decide_continuous2(
    des,
    data_continuous2(10, st$mean[1L, ], matrix(st$cov[1L, c(1, 2, 2, 3)], 2)),
    NULL
  )
  expect_equal(s$p_go[3L, 2L], res$prob_go, tolerance = 1e-12)
  expect_equal(s$p_nogo[3L, 2L], res$prob_nogo, tolerance = 1e-12)
})
