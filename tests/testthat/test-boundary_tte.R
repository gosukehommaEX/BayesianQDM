test_that("boundary_tte() matches the fixed-control reference", {
  # Reference: Python_reference/ref_fixed_tte_qdm.py, design T1
  b <- boundary_tte(ref_design_tte("T1"), 30, c(6, 3), c(150, 100), NULL)
  expect_identical(names(b$grid), c("events_t", "exposure_t"))
  expect_equal(b$grid$events_t, c(3, 6, 3, 6))
  expect_equal(b$grid$exposure_t, c(100, 100, 150, 150))
  expect_equal(b$p_go[4L], 0.686259165773, tolerance = 1e-10)
  expect_equal(b$p_nogo[4L], 0.00106102273767, tolerance = 1e-10)
  # Go: exp(-3 hazard) > 0.72 + 0.15; NoGo: exp(-3 hazard) < 0.72
  expect_equal(b$p_go, stats::pgamma(-log(0.87) / 3, 0.1 + b$grid$events_t,
                                     1 + b$grid$exposure_t))
  expect_equal(b$p_nogo, stats::pgamma(-log(0.72) / 3, 0.1 + b$grid$events_t,
                                       1 + b$grid$exposure_t,
                                       lower.tail = FALSE))
  expect_identical(b$axes$x, "exposure_t")
  expect_identical(b$axes$y, "events_t")
  expect_null(b$axes$caption)
  b0 <- boundary_tte(ref_design_tte("T1"), 30, NULL, 150, NULL)
  expect_equal(b0$grid$events_t, 0:30)
  b12 <- boundary_tte(ref_design_tte("T1"), 12, NULL, 150, NULL)
  expect_equal(b12$grid$events_t, 0:12)
})

test_that("boundary_tte() matches the references for random controls", {
  b <- boundary_tte(ref_design_tte("T2"), 20, 5, c(40, 60), data_tte(9, 35))
  expect_equal(c(b$p_go[1L], b$p_nogo[1L]), unname(ref_prob_tte$T2),
               tolerance = 1e-10)
  expect_match(b$axes$caption, "^Control data: ")
  for (name in c("T5", "T6")) {
    b <- boundary_tte(ref_design_tte(name), 30, 6, 150, NULL)
    expect_equal(c(b$p_go, b$p_nogo), unname(ref_prob_tte[[name]]),
                 tolerance = 1e-10, info = name)
  }
})

test_that("boundary_tte() estimates predictive probabilities", {
  # Reference: Python_reference/ref_fixed_tte_qdm.py, design T3 (Monte Carlo
  # estimate with 10000 simulated future studies)
  set.seed(11)
  b <- boundary_tte(ref_design_tte("T3"), 30, c(6, 20), 70, NULL)
  ref <- 0.917785289262
  expect_lt(abs(b$p_go[1L] - ref), 4.5 * sqrt(ref * (1 - ref) / 1e4) + 0.002)
  expect_equal(b$p_go + b$p_nogo, c(1, 1))
  expect_lt(b$p_go[2L], b$p_go[1L])
})

test_that("boundary_tte() validates its arguments", {
  des <- ref_design_tte("T1")
  expect_error(boundary_tte(des, 30, 1:3, NULL, NULL),
               "'exposure_t' must be given for a time-to-event endpoint.")
  expect_error(boundary_tte(des, 30, -1, 100, NULL),
               "'events_t' must be at least 0.")
  expect_error(boundary_tte(des, 30, 1, 0, NULL),
               "'exposure_t' must be greater than 0.")
  expect_error(boundary_tte(des, c(10, 20), NULL, 100, NULL),
               "'n_t' must have length 1.")
  expect_error(boundary_tte(ref_design_tte("T2"), 20, 1, 100,
                            data_binary(10, 2)),
               "'data_c' must be created by data_tte\\(\\).")
})
