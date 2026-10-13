# The references ref_cal and the arguments cal_args are defined in
# helper-calibration.R.
test_that("qdm_calibrate() matches the exact references", {
  for (name in names(ref_cal)) {
    ref <- ref_cal[[name]]
    expect_silent(cal <- run_cal(name))
    expect_s3_class(cal, "qdm_calibration", exact = TRUE)
    expect_equal(cal$gamma_go, ref$gamma_go, info = name)
    expect_equal(cal$gamma_nogo, ref$gamma_nogo, info = name)
    expect_equal(cal$prob_go, ref$prob_go, tolerance = 1e-10, info = name)
    expect_equal(cal$prob_nogo, ref$prob_nogo, tolerance = 1e-10,
                 info = name)
    expect_equal(cal$curve$gamma, seq(0.01, 0.99, by = 0.01), info = name)
    expect_equal(cal$curve$prob_go, ref$curve_go, tolerance = 1e-10,
                 info = name)
    expect_equal(cal$curve$prob_nogo, ref$curve_nogo, tolerance = 1e-10,
                 info = name)
    got <- cal$oc[, c("Go", "NoGo", "Gray", "Miss")]
    expect_equal(got, ref$oc, tolerance = 1e-9, info = name)
    expect_identical(cal$oc$scenario, c("scenario_go", "scenario_nogo"))
  }
})

test_that("qdm_calibrate() chooses the smallest threshold that meets the target", {
  cal <- run_cal("D1")
  k <- match(cal$gamma_go, cal$curve$gamma)
  expect_lte(cal$curve$prob_go[k], 0.05)
  expect_gt(cal$curve$prob_go[k - 1L], 0.05)
  k <- match(cal$gamma_nogo, cal$curve$gamma)
  expect_lte(cal$curve$prob_nogo[k], 0.2)
  expect_gt(cal$curve$prob_nogo[k - 1L], 0.2)
  # The probabilities of meeting the criteria do not increase with gamma
  expect_true(all(diff(cal$curve$prob_go) <= 1e-12))
  expect_true(all(diff(cal$curve$prob_nogo) <= 1e-12))
})

test_that("qdm_calibrate() returns the design with the chosen thresholds", {
  cal <- run_cal("D1")
  des <- ref_design("D1")
  expect_identical(cal$design$criterion$gamma_go, cal$gamma_go)
  expect_identical(cal$design$criterion$gamma_nogo, cal$gamma_nogo)
  cal$design$criterion$gamma_go <- des$criterion$gamma_go
  cal$design$criterion$gamma_nogo <- des$criterion$gamma_nogo
  expect_identical(cal$design, des)
  # Its operating characteristics are those reported in 'oc'
  cal <- run_cal("D5")
  oc <- qdm_oc(cal$design, rbind(cal_args$D5$scenario_go,
                                 cal_args$D5$scenario_nogo))
  expect_equal(oc$Go, cal$oc$Go, tolerance = 1e-12)
  expect_equal(oc$NoGo, cal$oc$NoGo, tolerance = 1e-12)
  expect_identical(cal$scenario_go, cal_args$D5$scenario_go)
  expect_identical(cal$target_nogo, 0.2)
  expect_identical(cal$method, "exact enumeration of the outcomes")
})

test_that("qdm_calibrate() accepts a custom grid", {
  cal <- run_cal("D1", grid = c(0.9, 0.5, 0.7, 0.5, 0.3))
  expect_equal(cal$curve$gamma, c(0.3, 0.5, 0.7, 0.9))
  expect_identical(cal$gamma_go, 0.7)
  expect_equal(cal$curve$prob_go,
               ref_cal$D1$curve_go[c(30, 50, 70, 90)], tolerance = 1e-10)
})

test_that("qdm_calibrate() warns when no threshold meets a target", {
  a <- cal_args$D1
  expect_warning(
    cal <- qdm_calibrate(ref_design("D1"), a$scenario_go, a$scenario_nogo,
                         target_go = 1e-6, target_nogo = 0.2),
    "No value in 'grid' gives a probability of meeting the Go criterion"
  )
  expect_true(is.na(cal$gamma_go))
  expect_true(is.na(cal$prob_go))
  expect_null(cal$design)
  expect_null(cal$oc)
  expect_equal(cal$gamma_nogo, ref_cal$D1$gamma_nogo)
  expect_warning(
    cal <- qdm_calibrate(ref_design("D1"), a$scenario_go, a$scenario_nogo,
                         target_go = 0.5, target_nogo = 1e-6,
                         grid = c(0.1, 0.2)),
    "No value in 'grid' gives a probability of meeting the NoGo criterion"
  )
  expect_true(is.na(cal$gamma_nogo))
  expect_equal(cal$gamma_go, 0.2)
})

test_that("qdm_calibrate() warns when Miss decisions are possible", {
  # Loose targets give small thresholds with gamma_go + gamma_nogo <= 1
  a <- cal_args$D5
  expect_warning(
    cal <- qdm_calibrate(ref_design("D5"), a$scenario_go, a$scenario_nogo,
                         target_go = 0.9, target_nogo = 0.9),
    "both met \\(Miss\\)"
  )
  expect_lte(cal$gamma_go + cal$gamma_nogo, 1)
  expect_true(any(cal$oc$Miss > 0))
})

test_that("qdm_calibrate() validates its arguments", {
  des <- ref_design("D1")
  a <- cal_args$D1
  expect_error(qdm_calibrate(des, rbind(a$scenario_go, a$scenario_go),
                             a$scenario_nogo, 0.05, 0.2),
               "'scenario_go' must have exactly one row.")
  expect_error(qdm_calibrate(des, a$scenario_go, data.frame(pi_t = 0.7),
                             0.05, 0.2),
               "'scenario_nogo' must have the column\\(s\\) pi_t, pi_c.")
  expect_error(qdm_calibrate(des, a$scenario_go, a$scenario_nogo, 0, 0.2),
               "'target_go' must be greater than 0 and less than 1.")
  expect_error(qdm_calibrate(des, a$scenario_go, a$scenario_nogo, 0.05, 1),
               "'target_nogo' must be greater than 0 and less than 1.")
  expect_error(qdm_calibrate(des, a$scenario_go, a$scenario_nogo, 0.05, 0.2,
                             grid = c(0, 0.5)),
               "'grid' must be greater than 0 and less than 1.")
  expect_error(qdm_calibrate(des, a$scenario_go, a$scenario_nogo, 0.05, 0.2,
                             nsim = -1),
               "'nsim' must be at least 1.")
  expect_error(qdm_calibrate(des, a$scenario_go, a$scenario_nogo, 0.05, 0.2,
                             seed = "a"),
               "'seed' must be numeric.")
  expect_error(qdm_calibrate(NULL, a$scenario_go, a$scenario_nogo, 0.05,
                             0.2),
               "'design' must be created by")
})

test_that("qdm_calibrate() with a seed leaves the random number stream unchanged", {
  set.seed(8)
  before <- .Random.seed
  cal <- run_cal("D3", seed = 1)
  expect_identical(.Random.seed, before)
  expect_equal(cal$gamma_go, ref_cal$D3$gamma_go)
})
