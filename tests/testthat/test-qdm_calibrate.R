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
               "'nsim' must be at least 1 and at most 2147483647.")
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

test_that("qdm_calibrate() matches the exact curves of a continuous design", {
  # Python_reference/ref_continuous_qdm.py, design C3 with method "MM":
  # probabilities that the Go criterion is met under mu_t = 10 and that the
  # NoGo criterion is met under mu_t = 20 (sd_t = 20), gamma = 0.1, ..., 0.9
  ex_go <- c(0.644148402974, 0.370239194278, 0.195251555187,
             0.0945586633576, 0.0416322583318, 0.0162867375315,
             0.00540389479526, 0.00136970089245, 0.000193194739424)
  ex_nogo <- c(0.33847199557, 0.123049979982, 0.0436537914186,
               0.0147880630576, 0.00468738422972, 0.00134979585112,
               0.000335015391109, 6.40013474872e-05, 6.72233598684e-06)
  nsim <- 20000
  grid <- seq(0.1, 0.9, by = 0.1)
  sc_go <- data.frame(mu_t = 10, sd_t = 20)
  sc_nogo <- data.frame(mu_t = 20, sd_t = 20)
  cal <- qdm_calibrate(ref_design_cont("C3", "MM"), sc_go, sc_nogo,
                       target_go = 0.1, target_nogo = 0.05, grid = grid,
                       nsim = nsim, seed = 5)
  tol <- function(p) 4.5 * sqrt(p * (1 - p) / nsim) + 0.002
  expect_true(all(abs(cal$curve$prob_go - ex_go) <= tol(ex_go)))
  expect_true(all(abs(cal$curve$prob_nogo - ex_nogo) <= tol(ex_nogo)))
  # The chosen thresholds are the smallest grid values whose simulated
  # probabilities meet the targets
  expect_equal(cal$gamma_go,
               grid[which(cal$curve$prob_go <= 0.1 + 1e-12)[1]])
  expect_equal(cal$gamma_nogo,
               grid[which(cal$curve$prob_nogo <= 0.05 + 1e-12)[1]])
  # The decision probabilities are those of qdm_oc() with the same seed
  oc <- qdm_oc(cal$design, rbind(sc_go, sc_nogo), nsim = nsim, seed = 5)
  expect_equal(cal$oc$Go, oc$Go)
  expect_equal(cal$oc$NoGo, oc$NoGo)
  expect_identical(cal$method, attr(oc, "method"))
})

test_that("qdm_calibrate() matches the exact curves for two binary endpoints", {
  # Python_reference/ref_binary2_qdm.py, design B1 (predictive, exact):
  # P(Go criterion met | scenario 1) and P(NoGo criterion met | scenario 2)
  ex_go <- c(0.593520788365099, 0.389565893365665, 0.240549158867837,
             0.159374861945474, 0.0931251408088504, 0.0268240062260297,
             0.0174331839434175, 0.00107659720288955, 0)
  ex_nogo <- c(0.6868067, 0.4085798, 0.2488433, 0.1210625, 0.0716579,
               0.0376649, 0.0119313, 0.0045342, 6.48e-05)
  grid <- seq(0.1, 0.9, by = 0.1)
  cal <- qdm_calibrate(ref_design_bin2("B1"), scen_bin2[1, ], scen_bin2[2, ],
                       target_go = 0.1, target_nogo = 0.05, grid = grid)
  expect_equal(cal$curve$prob_go, ex_go, tolerance = 1e-9)
  expect_equal(cal$curve$prob_nogo, ex_nogo, tolerance = 1e-9)
  expect_equal(cal$gamma_go, 0.5)
  expect_equal(cal$gamma_nogo, 0.6)
  expect_equal(cal$design$criterion$gamma_go, 0.5)
  expect_identical(cal$oc$Miss, c(0, 0))
})

test_that("qdm_calibrate() uses one threshold for both endpoints with a marginal rule", {
  des <- ref_design_bin2("B2")
  grid <- seq(0.2, 0.8, by = 0.1)
  cal <- qdm_calibrate(des, scen_bin2[1, ], scen_bin2[3, ], target_go = 0.1,
                       target_nogo = 0.1, grid = grid)
  expect_identical(length(cal$design$criterion$gamma_go), 2L)
  expect_equal(cal$design$criterion$gamma_go, rep(cal$gamma_go, 2))
  expect_equal(cal$design$criterion$gamma_nogo, rep(cal$gamma_nogo, 2))
  # The curve agrees with qdm_oc() for a design with the thresholds (g, g)
  g <- 0.4
  des_g <- des
  des_g$criterion$gamma_go <- c(g, g)
  oc <- suppressWarnings(qdm_oc(des_g, scen_bin2[1, ]))
  expect_equal(oc$Go + oc$Miss, cal$curve$prob_go[abs(grid - g) < 1e-9],
               tolerance = 1e-12)
  expect_true(all(diff(cal$curve$prob_go) <= 1e-12))
})

test_that("qdm_calibrate() requires feasible scenarios for two binary endpoints", {
  bad <- data.frame(pi_t1 = 0.35, pi_t2 = 0.15, rho_t = 0.6, pi_c1 = 0.2,
                    pi_c2 = 0.15, rho_c = 0.6)
  expect_error(qdm_calibrate(ref_design_bin2("B1"), bad, scen_bin2[2, ],
                             0.1, 0.1),
               "'scenario_go\\$rho_t' is outside the feasible range")
})

test_that("qdm_calibrate() matches the reference curves of design K1", {
  # Python_reference/ref_continuous2_qdm.py, section D2 (simulation of 1e5
  # studies per scenario)
  grid <- seq(0.05, 0.95, by = 0.05)
  ref_go <- c(0.56469, 0.38736, 0.27739, 0.20361, 0.15047, 0.11072, 0.08202,
              0.05972, 0.04304, 0.03054, 0.02119, 0.01411, 0.0089, 0.00552,
              0.00307, 0.00146, 0.00058, 0.00019, 0)
  ref_nogo <- c(0.7759, 0.63428, 0.52364, 0.43519, 0.36244, 0.30263, 0.25039,
                0.20524, 0.16722, 0.1336, 0.10499, 0.08139, 0.06218, 0.0451,
                0.031, 0.01991, 0.0115, 0.00557, 0.00147)
  sc <- scen_cont2("K1")
  # gamma_go + gamma_nogo < 1, so Miss decisions are possible
  cal <- withCallingHandlers(
    qdm_calibrate(ref_design_cont2("K1"), sc[1, ], sc[3, ], target_go = 0.05,
                  target_nogo = 0.2, grid = grid, nsim = 2000L, seed = 7),
    warning = function(w) {
      if (grepl("Miss", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  tol <- function(ref) 4.5 * sqrt(ref * (1 - ref) * (1 / 2000 + 1 / 1e5)) +
    0.002
  expect_true(all(abs(cal$curve$prob_go - ref_go) <= tol(ref_go)))
  expect_true(all(abs(cal$curve$prob_nogo - ref_nogo) <= tol(ref_nogo)))
  # The chosen thresholds meet the targets on the grid
  expect_lte(cal$prob_go, 0.05 + 1e-12)
  expect_lte(cal$prob_nogo, 0.2 + 1e-12)
  expect_s3_class(cal$design, "qdm_design_continuous2")
})
