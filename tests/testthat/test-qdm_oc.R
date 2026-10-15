# Exact references: Python_reference/ref_binary_qdm.py (rational arithmetic;
# the smallest distance between a probability of an outcome and a threshold
# is at least 3e-4, except for the intended ties of D7).
ref_oc <- list(
  D1 = data.frame(
    Go = c(0.012899294721250123, 0.14854199690661835,
           0.4890801931489531, 0.853069343654709),
    NoGo = c(0.8087444110220947, 0.41767794351065546,
             0.11426515535993896, 0.009738702385720657),
    Gray = c(0.1783562942566552, 0.4337800595827262,
             0.39665465149110796, 0.13719195395957035),
    Miss = c(0.0, 0.0,
             0.0, 0.0)
  ),
  D2 = data.frame(
    Go = c(0.008284858893991935, 0.13864155130259045,
           0.5111211700731249, 0.8841751124361872),
    NoGo = c(0.5935752005455708, 0.17640016449660517,
             0.023415970840182786, 0.000690132953071616),
    Gray = c(0.3981399405604372, 0.6849582842008044,
             0.4654628590866923, 0.11513475461074124),
    Miss = c(0.0, 0.0,
             0.0, 0.0)
  ),
  D3 = data.frame(
    Go = c(0.003903131648, 0.158212292608,
           0.665208557568, 0.980594720768),
    NoGo = c(0.79456894976, 0.22533728256,
             0.01526726656, 6.219776e-05),
    Gray = c(0.201527918592, 0.616450424832,
             0.319524175872, 0.019343081472),
    Miss = c(0.0, 0.0,
             0.0, 0.0)
  ),
  D4 = data.frame(
    Go = c(0.003903131648, 0.158212292608,
           0.665208557568, 0.980594720768),
    NoGo = c(0.55834574848, 0.08344332288,
             0.00281018368, 4.52608e-06),
    Gray = c(0.437751119872, 0.758344384512,
             0.331981258752, 0.019400753152),
    Miss = c(0.0, 0.0,
             0.0, 0.0)
  ),
  D5 = data.frame(
    Go = c(0.002175177111770235, 0.02824424826561626,
           0.1350357734374716, 0.42590770164544833),
    NoGo = c(0.776219456280734, 0.41002224966905626,
             0.12752251421483599, 0.014808410228795638),
    Gray = c(0.22160536660749583, 0.5617335020653275,
             0.7374417123476924, 0.559283888125756),
    Miss = c(0.0, 0.0,
             0.0, 0.0)
  ),
  D6 = data.frame(
    Go = c(0.0290828790595584, 0.1295783505690624,
           0.3524185448054784, 0.7117744587669504),
    NoGo = c(0.88082300796928, 0.600355270721536,
             0.283936647282688, 0.069097745711104),
    Gray = c(0.0, 0.0,
             0.0, 0.0),
    Miss = c(0.0900941129711616, 0.2700663787094016,
             0.3636448079118336, 0.2191277955219456)
  ),
  D7 = data.frame(
    Go = c(0.321296601088, 0.615339200512,
           0.841787707392, 0.965405052928),
    NoGo = c(0.20031995904, 0.06857588736,
             0.01526726656, 0.00126910464),
    Gray = c(0.478383439872, 0.316084912128,
             0.142945026048, 0.033325842432),
    Miss = c(0.0, 0.0,
             0.0, 0.0)
  )
)

test_that("qdm_oc() matches the exact references for all control forms", {
  for (name in names(ref_oc)) {
    des <- ref_design(name)
    scen <- if (des$n_c > 0) scen_c else scen_u
    if (name == "D6") {
      expect_warning(oc <- qdm_oc(des, scen), "both met \\(Miss\\)")
    } else {
      expect_silent(oc <- qdm_oc(des, scen))
    }
    got <- as.data.frame(oc)[, c("Go", "NoGo", "Gray", "Miss")]
    attr(got, "design") <- NULL
    attr(got, "method") <- NULL
    class(got) <- "data.frame"
    expect_equal(got, ref_oc[[name]], tolerance = 1e-9, info = name)
  }
})

test_that("qdm_oc() probabilities add up to 1", {
  oc <- qdm_oc(ref_design("D2"), data.frame(pi_t = c(0, 0.33, 1),
                                            pi_c = c(1, 0.5, 0)))
  expect_equal(oc$Go + oc$NoGo + oc$Gray + oc$Miss, rep(1, 3),
               tolerance = 1e-12)
})

test_that("qdm_oc() agrees with decisions for each outcome", {
  # Independent route: qdm_decide() for every outcome and binomial weights
  for (name in c("D2", "D3", "D5")) {
    des <- ref_design(name)
    pi_t <- 0.55
    pi_c <- 0.35
    y_c_all <- if (des$n_c > 0) 0:des$n_c else NA
    tab <- expand.grid(y_t = 0:des$n_t, y_c = y_c_all)
    tab$decision <- vapply(seq_len(nrow(tab)), function(k) {
      data_c <- if (des$n_c > 0) data_binary(des$n_c, tab$y_c[k]) else NULL
      qdm_decide(des, data_binary(des$n_t, tab$y_t[k]), data_c)$decision
    }, character(1))
    tab$w <- stats::dbinom(tab$y_t, des$n_t, pi_t) *
      (if (des$n_c > 0) stats::dbinom(tab$y_c, des$n_c, pi_c) else 1)
    oc <- qdm_oc(des, data.frame(pi_t = pi_t, pi_c = pi_c))
    for (d in c("Go", "NoGo", "Gray", "Miss")) {
      expect_equal(oc[[d]], sum(tab$w[tab$decision == d]), tolerance = 1e-12,
                   info = paste(name, d))
    }
  }
})

test_that("qdm_oc() keeps extra columns and records the design", {
  des <- ref_design("D1")
  scen <- data.frame(label = c("null", "alt"), pi_t = c(0.4, 0.7),
                     pi_c = 0.4, extra = 1:2)
  oc <- qdm_oc(des, scen)
  expect_s3_class(oc, c("qdm_oc", "data.frame"), exact = TRUE)
  expect_identical(names(oc), c("label", "pi_t", "pi_c", "extra", "Go",
                                "NoGo", "Gray", "Miss"))
  expect_identical(oc$label, c("null", "alt"))
  expect_identical(attr(oc, "design"), des)
  expect_identical(attr(oc, "method"), "exact enumeration of the outcomes")
  # The control rate is not needed without a concurrent control group
  oc_u <- qdm_oc(ref_design("D3"), data.frame(pi_t = 0.5, pi_c = 0.9))
  expect_equal(oc_u$Go, qdm_oc(ref_design("D3"),
                               data.frame(pi_t = 0.5))$Go)
})

test_that("qdm_oc() probability of Go increases with the treatment rate", {
  oc <- qdm_oc(ref_design("D1"),
               data.frame(pi_t = seq(0, 1, by = 0.05), pi_c = 0.3))
  expect_true(all(diff(oc$Go) >= -1e-12))
  expect_true(all(diff(oc$NoGo) <= 1e-12))
  expect_equal(oc$NoGo[1], 1 - oc$Gray[1] - oc$Go[1], tolerance = 1e-12)
})

test_that("qdm_oc() leaves the random number stream unchanged", {
  des <- ref_design("D1")
  set.seed(3)
  before <- .Random.seed
  a <- qdm_oc(des, scen_c, nsim = 10, seed = 99)
  expect_identical(.Random.seed, before)
  b <- qdm_oc(des, scen_c)
  expect_identical(a, b)
})

test_that("qdm_oc() validates its arguments", {
  des <- ref_design("D1")
  expect_error(qdm_oc(list(), scen_c), "'design' must be created by")
  expect_error(qdm_oc(des, c(pi_t = 0.4, pi_c = 0.4)),
               "'scenarios' must be a data frame with at least one row.")
  expect_error(qdm_oc(des, data.frame(pi_t = 0.4)),
               "'scenarios' must have the column\\(s\\) pi_t, pi_c.")
  expect_error(qdm_oc(des, data.frame(pi_t = 1.2, pi_c = 0.4)),
               "'scenarios\\$pi_t' must be at least 0 and at most 1.")
  expect_error(qdm_oc(des, data.frame(pi_t = 0.4, pi_c = 0.4, Go = 1)),
               "must not have columns named Go, NoGo, Gray or Miss.")
  expect_error(qdm_oc(des, scen_c, nsim = 0),
               "'nsim' must be at least 1 and at most 2147483647.")
  expect_error(qdm_oc(des, scen_c, nsim = 3e9),
               "'nsim' must be at least 1 and at most 2147483647.")
  expect_error(qdm_oc(des, scen_c, seed = 1.5),
               "'seed' must contain whole numbers only.")
  expect_error(qdm_oc(des, scen_c, seed = 1e10),
               "'seed' must be at least -2147483647 and at most 2147483647.")
})

test_that("qdm_oc() agrees with simulated trials", {
  skip_on_cran()
  des <- ref_design("D2")
  set.seed(42)
  nrep <- 4000
  y_t <- stats::rbinom(nrep, des$n_t, 0.6)
  y_c <- stats::rbinom(nrep, des$n_c, 0.4)
  dec <- vapply(seq_len(nrep), function(k) {
    qdm_decide(des, data_binary(des$n_t, y_t[k]),
               data_binary(des$n_c, y_c[k]))$decision
  }, character(1))
  oc <- qdm_oc(des, data.frame(pi_t = 0.6, pi_c = 0.4))
  for (d in c("Go", "NoGo", "Gray")) {
    p <- oc[[d]]
    expect_lt(abs(mean(dec == d) - p), 4 * sqrt(p * (1 - p) / nrep) + 1e-12)
  }
})

test_that("qdm_oc() reproduces Kang et al. (2026), Supplemental Table 2", {
  # Beta(0.5, 0.5) priors, n_t = n_c = 12, pi_c = 0.4 and treatment effects
  # -0.4, -0.2, 0, 0.4, 0.6. Published values (three decimals):
  # inst/validation/kang2026_supplement.csv, rows 1, 10, 19, 28 and 37 of
  # Table 2. The full comparison is in inst/validation/run_validation.R.
  jeff <- prior_beta(0.5, 0.5)
  scen <- data.frame(pi_t = c(0, 0.2, 0.4, 0.8, 1), pi_c = 0.4)
  d_post <- design_binary(12, 12, jeff, jeff, crit_post())
  d_pred <- design_binary(12, 12, jeff, jeff,
                          criterion_predictive(0, m_t = 50, gamma_go = 0.9,
                                               gamma_nogo = 0.7))
  pub_post <- rbind(c(0, 1, 0), c(0.013, 0.812, 0.175),
                    c(0.149, 0.418, 0.433), c(0.853, 0.010, 0.137),
                    c(0.997, 0, 0.003))
  pub_pred <- rbind(c(0, 0.998, 0.002), c(0.006, 0.677, 0.317),
                    c(0.076, 0.267, 0.657), c(0.739, 0.003, 0.258),
                    c(0.997, 0, 0.003))
  cols <- c("Go", "NoGo", "Gray")
  post <- as.matrix(qdm_oc(d_post, scen)[, cols])
  pred <- as.matrix(qdm_oc(d_pred, scen)[, cols])
  expect_lt(max(abs(post - pub_post)), 0.0005 + 1e-9)
  expect_lt(max(abs(pred - pub_pred)), 0.0005 + 1e-9)
})

# Continuous endpoint. Exact operating characteristics of the uncontrolled
# designs C3 and C6 (vague prior, hypothetical control):
# Python_reference/ref_continuous_qdm.py, where each criterion is a threshold
# on a statistic with a noncentral t distribution. The simulated
# probabilities are compared with a tolerance of 4.5 standard errors plus
# 0.002 (Python_reference/port_continuous_oc.py checks the same procedure).
exact_oc_cont <- list(
  C3 = list(
    scen = data.frame(mu_t = c(0, 20, 35), sd_t = c(20, 20, 25)),
    NI = list(Go = c(5.59007059741e-05, 0.339921868611, 0.951288822042),
              NoGo = c(0.502398874545, 0.000333627271368, 3.94025694328e-07)),
    MM = list(Go = c(5.60512055178e-05, 0.340178784483, 0.9513624542),
              NoGo = c(0.502915282701, 0.000335015391109, 3.96193102629e-07))
  ),
  C6 = list(
    scen = data.frame(mu_t = c(0, 5, 12), sd_t = c(12, 12, 15)),
    NI = list(Go = c(0.00299458547153, 0.12035425217, 0.722001478),
              NoGo = c(0.847885721764, 0.282128401915, 0.00881718335511)),
    MM = list(Go = c(0.00300601509893, 0.120621873617, 0.722483043333),
              NoGo = c(0.848072473479, 0.282387849867, 0.00883485467116))
  )
)

expect_oc_close <- function(oc, ref, nsim, info) {
  for (d in c("Go", "NoGo")) {
    p <- ref[[d]]
    expect_true(all(abs(oc[[d]] - p) <= 4.5 * sqrt(p * (1 - p) / nsim) +
                      0.002), info = paste(info, d))
  }
  expect_identical(oc$Miss, c(0, 0, 0), info = info)
}

test_that("qdm_oc() matches the exact references for continuous designs (MM)", {
  for (name in c("C3", "C6")) {
    ex <- exact_oc_cont[[name]]
    oc <- qdm_oc(ref_design_cont(name, "MM"), ex$scen, nsim = 20000,
                 seed = 20261013)
    expect_oc_close(oc, ex$MM, 20000, paste(name, "MM"))
  }
})

test_that("qdm_oc() matches the exact references for design C3 (NI)", {
  skip_on_cran()
  ex <- exact_oc_cont$C3
  oc <- qdm_oc(ref_design_cont("C3"), ex$scen, nsim = 1000, seed = 1)
  expect_oc_close(oc, ex$NI, 1000, "C3 NI")
})

test_that("qdm_oc() matches the exact references for design C6 (NI)", {
  skip_on_cran()
  ex <- exact_oc_cont$C6
  oc <- qdm_oc(ref_design_cont("C6"), ex$scen, nsim = 1000, seed = 2)
  expect_oc_close(oc, ex$NI, 1000, "C6 NI")
})

test_that("qdm_oc() simulates continuous designs with a concurrent control", {
  des <- ref_design_cont("C1", "MM")
  set.seed(5)
  before <- .Random.seed
  oc <- qdm_oc(des, cbind(label = c("a", "b", "c"), scen_cont_c),
               nsim = 500, seed = 7)
  expect_identical(.Random.seed, before)
  expect_s3_class(oc, c("qdm_oc", "data.frame"), exact = TRUE)
  expect_identical(names(oc), c("label", "mu_t", "sd_t", "mu_c", "sd_c",
                                "Go", "NoGo", "Gray", "Miss"))
  expect_equal(oc$Go + oc$NoGo + oc$Gray + oc$Miss, c(1, 1, 1))
  expect_identical(attr(oc, "method"),
                   paste("simulation of 500 studies; probabilities by",
                         "moment-matching approximation"))
  # Common random numbers: the probability of Go increases with the effect
  expect_true(all(diff(oc$Go) > 0))
  expect_true(all(diff(oc$NoGo) < 0))
  # The same seed gives the same result
  expect_identical(qdm_oc(des, scen_cont_c, nsim = 500, seed = 7)$Go, oc$Go)
})

test_that("qdm_oc() works with Monte Carlo probabilities", {
  # With the predictive criterion, the Monte Carlo probabilities of the Go
  # and NoGo regions of each simulated study add up to 1, so both criteria
  # cannot be met when gamma_go + gamma_nogo > 1 (0.7 + 0.6 in design C6)
  des <- ref_design_cont("C6", "MC", 2000L)
  oc <- qdm_oc(des, exact_oc_cont$C6$scen, nsim = 100, seed = 3)
  expect_identical(oc$Miss, c(0, 0, 0))
  expect_equal(oc$Go + oc$NoGo + oc$Gray, c(1, 1, 1))
  expect_lt(abs(oc$Go[3] - 0.722), 4.5 * sqrt(0.722 * 0.278 / 100) + 0.02)
})

test_that("qdm_oc() validates the scenarios of continuous designs", {
  des <- ref_design_cont("C1", "MM")
  expect_error(qdm_oc(des, data.frame(mu_t = 1, sd_t = 1)),
               paste("'scenarios' must have the column\\(s\\) mu_t, sd_t,",
                     "mu_c, sd_c."))
  expect_error(qdm_oc(des, data.frame(mu_t = 1, sd_t = 0, mu_c = 0,
                                      sd_c = 1)),
               "'scenarios\\$sd_t' must be greater than 0.")
  expect_error(qdm_oc(ref_design_cont("C3", "MM"), data.frame(mu_t = 1)),
               "'scenarios' must have the column\\(s\\) mu_t, sd_t.")
})

# Two binary endpoints. Exact references: Python_reference/ref_binary2_qdm.py
# (B: predictive criteria, rational arithmetic; D: posterior criteria on one
# endpoint at a time, differences of Beta variables). For designs B3 and B4
# no outcome has a probability within 5 Monte Carlo standard errors
# (nMC = 20000) of a threshold, so the Monte Carlo error does not change any
# decision. The published tables of Yamaguchi et al. (submitted) need the
# enumeration of 286 x 286 pairs of outcomes and are reproduced in
# inst/validation/validate_yamaguchi_binary2.R instead.
ref_oc_bin2 <- list(
  B1 = rbind(c(0.156029909046664, 0.558146336369747, 0.282478801684779,
               0.00334495289881015),
             c(0.3312351, 0.2440805, 0.4199216, 0.0047628),
             c(0.666686263299421, 0.0789967306310707, 0.249597442622606,
               0.00471956344690299)),
  B2 = rbind(c(0.0743517888004001, 0.326426105367674, 0.599222105831926, 0),
             c(0.2607633, 0.1627937, 0.576443, 0),
             c(0.535900183423335, 0.0307251043879517, 0.433374712188714, 0)),
  B3 = rbind(c(0.12022101, 0.406285269999998, 0.473493719999998, 0),
             c(0.426958559999999, 0.11187472, 0.461166719999999, 0),
             c(0.796917759999999, 0.01171712, 0.19136512, 0)),
  B4 = rbind(c(0.121502195183722, 0.424358895183721, 0.454138909632553, 0),
             c(0.465942411080158, 0.163130761080159, 0.370926827839682, 0),
             c(0.766733348737193, 0.0171192987371942, 0.216147352525611, 0))
)

test_that("qdm_oc() matches the exact references for two binary endpoints", {
  for (name in names(ref_oc_bin2)) {
    scen <- if (name %in% c("B1", "B2")) scen_bin2 else scen_bin2_d
    if (name == "B1") {
      expect_warning(oc <- qdm_oc(ref_design_bin2(name), scen, seed = 1),
                     "Miss")
    } else {
      oc <- qdm_oc(ref_design_bin2(name), scen, seed = 1)
    }
    expect_equal(as.matrix(oc[, c("Go", "NoGo", "Gray", "Miss")]),
                 ref_oc_bin2[[name]], tolerance = 1e-9,
                 ignore_attr = TRUE, info = name)
  }
})

test_that("qdm_oc() simulates studies with two binary endpoints", {
  nsim <- 4000
  expect_warning(oc <- qdm_oc(ref_design_bin2("B1", enumerate = FALSE),
                              scen_bin2, nsim = nsim, seed = 11), "Miss")
  ex <- ref_oc_bin2$B1
  tol <- 4.5 * sqrt(ex * (1 - ex) / nsim) + 0.002
  expect_true(all(abs(as.matrix(oc[, c("Go", "NoGo", "Gray", "Miss")]) -
                        ex) <= tol))
  expect_match(attr(oc, "method"), "^simulation of 4000 studies")
})

test_that("qdm_oc() enumerates the outcomes without a concurrent control", {
  # Recompute the operating characteristics from the decisions for all
  # outcomes of the treatment group
  for (name in c("B7", "B8")) {
    des <- ref_design_bin2(name)
    x <- outcomes_binary2(3)
    dec <- vapply(seq_len(nrow(x)), function(i) {
      qdm_decide(des, data_binary2(x[i, ]))$decision
    }, character(1))
    sc <- scen_bin2[2, c("pi_t1", "pi_t2", "rho_t")]
    w <- dmultinom_counts(x, cellprob_binary2(0.6, 0.5, 0))
    oc <- suppressWarnings(qdm_oc(des, sc))
    for (d in c("Go", "NoGo", "Gray", "Miss")) {
      expect_equal(oc[[d]], sum(w[dec == d]), tolerance = 1e-12,
                   info = paste(name, d))
    }
  }
})

test_that("qdm_oc() reports scenarios that are not evaluable", {
  sc <- rbind(scen_bin2[1, ], data.frame(pi_t1 = 0.35, pi_t2 = 0.15,
                                         rho_t = 0.6, pi_c1 = 0.2,
                                         pi_c2 = 0.15, rho_c = 0.6))
  for (enum in c(TRUE, FALSE)) {
    w <- character(0)
    oc <- withCallingHandlers(
      qdm_oc(ref_design_bin2("B2", enumerate = enum), sc, nsim = 100,
             seed = 1),
      warning = function(e) {
        w <<- c(w, conditionMessage(e))
        invokeRestart("muffleWarning")
      }
    )
    expect_true(any(grepl("Scenario\\(s\\) 2 are not evaluable", w)))
    expect_true(all(is.na(oc[2, c("Go", "NoGo", "Gray", "Miss")])))
    expect_false(anyNA(oc[1, c("Go", "NoGo", "Gray", "Miss")]))
  }
})

# Two continuous endpoints. References: Python_reference/ref_continuous2_qdm.py,
# section D (simulation of 1e5 studies per scenario from the raw data, moment
# matching). The tolerance is 4.5 standard errors of both simulations plus
# 0.002.
ref_oc_cont2 <- list(
  K1 = rbind(
    c(0.02988, 0.48772, 0.4824, 0),
    c(0.22913, 0.32022, 0.45065, 0),
    c(0.50278, 0.30381, 0.1934, 1e-05)
  ),
  K2 = rbind(
    c(0.04448, 0.64277, 0.31275, 0),
    c(0.44666, 0.29693, 0.25641, 0),
    c(0.76227, 0.1782, 0.03907, 0.02046)
  ),
  K3 = rbind(
    c(0.01007, 0.51147, 0.47846, 0),
    c(0.31993, 0.02817, 0.6519, 0),
    c(0.76923, 0.00082, 0.22995, 0)
  ),
  K4 = rbind(
    c(0.03128, 0.86085, 0.10725, 0.00062),
    c(0.55416, 0.14516, 0.29981, 0.00087),
    c(0.9525, 0.00022, 0.04727, 1e-05)
  ),
  K5 = rbind(
    c(0.04574, 0.41433, 0.53993, 0),
    c(0.33466, 0.05605, 0.60929, 0),
    c(0.66123, 0.00532, 0.33345, 0)
  ),
  K6 = rbind(
    c(0, 0.51155, 0.48845, 0),
    c(0.00027, 0.00192, 0.99781, 0),
    c(0.24199, 0, 0.75801, 0)
  )
)

oc_cont2_quiet <- function(...) {
  withCallingHandlers(qdm_oc(...), warning = function(w) {
    if (grepl("Miss", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}

check_oc_cont2 <- function(name, nsim = 2000L) {
  oc <- oc_cont2_quiet(ref_design_cont2(name), scen_cont2(name),
                       nsim = nsim, seed = 2026)
  got <- as.matrix(oc[, c("Go", "NoGo", "Gray", "Miss")])
  ref <- ref_oc_cont2[[name]]
  se <- sqrt(ref * (1 - ref) * (1 / nsim + 1 / 1e5))
  expect_true(all(abs(got - ref) <= 4.5 * se + 0.002), info = name)
  expect_equal(unname(rowSums(got)), rep(1, 3), tolerance = 1e-12)
}

test_that("qdm_oc() matches the references for two continuous endpoints", {
  check_oc_cont2("K1")
  check_oc_cont2("K6")
})

test_that("qdm_oc() matches the references for designs K2 and K3", {
  skip_on_cran()
  check_oc_cont2("K2")
  check_oc_cont2("K3")
})

test_that("qdm_oc() matches the references for designs K4 and K5", {
  skip_on_cran()
  check_oc_cont2("K4")
  check_oc_cont2("K5")
})

test_that("qdm_oc() simulates two continuous endpoints with Monte Carlo", {
  oc <- oc_cont2_quiet(ref_design_cont2("K1", "MC", 2000L),
                       scen_cont2("K1")[2, ], nsim = 400L, seed = 5)
  ref <- ref_oc_cont2$K1[2, ]
  got <- unlist(oc[, c("Go", "NoGo", "Gray", "Miss")])
  se <- sqrt(ref * (1 - ref) / 400)
  expect_true(all(abs(got - ref) <= 4.5 * se + 0.02))
  expect_match(attr(oc, "method"),
               "simulation of 400 studies; .*Monte Carlo simulation")
})

test_that("qdm_oc() matches the exact references with a fixed control", {
  # References: Python_reference/ref_fixed_tte_qdm.py (enumeration of the
  # outcomes)
  oc <- qdm_oc(ref_design("D8"), scen_u)
  expect_identical(attr(oc, "method"), "exact enumeration of the outcomes")
  expect_equal(oc$Go, c(0.00259482740067, 0.244662796684, 0.872478753853,
                        0.999436586302), tolerance = 1e-9)
  expect_equal(oc$NoGo, c(0.411448861957, 0.01596116279, 4.7344970692e-05,
                          7.9776710656e-10), tolerance = 1e-9)
  expect_equal(oc$Miss, rep(0, 4L))
  oc9 <- qdm_oc(ref_design("D9"), data.frame(pi_t = c(0.2, 0.4, 0.6)))
  expect_equal(oc9$Go, c(0.1208738816, 0.6177193984, 0.9452381184),
               tolerance = 1e-9)
  expect_equal(oc9$NoGo, c(0.6777995264, 0.1672897536, 0.0122945536),
               tolerance = 1e-9)
  expect_equal(oc9$Gray, c(0.201326592, 0.214990848, 0.042467328),
               tolerance = 1e-9)
  sc <- data.frame(pi_t1 = c(0.3, 0.5, 0.7), pi_t2 = c(0.4, 0.6, 0.7),
                   rho_t = c(0.2, 0, 0.3))
  ocb <- qdm_oc(ref_design_bin2("BF1"), sc)
  expect_equal(ocb$Go, c(0.0726080130921, 0.216, 0.579202246),
               tolerance = 1e-9)
  expect_equal(ocb$NoGo, c(0.721260675859, 0.338, 0.167444577),
               tolerance = 1e-9)
  expect_equal(ocb$Gray, c(0.206131311049, 0.446, 0.253353177),
               tolerance = 1e-9)
})

test_that("qdm_oc() agrees with the exact references for a fixed continuous control", {
  # References: Python_reference/ref_fixed_tte_qdm.py (noncentral t
  # distributions)
  nsim <- 20000L
  tol <- function(p) 4.5 * sqrt(p * (1 - p) / nsim) + 0.002
  oc7 <- qdm_oc(ref_design_cont("C7"), scen_cont_u, nsim = nsim, seed = 3)
  ref_go <- c(0.000102665519772, 0.4, 0.999325919193)
  ref_nogo <- c(0.63159677607, 0.000947098825199, 2.8762586985e-11)
  expect_true(all(abs(oc7$Go - ref_go) <= tol(ref_go)))
  expect_true(all(abs(oc7$NoGo - ref_nogo) <= tol(ref_nogo)))
  expect_identical(attr(oc7, "method"),
                   paste("simulation of 20000 studies; probabilities by",
                         "exact computation (fixed control)"))
  oc8 <- qdm_oc(ref_design_cont("C8"),
                data.frame(mu_t = c(0, 5, 12), sd_t = c(12, 12, 15)),
                nsim = nsim, seed = 4)
  ref_go <- c(0.0124614836728, 0.26087426553, 0.875577271217)
  ref_nogo <- c(0.903337893177, 0.378241335701, 0.0172420544481)
  expect_true(all(abs(oc8$Go - ref_go) <= tol(ref_go)))
  expect_true(all(abs(oc8$NoGo - ref_nogo) <= tol(ref_nogo)))
})

test_that("qdm_oc() simulates time-to-event designs", {
  # Reference: Python_reference/ref_fixed_tte_qdm.py, design T1
  nsim <- 20000L
  tol <- function(p) 4.5 * sqrt(p * (1 - p) / nsim) + 0.002
  oc <- qdm_oc(ref_design_tte("T1"), data.frame(median_t = c(6.3, 20)),
               nsim = nsim, seed = 2)
  expect_identical(names(oc), c("median_t", "Go", "NoGo", "Gray", "Miss"))
  ref_go <- c(2.48918966507e-09, 0.899733183632)
  ref_nogo <- c(0.205961280815, 2.15285582264e-09)
  expect_true(all(abs(oc$Go - ref_go) <= tol(ref_go)))
  expect_true(all(abs(oc$NoGo - ref_nogo) <= tol(ref_nogo)))
  oc2 <- qdm_oc(ref_design_tte("T2"),
                data.frame(median_t = c(6, 12), median_c = 6), nsim = 200,
                seed = 1)
  expect_identical(names(oc2), c("median_t", "median_c", "Go", "NoGo",
                                 "Gray", "Miss"))
  expect_identical(attr(oc2, "method"),
                   paste("simulation of 200 studies; probabilities by",
                         "numerical integration"))
  expect_equal(oc2$Go + oc2$NoGo + oc2$Gray + oc2$Miss, c(1, 1))
  expect_gt(oc2$Go[2L], oc2$Go[1L])
  expect_error(qdm_oc(ref_design_tte("T2"), data.frame(median_t = 6),
                      nsim = 10),
               "'scenarios' must have the column\\(s\\) median_t, median_c.")
  expect_error(qdm_oc(ref_design_tte("T1"), data.frame(median_t = 0),
                      nsim = 10),
               "'scenarios\\$median_t' must be greater than 0.")
})
