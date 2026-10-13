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
  expect_error(qdm_oc(des, scen_c, nsim = 0), "'nsim' must be at least 1.")
  expect_error(qdm_oc(des, scen_c, seed = 1.5),
               "'seed' must contain whole numbers only.")
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
