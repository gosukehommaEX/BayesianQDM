test_that("format.qdm_boundary() tabulates the intervals by condition", {
  b <- qdm_boundary(ref_design("D8"), n_t = c(10, 20))
  f <- format(b)
  expect_identical(f[1:3], c("Decision boundary of a design with one binary endpoint",
                             "", "Responders in the treatment group (y_t)"))
  expect_identical(f[4:6], c("  n_t   NoGo   Gray   Go",
                             "  10    0-1    2-5    6-10",
                             "  20    0-3    4-9    10-20"))
  expect_true("Observed response rate (y_t / n_t)" %in% f)
  expect_true(any(grepl("^  10 +0\\.0% to 10\\.0% +20\\.0% to 50\\.0% +60\\.0% to 100\\.0%$",
                        f)))
  expect_true(any(grepl("fixed control rate \\(y_t / n_t - 0.2\\)", f)))
  expect_true(any(grepl("^  20 +-0\\.2 to -0\\.05 +0 to 0\\.25 +0\\.3 to 0\\.8$",
                        f)))
  # At most n conditions
  f1 <- format(b, n = 1)
  expect_identical(f1[length(f1)], "... 1 more value of n_t (see the data frame)")
  f3 <- format(qdm_boundary(ref_design("D8"), n_t = c(10, 20, 30)), n = 1)
  expect_identical(f3[length(f3)],
                   "... 2 more values of n_t (see the data frame)")
  expect_error(format(b, n = 0), "'n' must be at least 1.")
  expect_error(format(b, digits = 0.5),
               "'digits' must contain whole numbers only.")
})

test_that("format.qdm_boundary() shows open intervals of observed means", {
  b <- qdm_boundary(ref_design_cont("C1"), sd_t = 25,
                    data_c = ref_data_cont$C1$c)
  f <- format(b, digits = 4)
  expect_identical(f[2:3], c("  sd_t = 25",
                             "  Control data: continuous, n = 12, mean = -12, sd = 25"))
  expect_true(any(grepl("^  12 +<= -2\\.686 +-2\\.686 to 20\\.74 +>= 20\\.74$",
                        f)))
  expect_true(any(grepl("^  n_t +NoGo +Gray +Go$", f)))
  # Intervals unbounded on both sides
  b$from[1L] <- -Inf
  b$to[1L] <- Inf
  expect_true(any(grepl("^  12 +all", format(b))))
})

test_that("format.qdm_boundary() orders the decisions by the observed quantity", {
  # More events are unfavorable: Go comes first
  b <- qdm_boundary(ref_design_hr("H1"), exposure_t = 200,
                    data_c = ref_data_hr$c)
  f <- format(b)
  expect_true(any(grepl("^  exposure_t +Go +Gray +NoGo$", f)))
  expect_true(any(grepl("^  200 +0-14 +15-25 +26-30$", f)))
  expect_true(any(grepl("hazard ratio", f)))
  bt <- qdm_boundary(ref_design_tte("T1"), n_t = 10, exposure_t = 100)
  ft <- format(bt)
  expect_true(any(grepl("^Observed milestone survival", ft)))
  expect_true(any(grepl("^  100 +0-4 +5-10$", ft)))
  # A decision that does not occur for a condition
  bd <- qdm_boundary(ref_design("D1"))
  expect_true(any(grepl("^  12 +0-12 +- +-$", format(bd, n = 13))))
})

test_that("format.qdm_boundary() formats a subset as a data frame", {
  b <- qdm_boundary(ref_design("D8"))
  sub <- b[1:2, ]
  attr(sub, "axes") <- NULL
  f <- format(sub)
  expect_s3_class(f, "data.frame")
  expect_identical(as.character(f$decision), c("NoGo", "Gray"))
  # format.data.frame() returns character columns of class "AsIs"
  expect_identical(as.character(f$to), c("3", "9"))
})
