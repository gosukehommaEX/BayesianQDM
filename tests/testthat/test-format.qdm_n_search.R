test_that("format.qdm_n_search() describes the chosen sample sizes", {
  res <- find_n_d1(n_t = seq(6, 20, by = 2))
  head <- c("Search of the sample size (exact enumeration of the outcomes)",
            paste("  Targets: P(Go criterion met | scenario_go) <= 0.1 and",
                  "P(NoGo criterion met | scenario_nogo) <= 0.1"))
  expect_identical(format(res), c(
    head,
    "  Smallest candidate that meets the targets: n_t = 14, n_c = 14",
    "    P(Go criterion met | scenario_go) = 0.08785",
    "    P(NoGo criterion met | scenario_nogo) = 0.03288",
    paste("  Smallest candidate from which all larger candidates meet the",
          "targets: n_t = 20, n_c = 20")
  ))
  expect_identical(format(res, digits = 2L)[4L],
                   "    P(Go criterion met | scenario_go) = 0.088")
})

test_that("format.qdm_n_search() covers the other outcomes of the search", {
  # Every larger candidate meets the targets: no further line
  expect_length(format(find_n_d1(n_t = c(20, 24))), 5L)
  # The largest candidate does not meet the targets
  expect_identical(format(find_n_d1(n_t = c(14, 16)))[6L],
                   paste("  The largest candidate does not meet the targets",
                         "(the probabilities are not monotone in the sample",
                         "size)."))
  none <- suppressWarnings(find_n_d1(n_t = c(6, 8)))
  expect_identical(format(none)[3L],
                   paste("  No candidate from n_t = 6, n_c = 6 to n_t = 8,",
                         "n_c = 8 meets the targets."))
})
