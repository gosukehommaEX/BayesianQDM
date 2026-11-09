test_that("report_n_search() describes the search for the sample size", {
  res <- find_n_d1(n_t = seq(6, 20, by = 2))
  blocks <- report_n_search(res, 3L)
  expect_identical(vapply(blocks, function(b) b$type, character(1)),
                   c("heading", "par", "table", "par", "table", "par"))
  expect_identical(blocks[[1]]$text, "Sample size")
  text <- blocks[[2]]$text
  expect_match(text, paste("The sample sizes $n_t = 14$ and $n_c = 14$ were",
                           "chosen from 8 candidates, with $n_t$ from 6 to 20",
                           "(exact enumeration of the outcomes)."),
               fixed = TRUE)
  expect_match(text, paste("They are those of the smallest candidate that",
                           "meets two targets. Under the scenario",
                           "$\\pi_t = 0.4$, $\\pi_c = 0.4$, under which a Go",
                           "decision would be wrong, the probability that the",
                           "Go criterion is met is at most 0.1 (it is 0.088);",
                           "under the scenario $\\pi_t = 0.7$, $\\pi_c = 0.4$,",
                           "under which a NoGo decision would be wrong, the",
                           "probability that the NoGo criterion is met is at",
                           "most 0.1 (it is 0.033)."),
               fixed = TRUE)
  expect_match(text, paste("all candidates from $n_t = 20$ and $n_c = 20$ on",
                           "meet the targets."), fixed = TRUE)
  expect_identical(blocks[[3]]$header,
                   c("$n_t$", "$n_c$",
                     "P(Go criterion met when Go would be wrong)",
                     "P(NoGo criterion met when NoGo would be wrong)",
                     "Targets met"))
  expect_identical(blocks[[3]]$rows[5L, ],
                   c("14", "14", "0.088", "0.033", "yes"))
  expect_identical(blocks[[3]]$rows[1L, 5L], "no")
  expect_identical(blocks[[5]]$header[1L], "Scenario")
  expect_identical(blocks[[5]]$rows[, 1L],
                   c("Go would be wrong", "NoGo would be wrong"))
  expect_identical(blocks[[6]]$text, report_legend(res$design))
})

test_that("report_n_search() describes the other outcomes of the search", {
  res <- find_n_d1(n_t = c(14, 16))
  expect_match(report_n_search(res, 3L)[[2]]$text,
               "The largest candidate does not meet the targets.",
               fixed = TRUE)
  res2 <- find_n_d1(n_t = c(20, 24))
  expect_false(grepl("larger|largest|monotone",
                     report_n_search(res2, 3L)[[2]]$text))
})
