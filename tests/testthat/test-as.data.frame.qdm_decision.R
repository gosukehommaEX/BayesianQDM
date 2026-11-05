test_that("as.data.frame() gives one row for a decision on one endpoint", {
  des <- ref_design("D1")
  dec <- qdm_decide(des, data_binary(12, 9), data_binary(12, 4))
  df <- as.data.frame(dec)
  expect_s3_class(df, "data.frame")
  expect_identical(names(df), c("decision", "prob_go", "gamma_go", "met_go",
                                "prob_nogo", "gamma_nogo", "met_nogo"))
  expect_identical(nrow(df), 1L)
  expect_identical(df$decision, unname(dec$decision))
  expect_type(df$decision, "character")
  expect_identical(df$prob_go, unname(dec$prob_go))
  expect_identical(df$gamma_go, des$criterion$gamma_go)
  expect_identical(df$met_go, unname(dec$met_go))
  expect_identical(df$prob_nogo, unname(dec$prob_nogo))
  expect_identical(df$gamma_nogo, des$criterion$gamma_nogo)
  expect_identical(df$met_nogo, unname(dec$met_nogo))
})

test_that("as.data.frame() of decisions can be combined with rbind()", {
  des <- ref_design("D1")
  d1 <- qdm_decide(des, data_binary(12, 9), data_binary(12, 4))
  d2 <- qdm_decide(des, data_binary(12, 3), data_binary(12, 6))
  df <- rbind(as.data.frame(d1, row.names = "first"),
              as.data.frame(d2, row.names = "second"))
  expect_identical(rownames(df), c("first", "second"))
  expect_identical(df$decision, unname(c(d1$decision, d2$decision)))
  expect_identical(df$prob_go, unname(c(d1$prob_go, d2$prob_go)))
  expect_identical(df$prob_nogo, unname(c(d1$prob_nogo, d2$prob_nogo)))
  # A changed probability is not hidden by the conversion
  d3 <- d2
  d3$prob_go <- d3$prob_go + 1e-3
  expect_false(identical(as.data.frame(d3)$prob_go,
                         as.data.frame(d2)$prob_go))
})

test_that("as.data.frame() adds the region probabilities of two endpoints", {
  dec <- qdm_decide(ref_design_cont2("K1"), ref_data_cont2$K1$t,
                    ref_data_cont2$K1$c)
  df <- as.data.frame(dec)
  expect_identical(names(df),
                   c("decision", "prob_go", "gamma_go", "met_go", "prob_nogo",
                     "gamma_nogo", "met_nogo", paste0("R", 1:9)))
  expect_identical(unlist(df[paste0("R", 1:9)]), dec$prob_regions)
  expect_identical(df$prob_go, unname(dec$prob_go))
})

test_that("as.data.frame() has one column per endpoint for marginal rules", {
  dec <- qdm_decide(ref_design_cont2("K5"), ref_data_cont2$K5$t)
  df <- as.data.frame(dec)
  expect_identical(names(df),
                   c("decision", "prob_go1", "prob_go2", "gamma_go1",
                     "gamma_go2", "met_go", "prob_nogo1", "prob_nogo2",
                     "gamma_nogo1", "gamma_nogo2", "met_nogo",
                     paste0("R", 1:9)))
  expect_identical(c(df$prob_go1, df$prob_go2), unname(dec$prob_go))
  expect_identical(c(df$prob_nogo1, df$prob_nogo2), unname(dec$prob_nogo))
  expect_identical(c(df$gamma_go1, df$gamma_go2), c(0.6, 0.5))
  expect_identical(c(df$gamma_nogo1, df$gamma_nogo2), c(0.5, 0.6))
  expect_identical(df$met_go, unname(dec$met_go))
})

test_that("as.data.frame() adds the posterior parameters on request", {
  dec <- qdm_decide(ref_design("D1"), data_binary(12, 9), data_binary(12, 4))
  df <- as.data.frame(dec, posterior = TRUE)
  expect_identical(names(df)[8:11],
                   c("post_t_a", "post_t_b", "post_c_a", "post_c_b"))
  expect_identical(c(df$post_t_a, df$post_t_b), c(10, 4))
  expect_identical(c(df$post_c_a, df$post_c_b), c(5, 9))
  fixed <- qdm_decide(ref_design("D8"), data_binary(20, 9))
  expect_identical(as.data.frame(fixed, posterior = TRUE)$post_c_value, 0.2)
  expect_identical(ncol(as.data.frame(dec)), 7L)
  expect_error(as.data.frame(dec, posterior = NA),
               "'posterior' must be TRUE or FALSE.", fixed = TRUE)
})
