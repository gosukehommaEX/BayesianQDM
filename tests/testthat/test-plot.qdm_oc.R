test_that("plot.qdm_oc() draws the decision probabilities", {
  oc <- qdm_oc(ref_design("D1"),
               data.frame(pi_t = c(0.6, 0.2, 0.4), pi_c = 0.4))
  p <- plot(oc)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  # Layers: threshold lines, their labels, curves, points
  expect_equal(sort(unique(b$data[[1]]$xintercept)), c(0.1, 0.2))
  lines <- b$data[[3]]
  expect_equal(nrow(lines), 9L)
  expect_equal(sort(unique(lines$x)), c(-0.2, 0, 0.2))
  # Curves for Go, NoGo and Gray, in increasing order of the effect
  expect_identical(sort(unique(lines$group)), 1:3)
  expect_equal(lines$y[lines$group == 1], oc$Go[order(oc$pi_t)])
  expect_equal(lines$y[lines$group == 3], oc$Gray[order(oc$pi_t)])
  expect_identical(unique(lines$colour[lines$group == 1]), "#658D1B")
})

test_that("plot.qdm_oc() adds Miss decisions and shifted thresholds", {
  expect_warning(oc <- qdm_oc(ref_design("D6"), scen_c), "Miss")
  lines <- ggplot2::ggplot_build(plot(oc))$data[[3]]
  expect_identical(sort(unique(lines$group)), 1:4)
  expect_equal(lines$y[lines$group == 4], oc$Miss)
  expect_identical(unique(lines$colour[lines$group == 4]), "#E69F00")
  oc3 <- qdm_oc(ref_design("D3"), scen_u)
  b3 <- ggplot2::ggplot_build(plot(oc3))
  expect_equal(sort(unique(b3$data[[1]]$xintercept)), 5 / 14 + c(0.1, 0.2))
  expect_equal(sort(unique(b3$data[[3]]$x)), scen_u$pi_t)
})

test_that("plot.qdm_oc() needs distinct values on the horizontal axis", {
  # Two scenarios with the same effect (0.2) but different control rates
  oc <- qdm_oc(ref_design("D1"), data.frame(pi_t = c(0.5, 0.6),
                                            pi_c = c(0.3, 0.4)))
  expect_error(plot(oc), "distinct values on the horizontal axis")
  oc_cont <- qdm_oc(ref_design_cont("C1", "MM"),
                    data.frame(mu_t = 20, sd_t = c(20, 30), mu_c = 0,
                               sd_c = 25), nsim = 20, seed = 1)
  expect_error(plot(oc_cont), "distinct values on the horizontal axis")
})

test_that("plot.qdm_oc() needs the design attribute", {
  oc <- qdm_oc(ref_design("D1"), scen_c)
  attr(oc, "design") <- NULL
  expect_error(plot(oc), "'x' must be an unmodified result of qdm_oc().")
})

test_that("plot.qdm_oc() draws continuous operating characteristics", {
  oc <- qdm_oc(ref_design_cont("C1", "MM"), scen_cont_c, nsim = 200,
               seed = 1)
  b <- ggplot2::ggplot_build(plot(oc))
  expect_equal(sort(unique(b$data[[1]]$xintercept)), c(15, 30))
  expect_equal(sort(unique(b$data[[3]]$x)), c(0, 20, 40))
  oc3 <- qdm_oc(ref_design_cont("C3", "MM"), scen_cont_u, nsim = 200,
                seed = 1)
  b3 <- ggplot2::ggplot_build(plot(oc3))
  expect_equal(sort(unique(b3$data[[1]]$xintercept)), c(5, 20))
  expect_equal(sort(unique(b3$data[[3]]$x)), scen_cont_u$mu_t)
})

test_that("plot.qdm_oc() draws heat maps for two endpoints", {
  sc <- expand.grid(th1 = c(0, 0.2), th2 = c(0, 0.2, 0.4))
  sc <- data.frame(pi_t1 = 0.3 + sc$th1, pi_t2 = 0.3 + sc$th2, rho_t = 0.5,
                   pi_c1 = 0.3, pi_c2 = 0.3, rho_c = 0.5)
  w <- character(0)
  oc <- withCallingHandlers(qdm_oc(ref_design_bin2("B2"), sc),
                            warning = function(e) {
                              w <<- c(w, conditionMessage(e))
                              invokeRestart("muffleWarning")
                            })
  # Scenario 5 (pi_t1 = 0.3, pi_t2 = 0.7, rho_t = 0.5) is infeasible
  expect_identical(w, paste("Scenario(s) 5 are not evaluable: the",
                            "correlation is outside the feasible range for",
                            "the response rates. Their probabilities are",
                            "NA."))
  expect_identical(which(is.na(oc$Go)), 5L)
  p <- plot(oc)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  tiles <- b$data[[1]]
  expect_identical(nrow(tiles), 3L * nrow(sc))
  labels <- b$data[[2]]$label
  expect_identical(sum(labels == "NE"), 3L * sum(is.na(oc$Go)))
  expect_true(all(labels[labels != "NE"] %in%
                    sprintf("%.2f", c(oc$Go, oc$NoGo, oc$Gray))))
  # Several correlations at the same effects cannot be drawn together
  oc2 <- qdm_oc(ref_design_bin2("B2"), rbind(scen_bin2[1, ], scen_bin2[1, ]))
  expect_error(plot(oc2), "with different correlations or control")
  # Without a concurrent control group the axes are the treatment rates
  oc7 <- suppressWarnings(qdm_oc(ref_design_bin2("B7"),
                                 sc[, c("pi_t1", "pi_t2", "rho_t")]))
  n_dec <- 3L + any(oc7$Miss > 0, na.rm = TRUE)
  expect_identical(nrow(ggplot2::ggplot_build(plot(oc7))$data[[1]]),
                   n_dec * nrow(sc))
})

test_that("plot.qdm_oc() draws heat maps for two continuous endpoints", {
  sc <- expand.grid(mu_t1 = c(-10, 10), mu_t2 = c(-1, 1, 3))
  sc$sd_t1 <- 30
  sc$sd_t2 <- 3
  sc$rho_t <- 0.5
  oc <- withCallingHandlers(
    qdm_oc(ref_design_cont2("K4"), sc, nsim = 50L, seed = 1),
    warning = function(w) {
      if (grepl("Miss", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  p <- plot(oc)
  expect_s3_class(p, "ggplot")
  n_dec <- 3L + any(oc$Miss > 0)
  expect_identical(nrow(ggplot2::ggplot_build(p)$data[[1]]),
                   n_dec * nrow(sc))
})

test_that("plot.qdm_oc() draws fixed-control and time-to-event designs", {
  oc <- qdm_oc(ref_design("D8"), scen_u)
  b <- ggplot2::ggplot_build(plot(oc))
  expect_equal(sort(unique(b$data[[1]]$xintercept)), c(0.25, 0.4))
  expect_equal(sort(unique(b$data[[3]]$x)), scen_u$pi_t)
  sc <- data.frame(median_t = c(6.3, 9, 12, 20))
  oc_t <- qdm_oc(ref_design_tte("T1"), sc, nsim = 200, seed = 1)
  b_t <- ggplot2::ggplot_build(plot(oc_t))
  expect_equal(sort(unique(b_t$data[[1]]$xintercept)), c(0.72, 0.87))
  expect_equal(sort(unique(b_t$data[[3]]$x)),
               exp(-log(2) * 3 / sc$median_t))
})
