test_that("plot_compare_two() draws heat maps for two endpoints", {
  d1 <- ref_design_cont2("K1")
  d2 <- qdm_sample_size(d1, 20, 20)
  scen <- scen_cont2("K1")
  cmp <- qdm_compare(list(a = d1, b = d2), scen, nsim = 200L, seed = 1)
  p <- plot(cmp)
  b <- ggplot2::ggplot_build(p)
  n_dec <- if (any(cmp$Miss > 0)) 4L else 3L
  # One row of panels per design and one column per decision
  expect_identical(nrow(b$layout$layout), 2L * n_dec)
  tiles <- b$data[[1]]
  expect_identical(nrow(tiles), 2L * n_dec * nrow(scen))
  text <- b$data[[2]]
  expect_identical(text$label[text$PANEL == 1],
                   sprintf("%.2f", cmp$Go[cmp$design == "a"]))
  expect_identical(text$label[text$PANEL == n_dec + 1L],
                   sprintf("%.2f", cmp$Go[cmp$design == "b"]))
  expect_identical(tiles$fill[tiles$PANEL == 1],
                   shade_color(cmp$Go[cmp$design == "a"], "#658D1B"))
  p2 <- plot(cmp, decisions = "NoGo", cell_values = FALSE, col_nogo = "blue")
  b2 <- ggplot2::ggplot_build(p2)
  expect_identical(unique(b2$data[[2]]$label), "")
  expect_identical(b2$data[[1]]$fill[b2$data[[1]]$PANEL == 1],
                   shade_color(cmp$NoGo[cmp$design == "a"], "blue"))
})

test_that("plot_compare_two() needs the same axes", {
  scen <- scen_cont2("K1")
  expect_warning(
    cmp <- qdm_compare(list(a = ref_design_cont2("K1"),
                            b = ref_design_cont2("K4")),
                       scen, nsim = 100L, seed = 1),
    "b: The Go and NoGo criteria are both met (Miss)", fixed = TRUE
  )
  expect_error(plot(cmp), "The designs do not have the same axes",
               fixed = TRUE)
})
