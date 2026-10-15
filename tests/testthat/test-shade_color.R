# Expected values from Python_reference/ref_plot_colors.py: linear
# interpolation in RGB between white (probability 0) and the color
# (probability 1), rounded to integers.
test_that("shade_color() interpolates between white and the color", {
  p <- c(0, 0.2, 0.55, 0.8, 1)
  expect_identical(shade_color(p, "#658D1B"),
                   c("#FFFFFF", "#E0E8D1", "#AAC082", "#84A449", "#658D1B"))
  expect_identical(shade_color(p, "#D91E49"),
                   c("#FFFFFF", "#F7D2DB", "#EA839B", "#E14B6D", "#D91E49"))
  expect_identical(shade_color(p, "#939597"),
                   c("#FFFFFF", "#E9EAEA", "#C4C5C6", "#A9AAAC", "#939597"))
  expect_identical(shade_color(p, "#E69F00"),
                   c("#FFFFFF", "#FAECCC", "#F1CA73", "#EBB233", "#E69F00"))
  expect_identical(shade_color(p, "black"),
                   c("#FFFFFF", "#CCCCCC", "#737373", "#333333", "#000000"))
  expect_identical(shade_color(p, "darkgreen"),
                   c("#FFFFFF", "#CCE0CC", "#73AA73", "#338333", "#006400"))
})

test_that("shade_color() takes one color per probability", {
  expect_identical(
    shade_color(c(0.2, NA, 0.8, 1),
                c("#658D1B", "#D91E49", "black", "darkgreen")),
    c("#E0E8D1", "#FFFFFF", "#333333", "#006400")
  )
  # Probabilities outside [0, 1] are clipped (without clipping, -0.2 and
  # 1.3 would give the channel values 306 and -77)
  expect_identical(shade_color(c(-0.2, 1.3), "black"),
                   c("#FFFFFF", "#000000"))
})
