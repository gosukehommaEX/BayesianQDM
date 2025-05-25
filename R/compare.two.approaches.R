# library(dplyr)
# library(tidyr)
# library(purrr)
# library(ggplot2)
# library(ggh4x)
# library(scales)
# library(patchwork)
# results = tibble(
#   Approach = c('Convolution', 'WS.approx'),
#   approx = c(FALSE, TRUE)
# ) %>%
#   group_by_all() %>%
#   reframe(
#     design = rep(c('controlled', 'uncontrolled'), each = 4),
#     prob = rep(rep(c('posterior', 'predictive'), each = 2), 2),
#     prior = rep(c('N-Inv-Chisq', 'vague'), 4)
#   ) %>%
#   group_by_all() %>%
#   mutate(
#     muk = list(
#       tibble(
#         mu1 = seq(0, 8, by = 0.2),
#         mu2 = 0
#       )
#     )
#   ) %>%
#   reframe(
#     GoNoGoGray = map(list(1), ~ {
#       muk[[1]] %>%
#         group_by_all() %>%
#         reframe(
#           BayesDecisionProbContinuous(
#             nsim = 10000, prob = prob, design = design, prior = prior, approx = approx,
#             theta.TV = 2, theta.MAV = 0, theta.NULL = 0.5, gamma1 = 0.8, gamma2 = 0.3,
#             n1 = 12, n2 = 12, m1 = 120, m2 = 120, kappa01 = 5, kappa02 = 5,
#             nu01 = 5, nu02 = 5, mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5),
#             mu1 = mu1, mu2 = mu2, sigma1 = 1, sigma2 = 1, r = 12, seed = 1
#           )
#         )
#     })
#   ) %>%
#   unnest(GoNoGoGray) %>%
#   mutate(
#     theta = mu1 - mu2
#   )
# Display figure for comparing results by two different approaches
# results %>%
#   pivot_longer(
#     cols = c(Go, NoGo, Gray), names_to = 'Decision', values_to = 'Prob'
#   ) %>%
#   mutate(
#     Approach = factor(Approach, levels = c('Convolution', 'WS.approx')),
#     design = factor(design, levels = c('controlled', 'uncontrolled')),
#     prob = factor(prob, levels = c('posterior', 'predictive')),
#     prior = factor(prior, levels = c('N-Inv-Chisq', 'vague')),
#     Decision = factor(Decision, levels = c('Go', 'Gray', 'NoGo'))
#   ) %>%
#   ggplot(aes(x = theta, y = Prob)) +
#   facet_nested(
#     prob ~ design + prior,
#     nest_line = element_line(colour = 'black')
#   ) +
#   geom_line(aes(colour = Decision, linetype = Approach), linewidth = 1) +
#   theme_bw() +
#   scale_color_manual(
#     values = c('Go' = '#658D1B', 'Gray' = '#939597', 'NoGo' = '#D91E49'),
#     labels =  c('Go', 'Gray', 'NoGo')
#   ) +
#   scale_x_continuous(
#     #expand=c(0, 0),
#     limits = c(0, 8),
#     breaks = seq(0, 8, l = 9)
#   ) +
#   scale_y_continuous(
#     #expand=c(0, 0),
#     limits = c(0, 1),
#     breaks = seq(0, 1, l = 11)
#   ) +
#   labs(
#     title = 'Probability of making Go/Gray/NoGo decision',
#     x = expression(theta),
#     y = 'Probability'
#   ) +
#   theme(
#     text = element_text(size = 10),
#     legend.background = element_rect(fill = 'white', color = 'black'),
#     legend.key.width = unit(2, 'cm'),
#     legend.text = element_text(size = 10),
#     legend.title = element_blank(),
#     legend.position = 'bottom'
#   )
