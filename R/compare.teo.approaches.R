library(dplyr)
library(tidyr)
library(ggplot2)
# Set scenarios
scenario = tibble(
  prior = c('N-Inv-Chisq', 'vague')
) %>%
  group_by_all() %>%
  reframe(
    mu1 = seq(2.5, 5, by = 0.1),
    mu2 = 0
  )
# Result using convolution (i.e., using integrate)
result.convolution = scenario %>%
  group_by_all() %>%
  reframe(
    Approach = 'Convolution',
    BayesDecisionProbContinuous(
      nsim = 10000, design = 'controlled', prob = 'posterior', prior = prior, approx = FALSE, theta0 = c(2, 0), gamma1 = 0.8, gamma2 = 0.3,
      n1 = 12, n2 = 12, m1 = NULL, m2 = NULL, kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5, mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5),
      mu1 = mu1, mu2 = mu2, sigma1 = 1, sigma2 = 1, r = NULL, seed = 1
    )
  )
# Result using Welch-Satterthwaite approximation
result.WSapprox = scenario %>%
  group_by_all() %>%
  reframe(
    Approach = 'WS.approx',
    BayesDecisionProbContinuous(
      nsim = 10000, design = 'controlled', prob = 'posterior', prior = prior, approx = TRUE, theta0 = c(2, 0), gamma1 = 0.8, gamma2 = 0.3,
      n1 = 12, n2 = 12, m1 = NULL, m2 = NULL, kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5, mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5),
      mu1 = mu1, mu2 = mu2, sigma1 = 1, sigma2 = 1, r = NULL, seed = 1
    )
  )
# Display figure for comparing results by two different approaches
result.convolution %>%
  bind_rows(result.WSapprox) %>%
  pivot_longer(
    cols = c(Go, NoGo, Gray), names_to = 'Decision', values_to = 'Prob'
  ) %>%
  mutate(
    theta = mu1 - mu2,
    prior = factor(prior, levels = c('N-Inv-Chisq', 'vague')),
    Approach = factor(Approach, levels = c('Convolution', 'WS.approx')),
    Decision = factor(Decision, levels = c('Go', 'Gray', 'NoGo'))
  ) %>%
  ggplot(aes(x = theta, y = Prob)) +
  geom_line(aes(colour = Decision, linetype = Approach), linewidth = 1) +
  theme_bw() +
  facet_grid(. ~ prior) +
  scale_color_manual(
    values = c('Go' = '#658D1B', 'Gray' = '#939597', 'NoGo' = '#D91E49'),
    labels =  c('Go', 'Gray', 'NoGo')
  ) +
  scale_x_continuous(
    #expand=c(0, 0),
    limits = c(2.5, 5),
    breaks = seq(2.5, 5, l = 6)
  ) +
  scale_y_continuous(
    #expand=c(0, 0),
    limits = c(0, 1),
    breaks = seq(0, 1, l = 11)
  ) +
  labs(
    title = 'Probability of making Go/Gray/NoGo decision',
    x = expression(theta),
    y = 'Probability'
  ) +
  theme(
    text = element_text(size = 10),
    legend.background = element_rect(fill = 'white', color = 'black'),
    legend.key.width = unit(2, 'cm'),
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = 'bottom'
  )


result.convolution %>%
  bind_rows(result.WSapprox) %>%
  group_by(prior, mu1) %>%
  mutate(
    Go.diff = abs(Go[Approach == 'Convolution'] - Go[Approach == 'WS.approx']),
    #NoGo.diff = abs(NoGo[Approach == 'Convolution'] - NoGo[Approach == 'WS.approx']),
    #Gray.diff = abs(Gray[Approach == 'Convolution'] - Gray[Approach == 'WS.approx']),
  ) %>%
  group_by(prior) %>%
  filter(
    Go.diff == max(Go.diff)
  ) %>%
  arrange(
    prior
  )
