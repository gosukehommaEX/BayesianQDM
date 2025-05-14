
# BayesianQDM

## Overview


`BayesianQDM` is an R package that provides methods to calculate posterior probabilities and posterior predictive probabilities,
    and Go, NoGo and Gray probabilities for quantitative decision-making framework under bayesian paradigm.

For technical details about the methodology, please refer to Kang et al.(20XX).


## Installation


You can install the development version of BayesianQDM from GitHub with:


```r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/BayesianQDM")
```


## Usage


```r
library(BayesianQDM)


# Calculate bayesian posterior probability for binary endpoints
BayesPostPredBinary(
  prob = 'posterior', external = TRUE, theta0 = 0.15,
  n1 = 12, n2 = 15, y1 = 7, y2 = 9, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  m1 = NULL, m2 = NULL, ne1 = 12, ne2 = 12, ye1 = 6, ye2 = 6, ae1 = 0.5, ae2 = 0.5
)

# Calculate Go, NoGo and Gray probabilities using posterior predictive probability for controlled trial with binary endpoints 
BayesDecisionProbBinary(
  prob = 'predictive', design = 'controlled', theta.TV = NULL, theta.MAV = NULL, theta.NULL = 0, gamma1 = 0.9, gamma2 = 0.3,
  pi1 = c(0.2, 0.4, 0.6, 0.8), pi2 = rep(0.2, 4), n1 = 12, n2 = 12, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5, z = NULL,
  m1 = 30, m2 = 30, ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
)
```


## References


Kang et al.(20XX)). Title
