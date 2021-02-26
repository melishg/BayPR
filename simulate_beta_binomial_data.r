# Simulate data for beta-binomial models
# Author: Josh Betz (https://github.com/jbetz-jhu)

library(optimx)
library(dplyr)

set.seed(5924850) # Set RNG seed to reproduce a result

# Path to `source-bayes_prevalence_ratio.r` from BayPR Github
source(file.path("path", "to", "source-bayes_prevalence_ratio.r"))

# Path to output - Dataset 1: fixed r_{k}, one component
output.path.1 <-
  file.path("set", "path", "for", "single_bb_example.csv")

# Path to output - Dataset 2: fixed r_{k}, two component mixture
output.path.2 <-
  file.path("set", "path", "for", "mixture_bb_example.csv")

### Simulation Dataset 1 - One Component, Fixed r_{k} ##########################
# Number of variants to simulate
n.sim.data <- 200
# Control variant sample size
n.0 <- 250000
# Relative sample size in cases (n.1) relative to controls (n.0)
r.k <- 0.002
# Case sample size
n.1 <- r.k*n.0

# Distribution of variant prevalences in controls: Beta distribution
mean.pi.0 <- 0.05
var.pi.0 <- 0.000025
m.pi.0 <- mean.pi.0*(1 - mean.pi.0)/var.pi.0 - 1

# Mean/Variance of Prevalence Ratios: (gamma_{k}) cases/controls
pr.mean.1 <- 1
pr.var.1 <- 0.025

# Convert mean/variance of control variant prevalence to alpha/beta
alpha.pi.0 <- mean.pi.0*m.pi.0
beta.pi.0 <- (1 - mean.pi.0)*m.pi.0

# Convert prevalence ratio means/variances to alpha/beta
theta.k.mean.1 <- pr.mean.1*r.k
theta.k.var.1 <- pr.var.1*r.k^2

beta.1 <- 2 + (theta.k.mean.1*(1 + theta.k.mean.1))/theta.k.var.1
alpha.1 <- theta.k.mean.1*(beta.1 - 1)


# Mean value of gamma

# Simulate event data
sim.data <-
  data.frame(event = 1:n.sim.data,
             n.0 = n.0,
             n.1 = n.1,
             # Transformed Prevalence Ratio
             theta = rbeta(n = n.sim.data,
                           shape1 = alpha.1,
                           shape2 = beta.1)) %>%
  dplyr::mutate(gamma = theta/((1-theta)*r.k), # Prevalence Ratio
                pi.0 = rbeta(n = n.sim.data,
                             shape1 = alpha.pi.0,
                             shape2 = beta.pi.0),
                pi.1 = gamma*pi.0,
                y.0 = rbinom(n = n.sim.data, size = n.0, prob = pi.0),
                y.1 = rbinom(n = n.sim.data, size = n.1, prob = pi.1),
                t.k = y.1 + y.0,
                r.k = n.1/n.0,
                theta.k.hat = y.1/t.k,
                pi.0.hat = y.0/n.0,
                pi.1.hat = y.1/n.1,
                gamma.hat = pi.1.hat/pi.0.hat,
                flag = "include")

if(sum(sim.data$pi.1 > 1) > 0)
  stop("Outcome probability in case population > 1:\n Check null distribution ",
       "and prevalence ratio distribution.")
if(sum(sim.data$t.k < 1)) {
  warning("At least one observation with no observed events. These will be ",
          "removed from the dataset.")
  sim.data <-
    sim.data %>%
    dplyr::filter(sim.data$t.k > 0)
}


### Simulation Dataset 2 - Two Component, Fixed r_{k} ##########################
# Population 1 - Mean/Variance of Prevalence Ratios: (gamma_{k}) cases/controls
pr.mean.2 <- 3
pr.var.2 <- 1.5

# Proportion of variants from population 2
pr.component.2 <- 0.25

# Convert prevalence ratio means/variances to alpha/beta
theta.k.mean.2 <- pr.mean.2*r.k
theta.k.var.2 <- pr.var.2*r.k^2

beta.2 <- 2 + (theta.k.mean.2*(1 + theta.k.mean.2))/theta.k.var.2
alpha.2 <- theta.k.mean.2*(beta.2 - 1)

# Simulate event data
sim.data.2 <-
  data.frame(event = 1:n.sim.data,
             n.0 = n.0,
             n.1 = n.1,
             # Transformed Prevalence Ratio
             theta.1 = rbeta(n = n.sim.data,
                             shape1 = alpha.1,
                             shape2 = beta.1),
             theta.2 = rbeta(n = n.sim.data,
                             shape1 = alpha.2,
                             shape2 = beta.2),
             is.component.2 = rbinom(n = n.sim.data,
                                     size = 1,
                                     prob =  pr.component.2)) %>%
  dplyr::mutate(theta = (1 - is.component.2)*theta.1 + is.component.2*theta.2,
                gamma = theta/((1 - theta)*r.k), # Prevalence Ratio
                pi.0 = rbeta(n = n.sim.data,
                             shape1 = alpha.pi.0,
                             shape2 = beta.pi.0),
                pi.1 = gamma*pi.0,
                y.0 = rbinom(n = n.sim.data, size = n.0, prob = pi.0),
                y.1 = rbinom(n = n.sim.data, size = n.1, prob = pi.1),
                t.k = y.1 + y.0,
                r.k = n.1/n.0,
                theta.k.hat = y.1/t.k,
                pi.0.hat = y.0/n.0,
                pi.1.hat = y.1/n.1,
                gamma.hat = pi.1.hat/pi.0.hat,
                flag = "include")

if(sum(sim.data.2$pi.1 > 1) > 0)
  stop("Outcome probability in case population > 1:\n Check null distribution ",
       "and prevalence ratio distribution.")
if(sum(sim.data.2$t.k < 1) > 0) {
  warning("At least one observation with no observed events. These will be ",
          "removed from the dataset.")
  sim.data.2 <-
    sim.data.2 %>%
    dplyr::filter(sim.data.2$t.k > 0)
}


write.csv(x = sim.data,
          file = output.path.1,
          row.names = FALSE)

write.csv(x = sim.data.2,
          file = output.path.2,
          row.names = FALSE)