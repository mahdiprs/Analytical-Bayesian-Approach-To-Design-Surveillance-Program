#########################################################
#
# R Examples for the survey block model
#  using numerical method to calculate statistics of interest
#
#########################################################
library(calculus)
library(ggplot2)

# import functions
source("./R/funcsNumerical.R")

# set the model parameters
p0 <- 0.3 # Allow different p0
p <- 0.5 # probability of pest presence
rho0<- -log(1 - p0) # assumed prior mean for X0, given p0
rho <- -log(1 - p) # assumed prior mean for X0, given p
delta <- 0.1 # survey sensitivity level
lambda <- 1.35 # population growth rate >1
tr <- 0.95 # target probability of pest absence

# estimate the block size kappa with target probability of absence tr (T)
k0 <- kappa(rho0, lambda, delta, tr)
k0

k <- kappa(rho, lambda, delta, tr)
k

# Expected number of mop-ups given block size k
ExM(k0,k, rho0, rho, lambda, delta)

# variance number of mop-ups given block size k
VarM(k0,k, rho0, rho, lambda, delta)

# Expected number of surveys given block size k
ExN(k0,k, rho0, rho, lambda, delta)

# variance number of surveys given block size k
VarN(k0,k, rho0, rho, lambda, delta)

