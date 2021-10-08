#########################################################
#
# R Examples for the survey block model
#  using symbolic algebra to calculate statistics of interest
#
#########################################################

library(symengine)
library(ggplot2)
library(magrittr)
library(gridExtra)

# import functions

source("R/funcsSymbolic.R")

# set the model parameters
p0 <- 0.3 # Allow different p0
p <- 0.5 # probability of pest presence
rho0<- -log(1 - p0) # assumed prior mean for X0, given p0
rho <- -log(1 - p) # assumed prior mean for X0, given p
delta <- 0.1 # survey sensitivity level
lambda <- 1.35 # population growth rate >1
tr <- 0.95 # target probability of pest absence


# estimate the block size kappa with target probability of absence tr (T)
# not necessary to do this symbolically, but for clarity.

# First survey block
k0 <- kappa(rho0, lambda, delta, tr, p0=TRUE)
k0

# subsequent survey blocks
k <- kappa(rho, lambda, delta, tr, p0=FALSE)
k

# Parms now has rho and rho0
parms<-list(rho=rho, rho0=rho0, lambda=lambda, delta=delta)

#--------------
# Expected number (and variance) of mop-ups given block size k

# p0 = p
pgf<- make_pgfM(k)
calc_moments(pgf, parms)

# Calculate the PMF for p0 = p
pmfM<- calc_pmf(pgf, parms, support=0:15)

# now for case p0 != p
pgf<- make_pgfM(k, k0)
calc_moments(pgf, parms)

# Calculate the PMF for p0 != p
pmfM0<- calc_pmf(pgf, parms, support=0:15)

#... and save plot
p1<- pmfM %>% ggplot(aes(as.factor(n), p, group = 1)) +
      geom_path() +
      geom_point() +
      geom_path(aes(as.factor(n), p, group = 1), linetype=2, data=pmfM0) +
      geom_point(aes(as.factor(n), p, group = 1), shape=1, data=pmfM0) +
      geom_path(aes(as.factor(X), prX,group = 1), linetype=1, colour ="red",data=subset(dataPaperPMF, val =="M" & col=="eq")) +
      geom_path(aes(as.factor(X), prX,group = 1), linetype=2, colour ="red",data=subset(dataPaperPMF, val =="M" & col=="nEq")) +
      labs(x = "Mop-ups (M)", y = "PMF") +
      scale_y_continuous(breaks = seq(0,1,0.2), limits=c(0,1)) +
      theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black")
    )
p1
#--------------
# Survey number of first detection given block size k

# p0 = p
pgf<- make_pgfK(k)
calc_moments(pgf, parms)

# Calculate the PMF
pmfK<- calc_pmf(pgf, parms, support=0:15)

# now for case p0 != p
pgf<- make_pgfK(k, k0)
calc_moments(pgf, parms)

# Calculate the PMF
pmfK0<- calc_pmf(pgf, parms, support=0:15)

#... and save plot
p2<- pmfK %>% ggplot(aes(as.factor(n), p, group = 1)) +
      geom_path() +
      geom_point() +
      geom_path(aes(as.factor(n), p, group = 1), linetype=2, data=pmfK0) +
      geom_point(aes(as.factor(n), p, group = 1), shape=1, data=pmfK0) +
      geom_path(aes(as.factor(X), prX,group = 1), linetype=2, colour ="red",data=subset(dataPaperPMF, val =="Y" & col=="pnote")) +
      geom_path(aes(as.factor(X), prX,group = 1), linetype=1, colour ="red",data=subset(dataPaperPMF, val =="Y" & col=="p")) +
      labs(x = "Survey-number of first detection (Y)", y = "PMF") +
      scale_y_continuous(breaks = seq(0,0.5,0.1), limits=c(0, 0.5)) +
  scale_x_discrete(limits=0:16, labels = 0:16)+
      theme(
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")
    )
p2
#--------------
# Expected number of surveys (and variance) given block size k

# case p0 = p
pgf<- make_pgfN(k)
calc_moments(pgf, parms)

# Calculate the PMF
pmfN<- calc_pmf(pgf, parms, support=0:15)

# now for case p0 != p
pgf<- make_pgfN(k, k0)
calc_moments(pgf, parms)

# Calculate the PMF
pmfN0<- calc_pmf(pgf, parms, support=0:15)

#... and plot
p3<- pmfN %>% ggplot(aes(as.factor(n), p, group = 1)) +
      geom_path() +
      geom_point() +
      geom_path(aes(as.factor(n), p, group = 1), linetype=2, data=pmfN0) +
      geom_point(aes(as.factor(n), p, group = 1), shape=1, data=pmfN0) +
      geom_path(aes(as.factor(X), prX,group = 1), linetype=2, colour ="red",data=subset(dataPaperPMF, val =="N" & col=="eq")) +
      geom_path(aes(as.factor(X), prX,group = 1), linetype=1, colour ="red",data=subset(dataPaperPMF, val =="N" & col=="nEq")) +
      labs(x = "number of surveys (N)", y = "PMF") +
      scale_y_continuous(breaks = seq(0,1,0.2), limits=c(0, 1)) +
      scale_x_discrete(limits=0:16, labels = 0:16)+
      theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black")
    )

p3
#win.graph(12, 5)
png("figure2.png",width=12,height=5,units="in",res=300)
grid.arrange(p1, p2, p3, nrow=1)
dev.off()



