# kappa with  T= tr
kappa <- function(rho, lambda, delta,tr) {
  t <- 1
  while ((psi_Xk(t, rho, lambda, delta, 0) < tr)) {
    t <- t + 1
  }
  return(t)
}
# PGFM
pgfM <- function(t0,t,rho0,rho,lambda,delta,s){
  t1 <- PND(t0, rho0, lambda, delta)
  t2 <- PND(t, rho, lambda, delta)
  return(t1 + (1 - t1) * s * t2 / (1 - (1 - t2) * s))
}

# PGF N
pgfN <- function(t0,t,rho0,rho,lambda,delta,s){
  t1 <- PND(t0, rho0, lambda, delta)
  t2 <- s ^ t0
  t5 <- s ^ t
  t7 <- pgfK(t0, rho0, lambda, delta, s)
  t8 <- pgfK(t, rho, lambda, delta, s)
  t9 <- pgfM(t, t, rho, rho, lambda, delta, t8)
  return(t1 * t2 + (1 - t1) * t5 * t7 * t9)
}

# pgf for K
pgfK <- function(t, rho, lambda, delta, s) {
  tmp <- 0
  for (i in 1:t) {
    tmp <- tmp + Pk(i, rho, lambda, delta) * s^i
  }
  val <- tmp / (1 - PND(t, rho, lambda, delta))
  return(val)
}


# PK
Pk <- function(t, rho, lambda, delta) {
  if (t == 1) {
    val <- 1 - PND(1, rho, lambda, delta)
  } else {
    val <- PND(t - 1, rho, lambda, delta) * (1 - psi_Xk(t - 1, rho, lambda, delta, psi_X(lambda, 1 - delta)))
  }
  return(val)
}

# PND - probability of not detecting
PND <- function(t, rho, lambda, delta) {
  if (t == 1) {
    val <- psi_Xk(0, rho, lambda, delta, psi_X(lambda, 1 - delta))
  } else {
    val <- psi_Xk(t - 1, rho, lambda, delta, psi_X(lambda, 1 - delta)) * PND(t - 1, rho, lambda, delta)
  }
}

#
psi_X <- function(lambda, s) {
  val <- exp(lambda * (s - 1))
}

# EXN
ExN <- function(t0, t, rho0,rho, lambda, delta) {
  f <- function(s) pgfN(t0, t, rho0,rho, lambda, delta, s)
  val = derivative(f = f, var = c(s=1), order = 1)
  return(val)
}

# var N
VarN <- function(t0, t, rho0,rho, lambda, delta) {
  f <- function(s) pgfN(t0, t, rho0,rho, lambda, delta, s)
  dev2 = derivative(f = f, var = c(s=1), order = 2)
  val = dev2+ExN(t0, t, rho0,rho, lambda, delta)- ExN(t0, t, rho0,rho, lambda, delta)^2
  return(val)
}

# ExM
ExM <- function(t0, t, rho0,rho, lambda, delta) {
  t1 <- PND(t0, rho0, lambda, delta)
  t2 <- PND(t, rho, lambda, delta)
  val <- (1 - t1) / (t2)
  return(val)
}

# var M
VarM <- function(t0, t, rho0,rho, lambda, delta) {
  f <- function(s) pgfM(t0, t, rho0,rho, lambda, delta, s)
  dev2 = derivative(f = f, var = c(s=1), order = 2)
  val = dev2+ExM(t0, t, rho0,rho, lambda, delta)- ExM(t0, t, rho0,rho, lambda, delta)^2
  return(val)
}

#  mean K
ExK <- function(t, rho, lambda, delta) {
  tmp <- 0
  for (i in 1:t) {
    tmp <- tmp + Pk(i, rho, lambda, delta) * i
  }
  val <- tmp / (1 - PND(t, rho, lambda, delta))
  return(val)
}

# var K
VarK <- function(t, rho, lambda, delta) {
  f <- function(s) pgfK(t, rho, lambda, delta, s)
  dev2 = derivative(f = f, var = c(s=1), order = 2)
  val = dev2+ExK(t, rho, lambda, delta)- ExK(t, rho, lambda, delta)^2
  return(val)
}

psi_Xk <- function(t, rho, lambda, delta, s) {
  if (t == 0) {
    val <- psi_X(rho, s)
  } else {
    val <- psi_Xk(t - 1, rho, lambda, delta, psi_X(lambda, (1 - delta) * s)) / psi_Xk(t - 1, rho, lambda, delta, psi_X(lambda, 1 - delta))
  }
  return(val)
}
