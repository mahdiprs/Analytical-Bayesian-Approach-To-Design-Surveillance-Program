library(symengine)  # for symbolic derivatives

# R functions for the survey block paper

calc_moments<- function(pgf, parms) {
  # calculate first two moments (i.e mean, variance)
  require(symengine)
  e<- S(pgf)
  parms$s<- 1
  nam<- names(parms)
  args<- list()
  for(i in 1: length(parms)) {
    assign(nam[i],parms[[i]])
    args[[i]]<- S(eval(nam[i]))
  }
  args<- Vector(args)
  d1<- symengine::D(e, "s")
  d2<- symengine::D(d1, "s")
  fd1<- as.function(d1, args=args)
  fd2<- as.function(d2, args=args)
  ex<- do.call(fd1, parms)
  vx<- do.call(fd2, parms) + ex - ex^2
  return(list(Ex=ex,Vx=vx))
}

calc_pmf<- function(pgf, parms, support) {
  # Given PGF, calculate PMF for the support
  require(symengine)
  probs<- rep(NA, length(support))
  parms$s<- 0
  nam<- names(parms)
  args<- list()
  for(i in 1: length(parms)) {
    assign(nam[i],parms[[i]])
    args[[i]]<- S(eval(nam[i]))
  }
  args<- Vector(args)
  dn<- list()
  e<- symengine::S(pgf)
  dn[[1]]<- e
  fe<- as.function(e, args=args)
  probs[1]<- do.call(fe, parms)
  for(i in 2:length(support)) {
    dn[[i]]<- symengine::D(dn[[i-1]], "s")
    fe<- as.function(dn[[i]], args=args)
    ev<- do.call(fe, parms)
    lprob<- log(ev) - lfactorial(support[i])
    probs[i]<- exp(lprob)
  }
  data.frame(n=support,p=probs)
}

# Construct the various PGF's as expressions to facilitate symbolic derivatives

# kappa
kappa <- function(rho, lambda, delta, tr, p0=FALSE) {
 # calculate kappa for a given probability of absence (tr)
  t <- 1
  if(p0){
    pgf<- make_psi_Xk(t, p0=TRUE)
    f<- as.function(S(pgf))
    P <- f(s=0, rho0=rho, lambda=lambda, delta=delta)
  while (P < tr) {
    t <- t + 1
    pgf<- make_psi_Xk(t, p0=TRUE)
    f<- as.function(S(pgf))
    P <- f(s=0, rho0=rho, lambda=lambda, delta=delta)
    }
  }
  else{
    pgf<- make_psi_Xk(t, p0=FALSE)
    f<- as.function(S(pgf))
    P <- f(s=0, rho=rho, lambda=lambda, delta=delta)
    while (P < tr) {
      t <- t + 1
      pgf<- make_psi_Xk(t, p0=FALSE)
      f<- as.function(S(pgf))
      P <- f(s=0, rho=rho, lambda=lambda, delta=delta)
    }
  }
  return(t)
}


# PGF helper functions

psi_Xk_num<- function(k, p0) {
  if(p0) {
    if(k==0){
      num<- "exp(rho0*(s-1))"
    }
    else{
      num<- paste0("exp(rho0*(",psi_X_num(k-1),"-1))")
    }
  }
  else{
    if(k==0){
      num<- "exp(rho*(s-1))"
    }
    else{
      num<- paste0("exp(rho*(",psi_X_num(k-1),"-1))")
    }
  }
  num
}

psi_X_num<- function(k) {
  if(k==0){
    num<- "exp(lambda*((1-delta)*s-1))"
  }
  else {
    num<- paste0("exp(lambda*((1-delta)*",psi_X_num(k-1),"-1))")
  }
  num
}

psi_Xk_den<- function(k, p0) {
  if(p0) {
    if(k==0){
      num<- NULL
    }
    else{
      num<- paste0("exp(rho0*(",psi_X_den(k-1),"))")
    }
  }
  else{
    if(k==0){
      num<- NULL
    }
    else{
      num<- paste0("exp(rho*(",psi_X_den(k-1),"))")
    }
  }
  num
}

psi_X_den<- function(k) {
  if(k==0){
    num<- "exp(-lambda*delta)-1"
  }
  else {
    num<- paste0("exp(lambda*((1-delta)*",psi_X_den(k-1),"))-1")
  }
  num
}

make_psi_Xk<- function(k, p0) {
  if(k==0) {
    pgf<- psi_Xk_num(k=k, p0=p0)
  }
  else{
    num<- psi_Xk_num(k=k, p0=p0)
    denom<- psi_Xk_den(k=k, p0=p0)
    pgf<- paste0(num,"/",denom)
  }
  pgf
}


make_PND<- function(k, p0) {
  if(p0) {
    if(k==1){
      num<- "exp(rho0*(exp(-lambda*delta)-1))"
    }
    else{
      num<- paste0("exp(rho0*(",psi_X_den(k-1),"))")
    }
  }
  else{
    if(k==1){
      num<- "exp(rho*(exp(-lambda*delta)-1))"
    }
    else{
      num<- paste0("exp(rho*(",psi_X_den(k-1),"))")
    }
  }
  num
}


psi_D_num<- function(k) {
  if(k==0){
    num<- "exp(lambda*((1-delta)*exp(-lambda*delta)-1))"
  }
  else {
    num<- paste0("exp(lambda*((1-delta)*",psi_D_num(k-1),"-1))")
  }
  num
}

psi_Dk_num<- function(k, p0) {
  if(p0) {
    if(k==0){
      num<- "exp(rho0*(exp(-lambda*delta)-1))"
    }
    else{
      num<- paste0("exp(rho0*(",psi_D_num(k-1),"-1))")
    }
  }
  else{
    if(k==0){
      num<- "exp(rho*(exp(-lambda*delta)-1))"
    }
    else{
      num<- paste0("exp(rho*(",psi_D_num(k-1),"-1))")
    }
  }
  num
}

psi_Dk<- function(k, p0) {
  if(k==0) {
    val<- psi_Dk_num(k=k, p0=p0)
  }
  else{
    num<- psi_Dk_num(k=k, p0=p0)
    denom<- psi_Xk_den(k=k, p0=p0)
    val<- paste0(num,"/",denom)
  }
  val
}


# PK
make_Pk <- function(k, p0) {
  if (k == 1) {
    val <- paste0("1 - ",make_PND(1, p0))
  } else {
    val <- paste0(make_PND(k-1, p0)," * ", "(1 - ",psi_Dk(k-1, p0),")")
  }
  return(val)
}

# Build the various PGF's

make_pgfK<- function(k, k0=NULL){
  tmp<- list()
  if(!is.null(k0)) {
    for(i in 1:k0){
      tmp[[i]]<- paste0("(",make_Pk(i, p0=TRUE),") * s^",i)
    }
    den<- paste0("1- ",make_PND(k0, p0=TRUE))
    num<- paste0(tmp, collapse="+")
  }else {
    for(i in 1:k){
      tmp[[i]]<- paste0("(",make_Pk(i, p0=FALSE),") * s^",i)
    }
      den<- paste0("1- ",make_PND(k, p0=FALSE))
      num<- paste0(tmp, collapse="+")
    }
  return(paste0("(",num,")/(",den,")"))
}


make_pgfM <- function(k, k0=NULL) {
  if(!is.null(k0)){
    t1 <- make_PND(k0, p0=TRUE)
    t2 <- paste0("(1-",t1,") * s")
    t3 <- make_PND(k, p0=FALSE)
    t4 <- paste0(t3," / ","(1 - (1 - ",t3,") * s)")
    val <- paste0(t1," + ",t2," * ",t4)
  }
  else{
    t1 <- make_PND(k, p0=FALSE)
    val <- paste0(t1," / ","(1 - (1 - ",t1,") * s)")
  }
  return(val)
}

make_pgfMk <- function(k) {
  t1 <- make_PND(k, p0=FALSE)
  t2<- make_pgfK(k)
  val <- paste0(t1," / ","(1 - (1 - ",t1,") * ",t2,")")
  return(val)
}

make_pgfN <- function(k, k0=NULL) {
  if(!is.null(k0)) {
    t1<- make_PND(k0, p0=TRUE)
    t2<- paste0("s^",k0)
    t3 <- paste0("s^",k)
    t4<- make_pgfK(k, k0)
    t5<- make_pgfMk(k)
    val <- paste0(t1,"*",t2,"+","(1-",t1,")","*",t4,"*",t3,"*",t5)
  } else {
    t1 <- paste0("s^",k)
    t2 <- make_pgfMk(k)
    val <- paste0("(",t1,")","*","(",t2,")")
  }
  return(val)
}
