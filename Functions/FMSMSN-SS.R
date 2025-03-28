
#=========================packages========================

library(mnormt)
library(sn)
library(TSA)
library(tmvtnorm)
library(coda)
library(MASS)
library(MCMCpack)
library(moments)
library(mixsmsn)
library(optimx)

#==========================function==========================
garchfss= function(n,ns){
  ht=yt=rep(0,0)
  zt <-rmix((n+ns), p=piip, family="Skew.slash", arg=list(arg1p,arg2p))
  ht[1] = wp/(1-alpha1p-beta1p)
  yt[1] = wp+zt[1]*sqrt(ht[1]) 
  for(i in 2:(n+ns)){
    ht[i]= mu0p+alpha1p*(yt[i-1]-wp)^2+beta1p*ht[i-1]
    yt[i] = wp+zt[i]*sqrt(ht[i])
  }
  list(yt,zt,ht)
}

H2hat=function(w,m,a,b){
  H2=rep(0,0)
  H2[1]=m+a*(x0-w)^2+b*H0
  for(i in 2:length(x)) H2[i]=m+a*(x[i-1]-w)^2+b*H2[i-1]
  H2
}

ERROR=function(w,m,a,b){
  H2=H2hat(w,m,a,b)
  (x-w)/sqrt(H2)
}


dSS <- function(y, mu, sigma2, shape,nu){
  resp <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y)) {
    f <- function(u) 2*nu*u^(nu - 1)*dnorm(y[i],mu,sqrt(sigma2/u))*pnorm(u^(1/2)*shape*(sigma2^(-1/2))*(y[i]-mu))
    resp[i] <- integrate(f,0,1)$value
  }
  return(resp)
}

d.mixedSS <- function(x, pi1, mu, sigma2, shape, nu){
  g <- length(pi1)
  dens <- 0
  for (j in 1:g) dens <- dens + pi1[j]*dSS(x, mu[j], sigma2[j], shape[j], nu)
  return(dens)
}

estss<-function(x,x0,g,error,iter.max){
  dSS <- function(y, mu, sigma2, shape,nu){
    resp <- vector(mode = "numeric", length = length(y))
    for (i in 1:length(y)) {
      f <- function(u) 2*nu*u^(nu - 1)*dnorm(y[i],mu,sqrt(sigma2/u))*pnorm(u^(1/2)*shape*(sigma2^(-1/2))*(y[i]-mu))
      resp[i] <- integrate(f,0,1,stop.on.error = FALSE)$value
    }
    return(resp)
  }
  d.mixedSS <- function(x, pi1, mu, sigma2, shape, nu){
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dSS(x, mu[j], sigma2[j], shape[j], nu)
    return(dens)
  }  
  n=len=length(x)
  HH=rep(0,0)
  for(i in 1:n) HH[i]=x[i]^2
  H0=HH[1]
  ER=ERROR(0.1,0.1,0.1,0.5)
  st=smsn.mix(ER, nu = 1, g = g, get.init = TRUE, criteria = TRUE,group = TRUE, family = "Skew.t", calc.im=FALSE)
  nu=st$nu
  p=st$pii[1]
  pii=st$pii;mu=st$mu;sigma2=st$sigma2;shape=st$shape
  LLh<-function(parm){
    htt=rep(0,0)
    htt[1]<-exp(parm[1])+(1/(1+exp(-parm[2])))*(x0-parm[4])^2+((1-(1/(1+exp(-parm[2]))))/(1+exp(-parm[3])))*H0
    SS=log(d.mixedSS((x[1]-parm[4])/sqrt(htt[1]),pii,mu,sigma2,shape,nu))
    for(i in 2:n){
      htt[i]<-exp(parm[1])+(1/(1+exp(-parm[2])))*(x[i-1]-parm[4])^2+((1-(1/(1+exp(-parm[2]))))/(1+exp(-parm[3])))*htt[i-1]
      SS=SS+log(d.mixedSS((x[i]-parm[4])/sqrt(htt[i]),pii,mu,sigma2,shape,nu))
    }
    -SS
  }
  LLhest=0
  try(LLhest<-optimx(c(rnorm(1,-10,0.5),rnorm(3,0,1)),LLh,method="L-BFGS-B"), silent = TRUE)
  while(length(LLhest)==1|is.na(LLhest[1])) try(LLhest<-optimx(c(rnorm(1,-10,0.5),rnorm(3,0,0.5)),LLh,method="L-BFGS-B"), silent = TRUE)
  parm0=c(LLhest$p1,LLhest$p2,LLhest$p3,LLhest$p4)
  mu0=exp(LLhest$p1)
  alpha1=(1/(1+exp(-LLhest$p2)))
  beta1=((1-(1/(1+exp(-LLhest$p2))))/(1+exp(-LLhest$p3)))
  w=LLhest$p4
  HH=H2hat(w,mu0,alpha1,beta1)
  ER=ERROR(w,mu0,alpha1,beta1)
  lk <- sum(log(d.mixedSS(ER, pii, mu, sigma2, shape, nu)))
  n <- length(ER)
  delta <- Delta <- Gama <- rep(0, g)
  for (k in 1:g) {
    delta[k] <- shape[k]/(sqrt(1 + shape[k]^2))
    Delta[k] <- sqrt(sigma2[k]) * delta[k]
    Gama[k] <- sigma2[k] - Delta[k]^2
  }
  teta <- c(mu, Delta, Gama, pii, nu)
  mu.old <- mu
  Delta.old <- Delta
  Gama.old <- Gama
  criterio <- 1
  count <- 0
  while ((criterio > error) & (count <= iter.max)) {
    count <- count + 1
    tal <- matrix(0, n, g)
    S1 <- matrix(0, n, g)
    S2 <- matrix(0, n, g)
    S3 <- matrix(0, n, g)
    for (j in 1:g) {
      u <- vector(mode = "numeric", length = n)
      E <- vector(mode = "numeric", length = n)
      dj <- ((ER - mu[j])/sqrt(sigma2[j]))^2
      if (length(which(dj < 1e-10)) > 0) 
        dj[which(dj < 1e-10)] <- 1e-10
      Mtij2 <- 1/(1 + (Delta[j]^2) * (Gama[j]^(-1)))
      Mtij <- sqrt(Mtij2)
      mutij <- Mtij2 * Delta[j] * (Gama[j]^(-1)) * 
        (ER - mu[j])
      A <- mutij/Mtij
      for (i in 1:n) {
        d <- dSS(ER[i], mu[j], sigma2[j], shape[j], 
                 nu)
        if (d < 1e-10) 
          d <- 1e-10
        E[i] <- exp((nu + 1) * log(2) + log(nu) + lgamma(nu + 
                                                           1) - (log(d) + log(pi) + 0.5 * log(sigma2[j])) - 
                      (nu + 1) * log(dj[i] + A[i]^2) + pgamma(1, 
                                                              nu + 1, (dj[i] + A[i]^2)/2, log.p = TRUE))
        faux <- function(u) exp((nu + 0.5) * log(u) - 
                                  u * dj[i]/2 + pnorm(u^(1/2) * A[i], log.p = TRUE))
        aux22 <- integrate(faux, 0, 1)$value
        u[i] <- ((sqrt(2) * nu)/(d * sqrt(pi) * sqrt(sigma2[j]))) * 
          aux22
      }
      d1 <- dSS(ER, mu[j], sigma2[j], shape[j], nu)
      if (length(which(d1 == 0)) > 0) 
        d1[which(d1 == 0)] <- .Machine$double.xmin
      d2 <- d.mixedSS(ER, pii, mu, sigma2, shape, nu)
      if (length(which(d2 == 0)) > 0) 
        d2[which(d2 == 0)] <- .Machine$double.xmin
      tal[, j] <- d1 * pii[j]/d2
      S1[, j] <- tal[, j] * u
      S2[, j] <- tal[, j] * (mutij * u + Mtij * E)
      S3[, j] <- tal[, j] * (mutij^2 * u + Mtij2 + 
                               Mtij * mutij * E)
      pii[j] <- (1/n) * sum(tal[, j])
      mu[j] <- sum(S1[, j] * ER - Delta.old[j] * S2[, 
                                                    j])/sum(S1[, j])
      Delta[j] <- sum(S2[, j] * (ER - mu[j]))/sum(S3[, 
                                                     j]);if(Delta[j]<0.001)
                                                       Delta[j]=runif(1,0.001,0.01)
      Gama[j] <- ifelse(sum(S1[, j] * (ER - mu[j])^2 - 
                              2 * (ER - mu[j]) * Delta[j] * S2[, j] + Delta[j]^2 * 
                              S3[, j])/sum(tal[, j]) == 0, .Machine$double.xmin^0.5, 
                        sum(S1[, j] * (ER - mu[j])^2 - 2 * (ER - mu[j]) * 
                              Delta[j] * S2[, j] + Delta[j]^2 * S3[, j])/sum(tal[, 
                                                                                 j]))
      sigma2[j] <- Gama[j] + Delta[j]^2
      shape[j] <- ((sigma2[j]^(-1/2)) * Delta[j])/(sqrt(1 - 
                                                          (Delta[j]^2) * (sigma2[j]^(-1))))
    }
    logvero.SS <- function(nu) sum(log(d.mixedSS(ER, pii, 
                                                 mu, sigma2, shape, nu)))
    nu <- optimize(logvero.SS, c(0, 100), tol = 1e-06, 
                   maximum = TRUE)$maximum
    pii[g] <- 1 - (sum(pii) - pii[g])
    zero.pos <- NULL
    zero.pos <- which(pii == 0)
    if (length(zero.pos) != 0) {
      pii[zero.pos] <- 1e-10
      pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
    }
    param <- teta
    teta <- c(mu, Delta, Gama, pii, nu)
    lk1 <- sum(log(d.mixedSS(ER, pii, mu, sigma2, shape, 
                             nu)))
    criterio <- abs(lk1/lk - 1)
    mu.old <- mu
    Delta.old <- Delta
    Gama.old <- Gama
    lk <- lk1
  }
  Ht=H2hat(w,mu0,alpha1,beta1)
  ER=ERROR(w,mu0,alpha1,beta1)
  LLhest<- optimx(parm0,LLh,method="L-BFGS-B")
  parm0=c(LLhest$p1,LLhest$p2,LLhest$p3,LLhest$p4)
  mu0=exp(LLhest$p1)
  alpha1=(1/(1+exp(-LLhest$p2)))
  beta1=((1-(1/(1+exp(-LLhest$p2))))/(1+exp(-LLhest$p3)))
  w=LLhest$p4
  Ht=H2hat(w,mu0,alpha1,beta1)
  ER=ERROR(w,mu0,alpha1,beta1)
  lk<- sum(log(d.mixedSS(ER, pii, mu, sigma2, shape,   nu))) 
  cl <- apply(tal, 1, which.max)
  icl <- 0
  for (j in 1:g) icl <- icl + sum(log(pii[j] * dSS(ER[cl == 
                                                        j], mu[j], sigma2[j], shape[j], nu)))
  d <- g * 3 + (g - 1)+1+4
  aic <- -2 * lk + 2 * d
  bic <- -2 * lk + log(n) * d
  edc <- -2 * lk + 0.2 * sqrt(n) * d
  icl <- -2 * icl + log(n) * d
  xp=w+ER*sqrt(Ht)
  obj.out <- list(mu = mu, sigma2 = sigma2, shape = shape,w=w,alpha1=alpha1,beta1=beta1,mu0=mu0, 
                  pii = pii, nu = nu, lk=lk, aic = aic, bic = bic, edc = edc, 
                  icl = icl, iter = count, n = length(ER), group = cl,ER=ER,Ht=Ht,xeq=sum(x!=xp),xdif=mean(abs(x-xp)))
  obj.out
}

goodss<-function(ER,g, pii, mu, sigma2, shape, nu){
  dSS <- function(y, mu, sigma2, shape,nu){
    resp <- vector(mode = "numeric", length = length(y))
    for (i in 1:length(y)) {
      f <- function(u) 2*nu*u^(nu - 1)*dnorm(y[i],mu,sqrt(sigma2/u))*pnorm(u^(1/2)*shape*(sigma2^(-1/2))*(y[i]-mu))
      resp[i] <- integrate(f,0,1,stop.on.error = FALSE)$value
    }
    return(resp)
  }
  d.mixedSS <- function(x, pi1, mu, sigma2, shape, nu){
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dSS(x, mu[j], sigma2[j], shape[j], nu)
    return(dens)
  }  
  lk<- sum(log(d.mixedSS(ER, pii, mu, sigma2, shape,   nu))) 
  d <- g * 3 + (g - 1)+1+4
  aic <- -2 * lk + 2 * d
  bic <- -2 * lk + log(n) * d
  edc <- -2 * lk + 0.2 * sqrt(n) * d
  obj.out <-c(lk=lk, aic = aic, bic = bic, edc = edc)
  obj.out
}


garchperss= function(n,xx,piip,arg,wp,mu0p,alpha1p,beta1p){
  ht=yt=rep(0,0)
  zt <-rmix((n), p=piip, family="Skew.slash", arg=arg)
  ht[1] = xx
  yt[1] = wp+zt[1]*sqrt(ht[1]) 
  for(i in 2:(n)){
    ht[i]= mu0p+alpha1p*(yt[i-1]-wp)^2+beta1p*ht[i-1]
    yt[i] = wp+zt[i]*sqrt(ht[i])
  }
  list(yt=yt,zt=zt,ht=ht)
}


