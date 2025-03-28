rm(list=ls())
start.time <- Sys.time()
ptm<-proc.time()

#==========================import function==========================
WD.PATH = paste(getwd(),"/Functions", sep = "")

source(paste(WD.PATH, '/FMSMSN-ST.r', sep = ""))
source(paste(WD.PATH, '/FMSMSN-SN.r', sep = ""))
source(paste(WD.PATH, '/FMSMSN-SS.r', sep = ""))
source(paste(WD.PATH, '/FMSMSN-SCN.r', sep = ""))

#==========================parameters st==========================
mu0p=10^(-6)
alpha1p=0.05
beta1p=0.9
wp=10^(-4)
mu2p <- -1
sigma2.2p <- 2
lambda1p <- 1
lambda2p <- -2
delta1p=lambda1p/sqrt(lambda1p^2+1);delta2p=lambda2p/sqrt(lambda2p^2+1)
nup = 5
piip <- c(0.6,0.4)
sigma2.1p<-(1-(piip[2]*((nup/(nup-2))-(nup/pi)*(gamma((nup-1)/2)/gamma((nup)/2))^2*delta2p^2)*sigma2.2p))/(piip[1]*((nup/(nup-2))-(nup/pi)*(gamma((nup-1)/2)/gamma((nup)/2))^2*delta1p^2))
sigma2p <- c(sigma2.1p,sigma2.2p)
shapep <- c(lambda1p,lambda2p)
deltap<-c(delta1p,delta2p)
mu1p<-(-sum(piip*sqrt(nup)*pi*(gamma((nup-1)/2)/gamma((nup)/2))*deltap*sqrt(sigma2p))-piip[2]*mu2p)/piip[1]
mup <- c(mu1p,mu2p)
arg1p = c(mu1p, sigma2.1p, lambda1p, nup)
arg2p = c(mu2p, sigma2.2p, lambda2p, nup)
Delta1p=sqrt(sigma2.1p)*delta1p;Delta2p=sqrt(sigma2.2p)*delta2p
Ro1p=sigma2.1p-Delta1p^2;Ro2p=sigma2.2p-Delta2p^2
ratp<-c(mu1p,mu2p,sigma2.1p,sigma2.2p,lambda1p,lambda2p,mu0p,alpha1p,beta1p,piip,nup,wp)

must<-function(m,s,d,n)  m+sqrt(n)*pi*(gamma((n-1)/2)/gamma(n/2))*s*d
ssst<-function(m,s,d,n)  ((n/(n-2))-(n/pi)*(gamma((n-1)/2)/gamma((n)/2))^2*d^2)*s^2

piip[1]*must(mu1p,sqrt(sigma2.1p),delta1p,nup)+(1-piip[1])*must(mu2p,sqrt(sigma2.2p),delta2p,nup)
piip[1]*ssst(mu1p,sqrt(sigma2.1p),delta1p,nup)+(1-piip[1])*ssst(mu2p,sqrt(sigma2.2p),delta2p,nup)

#==========================size==========================
R=500
n=250
ns=n*20/100

#====================================================
SS=matrix(rep(0,R*12*4),ncol=12*4)
resutfmsn_st1=resutfmsn_st2=resutfmsn_st3=resutfmsn_sn1=resutfmsn_sn2=resutfmsn_sn3=list()
resutfmsn_ss1=resutfmsn_ss2=resutfmsn_ss3=resutfmsn_scn1=resutfmsn_scn2=resutfmsn_scn3=list()

#====================program===================
for(r in 1:R){
  Data=garchfst(n,ns)
  x=Data[[1]][(ns+1):(n+ns)]
  x0=Data[[1]][ns]
  n=len=length(x)
  HH=rep(0,0)
  for(i in 1:n) HH[i]=x[i]^2
  H0=HH[1]
  resutfmsn_sn1[[r]]<-estsn(x,x0,1,10^(-5),100)
  resutfmsn_sn2[[r]]<-estsn(x,x0,2,10^(-5),100)
  resutfmsn_sn3[[r]]<-estsn(x,x0,3,10^(-5),100)
  resutfmsn_st1[[r]]<-estst(x,x0,1,10^(-5),100)
  resutfmsn_st2[[r]]<-estst(x,x0,2,10^(-5),100)
  resutfmsn_st3[[r]]<-estst(x,x0,3,10^(-5),100) 
  resutfmsn_ss1[[r]]<-estss(x,x0,1,10^(-5),100)
  resutfmsn_ss2[[r]]<-estss(x,x0,2,10^(-5),100)
  resutfmsn_ss3[[r]]<-estss(x,x0,3,10^(-5),100) 
  resutfmsn_scn1[[r]]<-estscn(x,x0,1,10^(-5),100) 
  resutfmsn_scn2[[r]]<-estscn(x,x0,2,10^(-5),100)
  resutfmsn_scn3[[r]]<-estscn(x,x0,3,10^(-5),100) 

res_sn1=c(resutfmsn_sn1[[r]]$lk, resutfmsn_sn1[[r]]$aic, resutfmsn_sn1[[r]]$bic, resutfmsn_sn1[[r]]$edc)
res_sn2=c(resutfmsn_sn2[[r]]$lk, resutfmsn_sn2[[r]]$aic, resutfmsn_sn2[[r]]$bic, resutfmsn_sn2[[r]]$edc)
res_sn3=c(resutfmsn_sn3[[r]]$lk, resutfmsn_sn3[[r]]$aic, resutfmsn_sn3[[r]]$bic, resutfmsn_sn3[[r]]$edc)

res_st1=c(resutfmsn_st1[[r]]$lk, resutfmsn_st1[[r]]$aic, resutfmsn_st1[[r]]$bic, resutfmsn_st1[[r]]$edc)
res_st2=c(resutfmsn_st2[[r]]$lk, resutfmsn_st2[[r]]$aic, resutfmsn_st2[[r]]$bic, resutfmsn_st2[[r]]$edc)
res_st3=c(resutfmsn_st3[[r]]$lk, resutfmsn_st3[[r]]$aic, resutfmsn_st3[[r]]$bic, resutfmsn_st3[[r]]$edc)


res_ss1=c(resutfmsn_ss1[[r]]$lk, resutfmsn_ss1[[r]]$aic, resutfmsn_ss1[[r]]$bic, resutfmsn_ss1[[r]]$edc)
res_ss2=c(resutfmsn_ss2[[r]]$lk, resutfmsn_ss2[[r]]$aic, resutfmsn_ss2[[r]]$bic, resutfmsn_ss2[[r]]$edc)
res_ss3=c(resutfmsn_ss3[[r]]$lk, resutfmsn_ss3[[r]]$aic, resutfmsn_ss3[[r]]$bic, resutfmsn_ss3[[r]]$edc)

res_scn1=c(resutfmsn_scn1[[r]]$lk, resutfmsn_scn1[[r]]$aic, resutfmsn_scn1[[r]]$bic, resutfmsn_scn1[[r]]$edc)
res_scn2=c(resutfmsn_scn2[[r]]$lk, resutfmsn_scn2[[r]]$aic, resutfmsn_scn2[[r]]$bic, resutfmsn_scn2[[r]]$edc)
res_scn3=c(resutfmsn_scn3[[r]]$lk, resutfmsn_scn3[[r]]$aic, resutfmsn_scn3[[r]]$bic, resutfmsn_scn3[[r]]$edc)

SS[r,]=c(res_sn1,res_sn2,res_sn3,res_st1,res_st2,res_st3,res_ss1,res_ss2,res_ss3,res_scn1,res_scn2,res_scn3)

}

write.csv(SS,"sim2_250.csv", sep = ",", col.names = NA,qmethod = "double")


proc.time()-ptm
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken  
