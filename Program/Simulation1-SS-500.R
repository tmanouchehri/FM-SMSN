rm(list=ls())
start.time <- Sys.time()
ptm<-proc.time()

#==========================import function=====================
WD.PATH = paste(getwd(),"/Functions", sep = "")

source(paste(WD.PATH, '/FMSMSN-SS.r', sep = ""))

#==========================parameters==========================
#Moderately Components
mu0p=10^(-6)
alpha1p=0.05
beta1p=0.9
wp=10^(-4)
mu2p <- -2
sigma2.2p <- 0.5
lambda1p <- 1
lambda2p <- -2
delta1p=lambda1p/sqrt(lambda1p^2+1);delta2p=lambda2p/sqrt(lambda2p^2+1)
nup =5
piip <- c(0.6,0.4)
sigma2.1p<-(1-piip[2]*((nup/(nup-2))-delta2p^2*(nup/(nup-1))^2*(2/pi))*sigma2.2p)/(piip[1]*((nup/(nup-2))-delta1p^2*(nup/(nup-1))^2*(2/pi)))
sigma2p <- c(sigma2.1p,sigma2.2p)
shapep <- c(lambda1p,lambda2p)
deltap<-c(delta1p,delta2p)
mu1p<-(-sum(piip*sqrt(2/pi)*deltap*sqrt(sigma2p)*(nup/(nup-1)))-piip[2]*mu2p)/piip[1]
mup <- c(mu1p,mu2p)
arg1p = c(mu1p, sigma2.1p, lambda1p, nup)
arg2p = c(mu2p, sigma2.2p, lambda2p, nup)
Delta1p=sqrt(sigma2.1p)*delta1p;Delta2p=sqrt(sigma2.2p)*delta2p
Ro1p=sigma2.1p-Delta1p^2;Ro2p=sigma2.2p-Delta2p^2
ratp<-c(mu1p,mu2p,sigma2.1p,sigma2.2p,lambda1p,lambda2p,mu0p,alpha1p,beta1p,piip,nup,wp)

must<-function(m,s,d,n)  m+sqrt(2/pi)*s*d*(n/(n-1))
ssst<-function(m,s,d,n)  ((n/(n-2))-(2/pi)*(n/(n-1))^2*d^2)*s^2

piip[1]*must(mu1p,sqrt(sigma2.1p),delta1p,nup)+(1-piip[1])*must(mu2p,sqrt(sigma2.2p),delta2p,nup)
piip[1]*ssst(mu1p,sqrt(sigma2.1p),delta1p,nup)+(1-piip[1])*ssst(mu2p,sqrt(sigma2.2p),delta2p,nup)

#==========================n=500==========================
R=500
n=500
ns=n*20/100
SS=matrix(rep(0,R*13),nrow=13)

#==========================histogram==========================
Data=garchfss(n,ns)
zss=Data[[2]][(ns+1):(n+ns)]
xfitss<-seq(min(zss),max(zss),length=400)
yfitss<-d.mixedSS(xfitss, piip, mup, sigma2p, shapep, nup)

X11()
hist(zss, breaks = seq(min(zss),max(zss),length=70), freq = FALSE, main = "FMSMSN-SS", xlab = "Z")
lines(xfitss,yfitss , lwd = 2)

#==========================program==========================
for(r in 1:R){
  Data=garchfss(n,ns)
  x=Data[[1]][(ns+1):(n+ns)]
  x0=Data[[1]][ns]
  n=len=length(x)
  HH=rep(0,0)
  for(i in 1:n) HH[i]=x[i]^2
  H0=HH[1]
  resutfmsn_ss=estss(x,x0,2,10^(-5),100)
  if(resutfmsn_ss$pii[1]>resutfmsn_ss$pii[2]){
    SS[1,r]=resutfmsn_ss$mu[1]
    SS[2,r]=resutfmsn_ss$mu[2]
    SS[3,r]=resutfmsn_ss$sigma2[1]
    SS[4,r]=resutfmsn_ss$sigma2[2]
    SS[5,r]=resutfmsn_ss$shape[1]
    SS[6,r]=resutfmsn_ss$shape[2]
    SS[7,r]=resutfmsn_ss$mu0
    SS[8,r]=resutfmsn_ss$alpha1
    SS[9,r]=resutfmsn_ss$beta1
    SS[10,r]=resutfmsn_ss$pii[1]
    SS[11,r]=resutfmsn_ss$pii[2]
    SS[12,r]=resutfmsn_ss$nu  
    SS[13,r]=resutfmsn_ss$w  
  }else
  {
    SS[1,r]=resutfmsn_ss$mu[2]
    SS[2,r]=resutfmsn_ss$mu[1]
    SS[3,r]=resutfmsn_ss$sigma2[2]
    SS[4,r]=resutfmsn_ss$sigma2[1]
    SS[5,r]=resutfmsn_ss$shape[2]
    SS[6,r]=resutfmsn_ss$shape[1]
    SS[7,r]=resutfmsn_ss$mu0
    SS[8,r]=resutfmsn_ss$alpha1
    SS[9,r]=resutfmsn_ss$beta1
    SS[10,r]=resutfmsn_ss$pii[2]
    SS[11,r]=resutfmsn_ss$pii[1]
    SS[12,r]=resutfmsn_ss$nu  
    SS[13,r]=resutfmsn_ss$w 
  }
}


write.csv(SS,"Sim1_ss500.csv", sep = ",", col.names = NA,qmethod = "double")


proc.time()-ptm
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken 
