rm(list=ls())
start.time <- Sys.time()
ptm<-proc.time()

#==========================import function==========================
WD.PATH = paste(getwd(),"/Functions", sep = "")

source(paste(WD.PATH, '/FMSMSN-ST.r', sep = ""))
source(paste(WD.PATH, '/FMSMSN-SN.r', sep = ""))
source(paste(WD.PATH, '/FMSMSN-SS.r', sep = ""))
source(paste(WD.PATH, '/FMSMSN-SCN.r', sep = ""))
source(paste(WD.PATH, '/FMSMSN-T.r', sep = ""))
source(paste(WD.PATH, '/FMSMSN-N.r', sep = ""))

#=======================import data and plot=======================
library(readxl)
mydata <-  read_excel("program/AAPL.xlsx")
head(mydata)
DJ <- data$Adj.Close
date <- as.Date(mydata[[1]], format = "%m/%d/%Y")
date <- date[!is.na(mydata[[6]])]
p=log(DJ)
rend.da=10^2*(p[2:1130]-p[1:1129])
date_rend <- date[-1]
Ntrain <- 904
Ntest <- 226
split_date <- date[Ntrain]
library(xts)
date <- mydata[[1]][AA:BB]
td <- as.Date(date, format="%Y-%m-%d")
Price_DJ <- xts(x=DJ, order.by=td)
Price_rend.da <- xts(x=rend.da, order.by=td)
df_price <- data.frame(Date = date, Price = DJ)
df_rend <- data.frame(Date = date_rend, Return = rend.da)


p1 <- ggplot(df_price, aes(x = Date, y = Price)) +
  geom_line(color = "black") +
  geom_vline(xintercept = as.numeric(split_date), color = "blue", linetype = "dashed", linewidth = 1) +
  labs(x = "Date", y = "Adjusted Close Price", title = "Apple Stock Price") +
  theme_minimal(base_size = 12) +
  scale_x_date(labels = date_format("%Y-%m"), date_breaks = "6 months") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


p2 <- ggplot(df_rend, aes(x = Date, y = Return)) +
  geom_line(color = "darkgreen") +
  geom_vline(xintercept = as.numeric(split_date), color = "blue", linetype = "dashed", linewidth = 1) +
  labs(x = "Date", y = "Log Return", title = "Daily Log Returns") +
  theme_minimal(base_size = 12) +
  scale_x_date(labels = date_format("%Y-%m"), date_breaks = "6 months") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

grid.arrange(p1, p2, ncol = 2)

op<-par(mfrow=c(3,1))
date <- mydata[[1]]
td <- as.Date(date, format="%m/%d/%Y")
Price_rend.da <- zoo(x=rend.da, order.by=td)
plot(Price_rend.da, xlab = "Time", ylab = "Returns", main = "Daily Returns of Apple Inc.", type = "l")
acf(ts(rend.da),200,main = "Daily Returns of Apple Inc.")
pacf(ts(rend.da),200,main = "Daily Returns of Apple Inc.")


par(mfcol=c(2,1))
hist(DJ, br=30, freq=F,xlab = "Adjust Closure Prices", ylab = "Frequency", main = "Apple Inc.")
qqnorm(DJ)
qqline(DJ)


par(mfcol=c(2,1))
hist(rend.da, br=200, freq=F,xlab = "Returns" ,main="Histogram Daily Returns of Apple Inc.")
curve(dnorm(x, mean=mean(rend.da), sd=sd(rend.da)), add=TRUE, min(rend.da), max(rend.da), col="red",lwd=2)
qqnorm(rend.da)
qqline(rend.da)

x=rend.da[2:Ntrain]
x0=rend.da[1]
n=len=length(x)
HH=rep(0,0)
for(i in 1:n) HH[i]=x[i]^2
H0=HH[1]

#=======================program=======================

resutfmsn_sn1<-estsn(x,x0,1,10^(-4),100)
resutfmsn_sn2<-estsn(x,x0,2,10^(-4),100)
resutfmsn_sn3<-estsn(x,x0,3,10^(-4),100)
resutfmsn_st1<-estst(x,x0,1,10^(-4),100)
resutfmsn_st2<-estst(x,x0,2,10^(-4),100)
resutfmsn_st3<-estst(x,x0,3,10^(-4),100)
resutfmsn_scn1<-estscn(x,x0,1,10^(-4),100) 
resutfmsn_scn2<-estscn(x,x0,2,10^(-4),100)
resutfmsn_scn3<-estscn(x,x0,3,10^(-4),100)
resutfmsn_ss1<-estss(x,x0,1,10^(-4),100)
resutfmsn_ss2<-estss(x,x0,2,10^(-4),100)
resutfmsn_ss3<-estss(x,x0,3,10^(-4),100)
resutfmsn_n1<-estn(x,x0,1,10^(-4),100)
resutfmsn_n2<-estn(x,x0,2,10^(-4),100)
resutfmsn_n3<-estn(x,x0,3,10^(-4),100)
resutfmsn_t1<-estt(x,x0,1,10^(-4),100)
resutfmsn_t2<-estt(x,x0,2,10^(-4),100)
resutfmsn_t3<-estt(x,x0,3,10^(-4),100)

PAR=cbind(c(resutfmsn_sn1$mu,0,0,resutfmsn_sn1$sigma2,0,0,resutfmsn_sn1$shape,0,0,resutfmsn_sn1$nu,0,resutfmsn_sn1$pii,0,0,resutfmsn_sn1$w,resutfmsn_sn1$mu0,resutfmsn_sn1$alpha1,resutfmsn_sn1$beta1),
          c(resutfmsn_sn2$mu,0,resutfmsn_sn2$sigma2,0,resutfmsn_sn2$shape,0,resutfmsn_sn2$nu,0,resutfmsn_sn2$pii,0,resutfmsn_sn2$w,resutfmsn_sn2$mu0,resutfmsn_sn2$alpha1,resutfmsn_sn2$beta1),
          c(resutfmsn_sn3$mu,resutfmsn_sn3$sigma2,resutfmsn_sn3$shape,resutfmsn_sn3$nu,0,resutfmsn_sn3$pii,resutfmsn_sn3$w,resutfmsn_sn3$mu0,resutfmsn_sn3$alpha1,resutfmsn_sn3$beta1),
          c(resutfmsn_st1$mu,0,0,resutfmsn_st1$sigma2,0,0,resutfmsn_st1$shape,0,0,resutfmsn_st1$nu,0,resutfmsn_st1$pii,0,0,resutfmsn_st1$w,resutfmsn_st1$mu0,resutfmsn_st1$alpha1,resutfmsn_st1$beta1),
          c(resutfmsn_st2$mu,0,resutfmsn_st2$sigma2,0,resutfmsn_st2$shape,0,resutfmsn_st2$nu,0,resutfmsn_st2$pii,0,resutfmsn_st2$w,resutfmsn_st2$mu0,resutfmsn_st2$alpha1,resutfmsn_st2$beta1),
          c(resutfmsn_st3$mu,resutfmsn_st3$sigma2,resutfmsn_st3$shape,resutfmsn_st3$nu,0,resutfmsn_st3$pii,resutfmsn_st3$w,resutfmsn_st3$mu0,resutfmsn_st3$alpha1,resutfmsn_st3$beta1),
          c(resutfmsn_ss1$mu,0,0,resutfmsn_ss1$sigma2,0,0,resutfmsn_ss1$shape,0,0,resutfmsn_ss1$nu,0,resutfmsn_ss1$pii,0,0,resutfmsn_ss1$w,resutfmsn_ss1$mu0,resutfmsn_ss1$alpha1,resutfmsn_ss1$beta1),
          c(resutfmsn_ss2$mu,0,resutfmsn_ss2$sigma2,0,resutfmsn_ss2$shape,0,resutfmsn_ss2$nu,0,resutfmsn_ss2$pii,0,resutfmsn_ss2$w,resutfmsn_ss2$mu0,resutfmsn_ss2$alpha1,resutfmsn_ss2$beta1),
          c(resutfmsn_ss3$mu,resutfmsn_ss3$sigma2,resutfmsn_ss3$shape,resutfmsn_ss3$nu,0,resutfmsn_ss3$pii,resutfmsn_ss3$w,resutfmsn_ss3$mu0,resutfmsn_ss3$alpha1,resutfmsn_ss3$beta1),
          c(resutfmsn_scn1$mu,0,0,resutfmsn_scn1$sigma2,0,0,resutfmsn_scn1$shape,0,0,resutfmsn_scn1$nu,resutfmsn_scn1$pii,0,0,resutfmsn_scn1$w,resutfmsn_scn1$mu0,resutfmsn_scn1$alpha1,resutfmsn_scn1$beta1),
          c(resutfmsn_scn2$mu,0,resutfmsn_scn2$sigma2,0,resutfmsn_scn2$shape,0,resutfmsn_scn2$nu,resutfmsn_scn2$pii,0,resutfmsn_scn2$w,resutfmsn_scn2$mu0,resutfmsn_scn2$alpha1,resutfmsn_scn2$beta1),
          c(resutfmsn_scn3$mu,resutfmsn_scn3$sigma2,resutfmsn_scn3$shape,resutfmsn_scn3$nu,resutfmsn_scn3$pii,resutfmsn_scn3$w,resutfmsn_scn3$mu0,resutfmsn_scn3$alpha1,resutfmsn_scn3$beta1),
          c(resutfmsn_n1$mu,0,0,resutfmsn_n1$sigma2,0,0,resutfmsn_n1$shape,0,0,resutfmsn_n1$nu,0,resutfmsn_n1$pii,0,0,resutfmsn_n1$w,resutfmsn_n1$mu0,resutfmsn_n1$alpha1,resutfmsn_n1$beta1),
          c(resutfmsn_n2$mu,0,resutfmsn_n2$sigma2,0,resutfmsn_n2$shape,0,resutfmsn_n2$nu,0,resutfmsn_n2$pii,0,resutfmsn_n2$w,resutfmsn_n2$mu0,resutfmsn_n2$alpha1,resutfmsn_n2$beta1),
          c(resutfmsn_n3$mu,resutfmsn_n3$sigma2,resutfmsn_n3$shape,resutfmsn_n3$nu,0,resutfmsn_n3$pii,resutfmsn_n3$w,resutfmsn_n3$mu0,resutfmsn_n3$alpha1,resutfmsn_n3$beta1),
          c(resutfmsn_t1$mu,0,0,resutfmsn_t1$sigma2,0,0,resutfmsn_t1$shape,0,0,resutfmsn_t1$nu,0,resutfmsn_t1$pii,0,0,resutfmsn_t1$w,resutfmsn_t1$mu0,resutfmsn_t1$alpha1,resutfmsn_t1$beta1),
          c(resutfmsn_t2$mu,0,resutfmsn_t2$sigma2,0,resutfmsn_t2$shape,0,resutfmsn_t2$nu,0,resutfmsn_st2$pii,0,resutfmsn_t2$w,resutfmsn_t2$mu0,resutfmsn_t2$alpha1,resutfmsn_t2$beta1),
          c(resutfmsn_t3$mu,resutfmsn_t3$sigma2,resutfmsn_t3$shape,resutfmsn_t3$nu,0,resutfmsn_t3$pii,resutfmsn_t3$w,resutfmsn_t3$mu0,resutfmsn_t3$alpha1,resutfmsn_t3$beta1))

write.csv(PAR,"Program/realdata_parameterestimation.csv", sep = ",", col.names = NA,qmethod = "double")


#=================================Criteria============================

res_sn1=goodsn(resutfmsn_sn1$ER,1, resutfmsn_sn1$pii, resutfmsn_sn1$mu, resutfmsn_sn1$sigma2, resutfmsn_sn1$shape, resutfmsn_sn1$nu)
res_sn2=goodsn(resutfmsn_sn2$ER,2, resutfmsn_sn2$pii, resutfmsn_sn2$mu, resutfmsn_sn2$sigma2, resutfmsn_sn2$shape, resutfmsn_sn2$nu)
res_sn3=goodsn(resutfmsn_sn3$ER,3, resutfmsn_sn3$pii, resutfmsn_sn3$mu, resutfmsn_sn3$sigma2, resutfmsn_sn3$shape, resutfmsn_sn3$nu)

res_st1=goodst(resutfmsn_st1$ER,1, resutfmsn_st1$pii, resutfmsn_st1$mu, resutfmsn_st1$sigma2, resutfmsn_st1$shape, resutfmsn_st1$nu)
res_st2=goodst(resutfmsn_st2$ER,2, resutfmsn_st2$pii, resutfmsn_st2$mu, resutfmsn_st2$sigma2, resutfmsn_st2$shape, resutfmsn_st2$nu)
res_st3=goodst(resutfmsn_st3$ER,3, resutfmsn_st3$pii, resutfmsn_st3$mu, resutfmsn_st3$sigma2, resutfmsn_st3$shape, resutfmsn_st3$nu)


res_ss1=goodss(resutfmsn_ss1$ER,1, resutfmsn_ss1$pii, resutfmsn_ss1$mu, resutfmsn_ss1$sigma2, resutfmsn_ss1$shape, resutfmsn_ss1$nu)
res_ss2=goodss(resutfmsn_ss2$ER,2, resutfmsn_ss2$pii, resutfmsn_ss2$mu, resutfmsn_ss2$sigma2, resutfmsn_ss2$shape, resutfmsn_ss2$nu)
res_ss3=goodss(resutfmsn_ss3$ER,3, resutfmsn_ss3$pii, resutfmsn_ss3$mu, resutfmsn_ss3$sigma2, resutfmsn_ss3$shape, resutfmsn_ss3$nu)

res_scn1=goodscn(resutfmsn_scn1$ER,1, resutfmsn_scn1$pii, resutfmsn_scn1$mu, resutfmsn_scn1$sigma2, resutfmsn_scn1$shape, resutfmsn_scn1$nu)
res_scn2=goodscn(resutfmsn_scn2$ER,2, resutfmsn_scn2$pii, resutfmsn_scn2$mu, resutfmsn_scn2$sigma2, resutfmsn_scn2$shape, resutfmsn_scn2$nu)
res_scn3=goodscn(resutfmsn_scn3$ER,3, resutfmsn_scn3$pii, resutfmsn_scn3$mu, resutfmsn_scn3$sigma2, resutfmsn_scn3$shape, resutfmsn_scn3$nu)

res_n1=goodsn(resutfmsn_n1$ER,1, resutfmsn_n1$pii, resutfmsn_n1$mu, resutfmsn_n1$sigma2, resutfmsn_n1$shape, resutfmsn_n1$nu)
res_n2=goodsn(resutfmsn_n2$ER,2, resutfmsn_n2$pii, resutfmsn_n2$mu, resutfmsn_n2$sigma2, resutfmsn_n2$shape, resutfmsn_n2$nu)
res_n3=goodsn(resutfmsn_n3$ER,3, resutfmsn_n3$pii, resutfmsn_n3$mu, resutfmsn_n3$sigma2, resutfmsn_n3$shape, resutfmsn_n3$nu)

res_t1=goodst(resutfmsn_t1$ER,1, resutfmsn_t1$pii, resutfmsn_t1$mu, resutfmsn_t1$sigma2, resutfmsn_t1$shape, resutfmsn_t1$nu)
res_t2=goodst(resutfmsn_t2$ER,2, resutfmsn_t2$pii, resutfmsn_t2$mu, resutfmsn_t2$sigma2, resutfmsn_t2$shape, resutfmsn_t2$nu)
res_t3=goodst(resutfmsn_t3$ER,3, resutfmsn_t3$pii, resutfmsn_t3$mu, resutfmsn_t3$sigma2, resutfmsn_t3$shape, resutfmsn_t3$nu)

#If the parameter estimation information is stored in a file.
PAR=as.matrix(read.csv("Program/realdata_parameter_estimation.csv"),ncol=18)[,-1]

a1=1;a2=2;a3=3
res_sn1=goodsn(ERROR(PAR[a1,15],PAR[a1,16],PAR[a1,17],PAR[a1,18]),1, PAR[a1,12], PAR[a1,1], PAR[a1,4], PAR[a1,7], PAR[a1,10])
res_sn2=goodsn(ERROR(PAR[a2,15],PAR[a2,16],PAR[a2,17],PAR[a2,18]),2, PAR[a2,12:13], PAR[a2,1:2], PAR[a2,4:5], PAR[a2,7:8], PAR[a2,10])
res_sn3=goodsn(ERROR(PAR[a3,15],PAR[a3,16],PAR[a3,17],PAR[a3,18]),3, PAR[a3,12:14], PAR[a3,1:3], PAR[a3,4:6], PAR[a3,7:9], PAR[a3,10])

a1=4;a2=5;a3=6
res_st1=goodst(ERROR(PAR[a1,15],PAR[a1,16],PAR[a1,17],PAR[a1,18]),1, PAR[a1,12], PAR[a1,1], PAR[a1,4], PAR[a1,7], PAR[a1,10])
res_st2=goodst(ERROR(PAR[a2,15],PAR[a2,16],PAR[a2,17],PAR[a2,18]),2, PAR[a2,12:13], PAR[a2,1:2], PAR[a2,4:5], PAR[a2,7:8], PAR[a2,10])
res_st3=goodst(ERROR(PAR[a3,15],PAR[a3,16],PAR[a3,17],PAR[a3,18]),3, PAR[a3,12:14], PAR[a3,1:3], PAR[a3,4:6], PAR[a3,7:9], PAR[a3,10])

a1=7;a2=8;a3=9
res_ss1=goodss(ERROR(PAR[a1,15],PAR[a1,16],PAR[a1,17],PAR[a1,18]),1, PAR[a1,12], PAR[a1,1], PAR[a1,4], PAR[a1,7], PAR[a1,10])
res_ss2=goodss(ERROR(PAR[a2,15],PAR[a2,16],PAR[a2,17],PAR[a2,18]),2, PAR[a2,12:13], PAR[a2,1:2], PAR[a2,4:5], PAR[a2,7:8], PAR[a2,10])
res_ss3=goodss(ERROR(PAR[a3,15],PAR[a3,16],PAR[a3,17],PAR[a3,18]),3, PAR[a3,12:14], PAR[a3,1:3], PAR[a3,4:6], PAR[a3,7:9], PAR[a3,10])

a1=10;a2=11;a3=12
res_scn1=goodscn(ERROR(PAR[a1,15],PAR[a1,16],PAR[a1,17],PAR[a1,18]),1, PAR[a1,12], PAR[a1,1], PAR[a1,4], PAR[a1,7], PAR[a1,10:11])
res_scn2=goodscn(ERROR(PAR[a2,15],PAR[a2,16],PAR[a2,17],PAR[a2,18]),2, PAR[a2,12:13], PAR[a2,1:2], PAR[a2,4:5], PAR[a2,7:8], PAR[a2,10:11])
res_scn3=goodscn(ERROR(PAR[a3,15],PAR[a3,16],PAR[a3,17],PAR[a3,18]),3, PAR[a3,12:14], PAR[a3,1:3], PAR[a3,4:6], PAR[a3,7:9], PAR[a3,10:11])

a1=13;a2=14;a3=15
res_n1=goodn(ERROR(PAR[a1,15],PAR[a1,16],PAR[a1,17],PAR[a1,18]),1, PAR[a1,12], PAR[a1,1], PAR[a1,4], PAR[a1,7], PAR[a1,10])
res_n2=goodn(ERROR(PAR[a2,15],PAR[a2,16],PAR[a2,17],PAR[a2,18]),2, PAR[a2,12:13], PAR[a2,1:2], PAR[a2,4:5], PAR[a2,7:8], PAR[a2,10])
res_n3=goodn(ERROR(PAR[a3,15],PAR[a3,16],PAR[a3,17],PAR[a3,18]),3, PAR[a3,12:14], PAR[a3,1:3], PAR[a3,4:6], PAR[a3,7:9], PAR[a3,10])

a1=16;a2=17;a3=18
res_t1=goodt(ERROR(PAR[a1,15],PAR[a1,16],PAR[a1,17],PAR[a1,18]),1, PAR[a1,12], PAR[a1,1], PAR[a1,4], PAR[a1,7], PAR[a1,10])
res_t2=goodt(ERROR(PAR[a2,15],PAR[a2,16],PAR[a2,17],PAR[a2,18]),2, PAR[a2,12:13], PAR[a2,1:2], PAR[a2,4:5], PAR[a2,7:8], PAR[a2,10])
res_t3=goodt(ERROR(PAR[a3,15],PAR[a3,16],PAR[a3,17],PAR[a3,18]),3, PAR[a3,12:14], PAR[a3,1:3], PAR[a3,4:6], PAR[a3,7:9], PAR[a3,10])

SS=rbind(res_sn1,res_sn2,res_sn3,res_st1,res_st2,res_st3,res_ss1,res_ss2,res_ss3,res_scn1,res_scn2,res_scn3,res_n1,res_n2,res_n3,res_t1,res_t2,res_t3)

write.csv(SS,"Program/Criteria.csv", sep = ",", col.names = NA,qmethod = "double")

#===========================Chose best model========================
which.max(SS[,1])
which.min(SS[,2])

#=======================residuals plot=================
zz <- read_excel("Program/pre1st2$zt.xlsx")
zz=zz[[2]]

X11()
hist(zz, breaks = 100, probability = T,  main = "ST_2")
xx = seq(min(zz), max(zz), (max(zz) - min(zz))/1000)
lines(xx, d.mixedST(xx, c(0.409864,0.590136),c(-0.4672731,.8179275),c(0.6719605,0.9749524),c(1.302494,-1.296279),6.14699))


X11()
qqnorm(zz)
qqline(zz)

X11()
yyy=quantile(zz, probs = seq(0, 1, 0.01), na.rm = FALSE,names = TRUE, type = 7, digits = 7)
library(sn)
zzz1=rst(n=length(yyy), xi=-.467, omega=.672, alpha=1.302, nu=3.747, dp=NULL)
zzz2=rst(n=length(yyy), xi=.818, omega=.975, alpha=-1.296, nu=3.747, dp=NULL)
xxx=0.41*zzz1+.59*zzz2
qqplot(xxx, yyy,main = expression("Q-Q plot for FM-ST_2"),xlab="Theoretical Quantiles",ylab="Sample Quantiles")
qqline(y <- rchisq(length(yyy), df = 3.7), distribution = function(p) qchisq(p, df = 3),probs = c(0.1, 0.6), col = "black")


proc.time()-ptm
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken  




