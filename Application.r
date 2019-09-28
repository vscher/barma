#Application for energy stored rates data
#Version: 1.2 
#Author: Vinícius T. Scher (vinitscher@gmail.com)
#Beta autoregressive moving average models
#Last update: august, 2019
install.packages("gdata")
#Required packages
library(gdata) #allows to read excel files.  
library(e1071) #allows to calculate asymetry and kustosis.
library(forecast) #allows to fit ARMA models.

#Required functions
source("simu_barma.r") #function to generates betaARMA observations.
source("barma.r")      #function to estimate parameters.
source("barma.fit.r")  #function to fit the betaARMA model.
source("karma.r")      #function to estimate KARMA model.
source("kum-mu-phi.r") #function for density function and cumulative distribution function for KARMA model. 
source("karma.fit.r")  #function to fit KARMA model.
source("ljung_box.r")  #function to compute Ljung-Box test.
source("dufour.r")     #function to compute Dufor and Roy test.
source("monti.r")      #function to compute Monti test.
source("Kwan_Chest.R") #function to compute Kwan-Sim tests.

#Program start
#Data energy, 196 months, for jan/2001 to oct/2016.
#read the data energy. (V1 = month, V2 = stored energy rate in %).
data = read.xls("hydrologydata.xls",head=FALSE) 

#enable to accesses variable giving their names.
attach(data)  

#trasform variable V2 (stored energy rate) for unit interval.
y1<-V2/100 

#convert data in time series object.
Y<-ts(y1,start = c(2001,1),frequency = 12)

#sample size
n<-length(Y)
#number of forecast steps.
h1<-6 

#taking off the last 6 observations.
n<-n-h1
y<-ts(y1[1:n],start = c(2001,1),frequency = 12)

#Data description.
summary(y)  #resume
var(y)      #variance.
skewness(y) #skewness.
kurtosis(y) #kurtosis.

#Some graphics, part I
#autocorrelation function.
Acf(y) 

#partial autocorrelation function.
Pacf(y)

#Fit betaARMA model
fit1<-barma(y,ar=c(1),ma=c(1),h=h1,diag=0,resid=1)                       #performing parameter estimation and forecast values.
plot(y,type="l",ylab="energy stored rates",xlab="years")                 #data energy graphics.
lines(fit1$fitted,lty=2,col="red")                                       #predicted values.
legend(x=2000.5,y=0.5,c("observed","predicted"),cex=1, lty=1:2,bty="n")  #legend, observed and predicted.
res<-fit1$resid1                                                         #save residuals.
fc1<-fit1$forecast                                                       #save the forecast values.

#Portmanteau tests
#lag m vector
vm<-c(seq(from = 3, to = 30, by = 1))
f<-2 #degrees of freedom

#initializing response vector
QLB<-matrix(rep(NA,30),nrow=1, ncol=30) 
QMO<-matrix(rep(NA,30),nrow=1, ncol=30) 
QDR<-matrix(rep(NA,30),nrow=1, ncol=30) 
QKW1<-matrix(rep(NA,30),nrow=1, ncol=30) 
QKW4<-matrix(rep(NA,30),nrow=1, ncol=30) 
Q4<-matrix(rep(NA,30),nrow=1, ncol=30) 

#getting p-values
for(m in vm)
{
  QLB[,m]<-Ljung.Box(res, lag=m, fitdf =f)$p.value                                  #Ljung-Box test
  QMO[,m]<-Monti.test(res, lag=m, fitdf = f)$p.value                                #Monti test
  QDR[,m]<-Dufour.test(res, lag=m, fitdf = f)$p.value                               #Dufor and Roy test 
  QKW1[,m]<-Kwan.sim.chest(res, lag=m, fitdf = f,test=1)$p.value                    #Kwan and Sim test(1) 
  QKW4[,m]<-Kwan.sim.chest(res, lag=m, fitdf = f,type="correlation",test=4)$p.value #Kwan and Sim test(4)  
  Q4[,m]<-Kwan.sim.chest(res, lag=m, fitdf = f,type="partial",test=4)$p.value       #New portmanteau test statistic  
}

#graphics for portmanteau tests
nlag <- 30

#Ljung-Box
plot(1L:nlag, QLB, xlab = "lag",ylab = "p-value",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#Monti
plot(1L:nlag, QMO, xlab = "lag",ylab = "p-values",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#Dufour and Roy
plot(1L:nlag, QDR, xlab = "lag",ylab = "p-values",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#Kwan-Sim(1)
plot(1L:nlag, QKW1, xlab = "lag",ylab = "p-values",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#Kwan-Sim(4)
plot(1L:nlag, QKW4, xlab = "lag",ylab = "p-values",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#New portmanteau test statistic 
plot(1L:nlag, Q4, xlab = "lag",ylab = "p-values",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)


#Misspecification when applied to a diferent model.
#beta-AR(1)
res2<-barma(y,ar=c(1),ma=NA,h=0,diag=0,resid=1)$resid1
f<-1
#getting p-values
for(m in vm)
{
  QLB[,m]<-Ljung.Box(res2, lag=m, fitdf =f)$p.value                                  #Ljung-Box test
  QMO[,m]<-Monti.test(res2, lag=m, fitdf = f)$p.value                                #Monti test
  QDR[,m]<-Dufour.test(res2, lag=m, fitdf = f)$p.value                               #Dufor and Roy test 
  QKW1[,m]<-Kwan.sim.chest(res2, lag=m, fitdf = f,test=1)$p.value                    #Kwan and Sim test(1) 
  QKW4[,m]<-Kwan.sim.chest(res2, lag=m, fitdf = f,type="correlation",test=4)$p.value #Kwan and Sim test(4)  
  Q4[,m]<-Kwan.sim.chest(res2, lag=m, fitdf = f,type="partial",test=4)$p.value       #New portmanteau test statistic  
}

#initializing the result matrix for p-values AR(1) is fitted
table4<-matrix(rep(NA,36),nrow=6,ncol=6) #declaring result matrix

for(i in 1:6)
{
  table4[1,i]<-round(QLB[,i*5],6)
  table4[2,i]<-round(QMO[,i*5],6)
  table4[3,i]<-round(QDR[,i*5],6)
  table4[4,i]<-round(QKW1[,i*5],6)
  table4[5,i]<-round(QKW4[,i*5],6)
  table4[6,i]<-round(Q4[,i*5],6)
}
table4

#Ljung-Box
plot(1L:nlag, QLB, xlab = "lag",ylab = "p-value",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#Monti
plot(1L:nlag, QMO, xlab = "lag",ylab = "p-values",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#Dufour and Roy
plot(1L:nlag, QDR, xlab = "lag",ylab = "p-values",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#Kwan-Sim(1)
plot(1L:nlag, QKW1, xlab = "lag",ylab = "p-values",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#Kwan-Sim(4)
plot(1L:nlag, QKW4, xlab = "lag",ylab = "p-values",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#New portmanteau test statistic 
plot(1L:nlag, Q4, xlab = "lag",ylab = "p-values",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#beta-MA(1)
res3<-barma(y,ar=NA,ma=c(1),h=0,diag=0,resid=1)$resid1
f<-1
#getting p-values
for(m in vm)
{
  QLB[,m]<-Ljung.Box(res3, lag=m, fitdf =f)$p.value                                  #Ljung-Box test
  QMO[,m]<-Monti.test(res3, lag=m, fitdf = f)$p.value                                #Monti test
  QDR[,m]<-Dufour.test(res3, lag=m, fitdf = f)$p.value                               #Dufor and Roy test 
  QKW1[,m]<-Kwan.sim.chest(res3, lag=m, fitdf = f,test=1)$p.value                    #Kwan and Sim test(1) 
  QKW4[,m]<-Kwan.sim.chest(res3, lag=m, fitdf = f,type="correlation",test=4)$p.value #Kwan and Sim test(4)  
  Q4[,m]<-Kwan.sim.chest(res3, lag=m, fitdf = f,type="partial",test=4)$p.value       #New portmanteau test statistic  
}

#Ljung-Box
plot(1L:nlag, QLB, xlab = "lag",ylab = "p-value",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#Monti
plot(1L:nlag, QMO, xlab = "lag",ylab = "p-values",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#Dufour and Roy
plot(1L:nlag, QDR, xlab = "lag",ylab = "p-values",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#Kwan-Sim(1)
plot(1L:nlag, QKW1, xlab = "lag",ylab = "p-values",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#Kwan-Sim(4)
plot(1L:nlag, QKW4, xlab = "lag",ylab = "p-values",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)

#New portmanteau test statistic 
plot(1L:nlag, Q4, xlab = "lag",ylab = "p-values",las=1, ylim = c(0,1),lwd=1.5)
abline(h = 0.05, lty = 2, col = "blue",lwd=1.5)


#fit usual ARMA models using functions of forecast package
auto.arima(y)
y_arma<-arima(y,order = c(1,0,1))                                  #performing parameter estimation for ARMA(1,1).
fc2<-forecast(y_arma,h=h1)$mean[1:6]                               #performing forecast values.
y_ar<-arima(y,order = c(2,0,0))                                    #performing paramter estimation for AR(2).
fc3<-forecast(y_ar,h=h1)$mean[1:6]                                 #performing forecast values.

#fit KARMA models. 
fc4<-karma(y,ar=c(1),ma=c(1),h=h1,diag=0,X = NA,X_hat=NA)$forecast #selecting values to forecast.

#forecast with Holt-Winters algoritm.
fc5<-holt(y,h=h1)$mean[1:6]                                        #using hw function without seasonal component for forecast package.
                 
#Mean Absolute error(MAE)
#initializing the result matrix
table5<-matrix(rep(NA,30),nrow=5,ncol=6)                           #declaring result matrix

#save the observations to evaluate forecasting accuracy.
obs<-Y[191:196]

#loop for results.
for(i in 1:6)
{
  table5[1,i]<-round(mean(abs((obs[1:i] - fc1[1:i]))),6)
  table5[2,i]<-round(mean(abs((obs[1:i] - fc2[1:i]))),6)
  table5[3,i]<-round(mean(abs((obs[1:i] - fc3[1:i]))),6)
  table5[4,i]<-round(mean(abs((obs[1:i] - fc4[1:i]))),6)
  table5[5,i]<-round(mean(abs((obs[1:i] - fc5[1:i]))),6)
}

table5


