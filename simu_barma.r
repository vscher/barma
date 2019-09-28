# n: tamanho da amostra
# freq: frequencia anual de observacoes (12 para observacoes mensais)

simu.barma <- function(n,phi=NA,theta=NA,
                       alpha=0.0,prec=100,freq=12,link="logit")
{
  ar<-NA
  ma<-NA
  
  if(any(is.na(phi)==F))
  {
    ar <- 1:length(phi)
  }
  
  if(any(is.na(theta)==F))
  {
    ma <- 1:length(theta)
  }
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))
  {  
    stats <- make.link(linktemp)
  }else{
    stop(paste(linktemp, "link not available, available links are \"logit\", ",
               "\"probit\" and \"cloglog\""))
  } 
  
  link <- structure(list(link = linktemp, 
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv
  )
  )
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  
  
  ###### ARMA model
  if(any(is.na(phi)==F) && any(is.na(theta)==F))
  {
    #print("ARMA model")
    # seasonal part
    
    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)
    
    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)
    
    error<-rep(0,n+m) # E(error)=0 
    eta<- y <- rep(NA,n+m)
    
    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + as.numeric(phi%*%ynew[i-ar]) + as.numeric(theta%*%error[i-ma])
      mu[i]   <- linkinv(eta[i])
      y[i]    <- rbeta(1, mu[i]*prec, (1-mu[i])*prec)
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-eta[i]   
      #error[i]<- y[i]-mu[i]
    }
    
    return(ts(y[(m+1):(n+m)],frequency=freq) )
    
  } # ARMA model
  
  
  ###### AR model
  if(any(is.na(phi)==F) && any(is.na(theta)==T))
  {
    #print("AR model")    
    
    p <- max(ar)
    m <- 2*p
    
    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)
    
    eta<- y <- rep(NA,n+m)
    
    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + (phi%*%ynew[i-ar]) 
      mu[i]   <- linkinv(eta[i])
      y[i]    <- rbeta(1, mu[i]*prec, (1-mu[i])*prec)
      ynew[i] <- linkfun(y[i])
    }
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } # AR model
  
  
  ###### MA model
  if(any(is.na(phi)==T) && any(is.na(theta)==F))
  {
    #print("MA model")    
    
    q <- max(ma)
    m <- 2*q
    
    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)
    
    eta <- y <- error <- rep(0,n+m) # E(error)=0 
    
    #print(ma)
    
    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + (theta%*%error[i-ma])
      mu[i]   <- linkinv(eta[i])
      y[i]    <- rbeta(1, mu[i]*prec, (1-mu[i])*prec)
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-eta[i]   
      #error[i]<- y[i]-mu[i]
    }
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } # fim MA model  
  
}

# y<-simu.barma(50,phi=c(0.45,0.3),theta=c(0.5),prec=100)#,AR="NA",MA="NA",alpha=0.0,prec=100,freq=12,link="logit")
# par(mfrow=c(1,1))
# plot(y)
# print(summary(y) )
# print(length(y) )
# 
# fit1<-barma(y,ar=c(1,2),ma=c(1))#,AR=c(1,2),MA=c(1,2))#, link = "probit")
# print(fit1$coef)
# 
# simu.barma(50,phi=c(0.4,-0.3))
# simu.barma(50,phi=c(0.4,-0.3),theta=c(-0.2))
# simu.barma(50,phi=c(0.4,-0.3),PHI=c(0.3))
# simu.barma(50,phi=c(0.4,-0.3),THETA=c(0.3))
# simu.barma(50,phi=c(0.4,-0.3),THETA=c(0.3),PHI=c(-0.2))
# simu.barma(50,phi=c(0.4,-0.3),theta=c(-0.2),THETA=c(0.3),PHI=c(-0.2))
# 
# simu.barma(50,theta=c(0.4,-0.3))
# simu.barma(50,theta=c(0.4,-0.3),PHI=c(0.3))
# simu.barma(50,theta=c(0.4,-0.3),THETA=c(0.3))
# simu.barma(50,theta=c(0.4,-0.3),THETA=c(0.3),PHI=c(-0.2))

