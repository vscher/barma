# Created by Fabio M Bayer (bayer@ufsm.br), july/2017

karma.fit<- function (y, ar, ma, a, b, link, names_phi,names_theta,names_beta,diag,h1,X,X_hat,resid)
{
  maxit1<-100
  
  z <- c()
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta <-  link$mu.eta
  diflink <- link$diflink
  
  ynew = linkfun(y)
  ystar = log(y/(1-y))
  
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)
  
  y_prev <- c(rep(NA,(n+h1)))
  
  
  # initial values
  if(any(is.na(ar)==F)) # with AR
  {
    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
    
    for(i in 1:(n-m))
    {
      P[i,] <- ynew[i+m-ar]
    }
    
    Z <- cbind(rep(1,(n-m)),P)
  }else{
    Z <- as.matrix(rep(1,(n-m)))
  }
  
  if(any(is.na(X)==T)) # without regressors
  {
    x <- as.matrix(Z)
    Y <- y[(m+1):n]
    Ynew = linkfun(Y)
    ajuste = lm.fit(x, Ynew)
    mqo = c(ajuste$coef)
    k = length(mqo)
    n1 = length(Y)
    mean = fitted(ajuste)
    mean = linkinv(mean)
    
    prec<- 10
    
  }else{ # com regressores
    X_hat<-as.matrix(X_hat)
    X<-as.matrix(X)
    x <- cbind(as.matrix(Z),X[(m+1):n,])
    Y <- y[(m+1):n]
    Ynew = linkfun(Y)
    Ystar = log(Y/(1-Y))
    ajuste = lm.fit(x, Ynew)
    mqo = c(ajuste$coef)
    k = length(mqo)
    n1 = length(Y)
    mean = fitted(ajuste)
    mean = linkinv(mean)
    dlink = diflink(mean)
    er = residuals(ajuste)
    sigma2 = sum(er^2)/((n1 - k) * (dlink)^2)
    
    prec<- 10

  }
  
  
  ############
  ######### ARMA model
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==T))
  { 
    #print("KARMA model",quote=F)
    reg <- c(mqo, rep(0,q1), prec) # initializing the parameter values
    
    loglik <- function(z) 
    {
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      theta = z[(p1+2):(p1+q1+1)]
      prec <- z[p1+q1+2] # precision parameter
      
      error<-rep(0,n) # E(error)=0 
      eta<-rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + (phi%*%ynew[i-ar]) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i] 
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings( log( dkum(y1,mu,prec) ) )
      sum(ll)
    } 
    
    escore <- function(z) 
    {
      alpha <- z[1]
      phi <- z[2:(p1+1)]
      theta <- z[(p1+2):(p1+q1+1)]
      prec <- z[p1+q1+2] # precision parameter
      
      error<- rep(0,n) # E(error)=0 
      eta<- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i]<- alpha + (phi%*%ynew[i-ar]) + (theta%*%error[i-ma])
        error[i]<- ynew[i]-eta[i] 
      }
      
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m))
      {
        R[i,] <- error[i+m-ma]
      }
      
      ###FB  recorrences
      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      
      for(i in (m+1):n)
      {
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      }
      
      v <- deta.dalpha[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rR <- deta.dtheta[(m+1):n,]
      
      mT <- diag(mu.eta(eta[(m+1):n]))
      
      delta <- log(0.5)/log(1-mu^prec)
      c1 <- mu^(prec-1)/( (1-mu^prec)*log(1-mu^prec) )
      c2 <- delta*log(1-y1^prec)+1
      vc <- prec*c1*c2
      
      a1 <- 1/prec
      a2 <- log(y1)
      a3 <- vc*mu*log(mu)/prec
      a4 <- delta-1
      a5 <- ((y1^prec)*log(y1))/(1-y1^prec)
      a <- as.vector(a1+a2+a3-(a4*a5)) # derivada em relação a phi
      
      #vc <- as.vector((mu^(prec-1))/((1-mu^prec)*log(1-mu^prec)) * ((log(0.5)*log(1-y1^prec) )/log(1-mu^prec) +1))
      
      Ualpha <- t(v) %*% mT %*% vc
      Uphi <-   t(rP) %*% mT %*% vc
      Utheta <- t(rR) %*% mT %*% vc
      #Uprec <-   sum(1/prec + log(y1) + ((mu^prec)*log(mu))/( (1-mu^prec)*log(1-mu^prec) ) * ((log(0.5)*log(1-y1^prec) )/(log(1-mu^prec))  + 1) 
      #               - (y1^prec * log(y1))/(1-y1^prec) * ( log(0.5)/log(1-mu^prec) -1))
      Uprec <- sum(a)               
      
      rval <- c(Ualpha,Uphi,Utheta,Uprec)
    }
    names_par <- c("alpha",names_phi,names_theta,"precision")
    
    opt <- optim(reg, loglik, escore, 
                 method = "BFGS", 
                 control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      opt <- optim(reg, loglik, #escore, 
                   method = "BFGS", 
                   control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
    
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+q1+2)]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    prec <- coef[p1+q1+2] # precision parameter
    
    z$alpha <- alpha
    z$phi <- phi
    z$theta <- theta
    z$prec <- prec
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + (phi%*%ynew[i-ar]) + (theta%*%errorhat[i-ma])
      errorhat[i]<- ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      R[i,] <- errorhat[i+m-ma]
    }
    
    ###FB  recorrences
    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    
    for(i in (m+1):n)
    {
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
    }
    
    v <- matrix(as.vector(deta.dalpha[(m+1):n]),ncol=1)
    rP <- deta.dphi[(m+1):n,]
    rR <- deta.dtheta[(m+1):n,]
    
    mT <- diag(mu.eta(etahat[(m+1):n]))
    c <- as.vector((muhat^(prec-1))/((1-muhat^prec)*log(1-muhat^prec)) * ((log(0.5)*log(1-y1^prec) )/log(1-muhat^prec) +1))
    vc <- prec*c
    
    lambda1 <- (muhat^(prec-2))/((1-muhat^prec) * (log(1-muhat^prec)))
    
    lambda2 <- (muhat^(2*prec-2))/((1-muhat^prec)^2 * (log(1-muhat^prec))^2)
    
    W <- diag(as.vector(- prec^2 * lambda2))
    vI <- matrix(rep(1,(n-m)),ncol=1)
    
    euler <- 0.5772156649
    delta <- log(0.5)/(log(1-muhat^prec))
    d <- -prec*muhat*log(muhat)*lambda2-prec*muhat*delta*lambda1*((1-digamma(delta+1)-euler)/((delta-1)*prec))
    D <- diag(as.vector(d))
    L1 <- 1/(prec^2) 
    L2<- delta*(muhat^2)*lambda2*(log(muhat)^2)*log(1-y1^prec)
    L3<- 2*delta*(muhat^2)*log(muhat)*lambda1*(( (y1^prec) * log(y1))/(1-y1^prec))
    L4<- (delta-1)*((y1^prec)*log(y1)^2)/((1-y1^prec)^2)   
    L5<- c*muhat*(log(muhat)^2)*(lambda1*(muhat^2)+1/(1-muhat^prec))
    
    L<- -L1+L2-L3-L4+L5
    
    Kaa <- t(v) %*% W %*% mT^2 %*% v
    Kap <- t(v) %*% W %*% mT^2 %*% rP
    Kpa <- t(Kap)
    Kat <- t(v) %*% W %*% mT^2 %*% rR
    Kta <- t(Kat)
    Kaprec <- t(v) %*% D %*% mT %*% vI
    Kpreca <- t(Kaprec)
    
    Kpp <- t(rP) %*% W %*% mT^2 %*% rP
    Kpt <- t(rP) %*% W %*% mT^2 %*% rR
    Ktp <- t(Kpt)
    Kpprec <- t(rP) %*% D %*% mT %*% vI 
    Kprecp <- t(Kpprec)
    
    Ktt <- t(rR) %*% W %*% mT^2 %*% rR
    Ktprec <- t(rR) %*% D %*% mT %*% vI
    Kprect <- t(Ktprec)
    
    Kprecprec <- sum(L)
    
    K <- -rbind(
      cbind(Kaa,Kap,Kat,Kaprec),
      cbind(Kpa,Kpp,Kpt,Kpprec),
      cbind(Kta,Ktp,Ktt,Ktprec),
      cbind(Kpreca,Kprecp,Kprect,Kprecprec)
    )
    
    z$K <- K
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar]) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 
    }
  }  
  
  
  #   ##### only AR model
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==T))
  {
    #print("KAR model",quote=F)
    q1<-0
    reg <- c(mqo, prec) # initializing the parameter values
    
    loglik <- function(z)
    {
      alpha <- z[1]
      phi = z[2:(p1+1)]
      prec <- z[p1+2] # precision parameter
      
      eta<-rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i]<-alpha + (phi%*%ynew[i-ar])
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings( log( dkum(y1,mu,prec) ) )
      sum(ll)
    }
    
    escore <- function(z)
    {
      alpha <- z[1]
      phi <- z[2:(p1+1)]
      prec <- z[p1+2] # precision parameter
      
      eta<- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i]<- alpha + (phi%*%ynew[i-ar])
      }
      
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      vI <- as.vector(rep(1,n-m))

      mT <- diag(mu.eta(eta[(m+1):n]))
      
      delta <- log(0.5)/log(1-mu^prec)
      c1 <- mu^(prec-1)/( (1-mu^prec)*log(1-mu^prec) )
      c2 <- delta*log(1-y1^prec)+1
      vc <- prec*c1*c2
      
      a1 <- 1/prec
      a2 <- log(y1)
      a3 <- vc*mu*log(mu)/prec
      a4 <- delta-1
      a5 <- ((y1^prec)*log(y1))/(1-y1^prec)
      a <- as.vector(a1+a2+a3-(a4*a5)) # derivada em relação a phi
      
      Ualpha <- t(vI) %*% mT %*% vc
      Uphi <-   t(P) %*% mT %*% vc
      Uprec <-   sum(a)
      
      rval <- c(Ualpha,Uphi,Uprec)
      
    }
    names_par <- c("alpha",names_phi,"precision")
    
    opt <- optim(reg, loglik, escore, 
                 method = "BFGS", 
                 control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      opt <- optim(reg, loglik, #escore, 
                   method = "BFGS", 
                   control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+2)]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    prec <- coef[p1+2] # precision parameter
    
    z$alpha <- alpha
    z$phi <- phi
    z$prec <- prec
    
    errorhat<-rep(0,n) # E(error)=0
    etahat<-rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + (phi%*%ynew[i-ar])
      errorhat[i]<-ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    v <- matrix(rep(1,(n-m)),ncol=1)
    rP <- P
    
    mT <- diag(mu.eta(etahat[(m+1):n]))
    c <- as.vector((muhat^(prec-1))/((1-muhat^prec)*log(1-muhat^prec)) * ((log(0.5)*log(1-y1^prec) )/log(1-muhat^prec) +1))
    vc <- prec*c
    lambda2 <- (muhat^(2*prec-2))/((1-muhat^prec)^2 * (log(1-muhat^prec))^2)
    W <- diag(as.vector(- prec^2 * lambda2))
    vI <- matrix(rep(1,(n-m)),ncol=1)
    lambda1 <- (muhat^(prec-2))/((1-muhat^prec) * (log(1-muhat^prec)))
    euler <- 0.5772156649
    delta <- log(0.5)/(log(1-muhat^prec))
    d <- -prec*muhat*log(muhat)*lambda2-prec*muhat*delta*lambda1*((1-digamma(delta+1)-euler)/((delta-1)*prec))
    D <- diag(as.vector(d))
    L1 <- 1/(prec^2) 
    L2<- delta*(muhat^2)*lambda2*(log(muhat)^2)*log(1-y1^prec)
    L3<- 2*delta*(muhat^2)*log(muhat)*lambda1*(( (y1^prec) * log(y1))/(1-y1^prec))
    L4<- (delta-1)*((y1^prec)*log(y1)^2)/((1-y1^prec)^2)   
    L5<- c*muhat*(log(muhat)^2)*(lambda1*(muhat^2)+1/(1-muhat^prec))
    
    L<- -L1+L2-L3-L4+L5
    
    Kaa <- t(v) %*% W %*% mT^2 %*% v
    Kap <- t(v) %*% W %*% mT^2 %*% rP
    Kpa <- t(Kap)
    Kaprec <- t(v) %*% D %*% mT %*% vI
    Kpreca <- t(Kaprec)
    
    Kpp <- t(rP) %*% W %*% mT^2 %*% rP
    Kpprec <- t(rP) %*% D %*% mT %*% vI 
    Kprecp <- t(Kpprec)
    
    Kprecprec <- sum(L)
    
    K <- -rbind(
      cbind(Kaa,Kap,Kaprec),
      cbind(Kpa,Kpp,Kpprec),
      cbind(Kpreca,Kprecp,Kprecprec)
    )
    
    z$K <- K
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
    }
    
  }
  
  

  #   ######### MA model
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==T))
  {
    p1<-0
    #print("KMA model",quote=F)
    reg <- c(mqo,rep(0,q1), prec) # initializing the parameter values
    
    loglik <- function(z)
    {
      alpha <- z[1]
      theta <- z[2:(q1+1)]
      prec <- z[q1+2] # precision parameter
      
      eta <- error <- rep(0,n) # E(error)=0
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + (theta%*%error[i-ma])
        error[i]<-ynew[i]-eta[i]
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings( log( dkum(y1,mu,prec) ) )
      sum(ll)
    }
    
    escore <- function(z)
    {
      alpha <- z[1]
      theta <- z[2:(q1+1)]
      prec <- z[q1+2] # precision parameter
      
      error<- rep(0,n) # E(error)=0
      eta<- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + (theta%*%error[i-ma])
        error[i]<- ynew[i]-eta[i]
      }
      
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m))
      {
        R[i,] <- error[i+m-ma]
      }
      
      ###FB  recorrences
      deta.dalpha <- rep(0,n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      
      for(i in (m+1):n)
      {
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      }
      
      v <- deta.dalpha[(m+1):n]
      rR <- deta.dtheta[(m+1):n,]
      
      mT <- diag(mu.eta(eta[(m+1):n]))
      
      delta <- log(0.5)/log(1-mu^prec)
      c1 <- mu^(prec-1)/( (1-mu^prec)*log(1-mu^prec) )
      c2 <- delta*log(1-y1^prec)+1
      vc <- prec*c1*c2
      
      a1 <- 1/prec
      a2 <- log(y1)
      a3 <- vc*mu*log(mu)/prec
      a4 <- delta-1
      a5 <- ((y1^prec)*log(y1))/(1-y1^prec)
      a <- as.vector(a1+a2+a3-(a4*a5)) # derivada em relação a phi
      
      Ualpha <- t(v) %*% mT %*% vc
      Utheta <- t(rR) %*% mT %*% vc
      Uprec <-   sum(a)      
      
      rval <- c(Ualpha,Utheta,Uprec)
    }
    names_par <- c("alpha",names_theta,"precision")
    
    opt <- optim(reg, loglik, escore, 
                 method = "BFGS", 
                 control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      opt <- optim(reg, loglik, #escore, 
                   method = "BFGS", 
                   control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(q1+2)]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    theta <- coef[2:(q1+1)]
    prec <- coef[q1+2] # precision parameter
    
    z$alpha <- alpha
    z$theta <- theta
    z$prec <- prec
    
    errorhat<-rep(0,n) # E(error)=0
    etahat<-rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + (theta%*%errorhat[i-ma])
      errorhat[i]<-ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      R[i,] <- errorhat[i+m-ma]
    }
    
    ###FB  recorrences
    deta.dalpha <- rep(0,n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    
    for(i in (m+1):n)
    {
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
    }
    
    v <- deta.dalpha[(m+1):n]
    rR <- deta.dtheta[(m+1):n,]
    
    mT <- diag(mu.eta(etahat[(m+1):n]))
    c <- as.vector((muhat^(prec-1))/((1-muhat^prec)*log(1-muhat^prec)) * ((log(0.5)*log(1-y1^prec) )/log(1-muhat^prec) +1))
    vc <- prec*c
    lambda2 <- (muhat^(2*prec-2))/((1-muhat^prec)^2 * (log(1-muhat^prec))^2)
    W <- diag(as.vector(- prec^2 * lambda2))
    vI <- matrix(rep(1,(n-m)),ncol=1)
    lambda1 <- (muhat^(prec-2))/((1-muhat^prec) * (log(1-muhat^prec)))
    euler <- 0.5772156649
    delta <- log(0.5)/(log(1-muhat^prec))
    d <- -prec*muhat*log(muhat)*lambda2-prec*muhat*delta*lambda1*((1-digamma(delta+1)-euler)/((delta-1)*prec))
    D <- diag(as.vector(d))
    L1 <- 1/(prec^2) 
    L2<- delta*(muhat^2)*lambda2*(log(muhat)^2)*log(1-y1^prec)
    L3<- 2*delta*(muhat^2)*log(muhat)*lambda1*(( (y1^prec) * log(y1))/(1-y1^prec))
    L4<- (delta-1)*((y1^prec)*log(y1)^2)/((1-y1^prec)^2)   
    L5<- c*muhat*(log(muhat)^2)*(lambda1*(muhat^2)+1/(1-muhat^prec))
    
    L<- -L1+L2-L3-L4+L5
    
    Kaa <- t(v) %*% W %*% mT^2 %*% v
    Kat <- t(v) %*% W %*% mT^2 %*% rR
    Kta <- t(Kat)
    Kaprec <- t(v) %*% D %*% mT %*% vI
    Kpreca <- t(Kaprec)
    
    Ktt <- t(rR) %*% W %*% mT^2 %*% rR
    Ktprec <- t(rR) %*% D %*% mT %*% vI
    Kprect <- t(Ktprec)
    
    Kprecprec <- sum(L)
    
    
    K <- -rbind(
      cbind(Kaa,Kat,Kaprec),
      cbind(Kta,Ktt,Ktprec),
      cbind(Kpreca,Kprect,Kprecprec)
    )
    
    z$K <- K
    
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 # original scale
    }
  } # fim KMA
  
  

  ######### KARMAX model
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==F))
  { 
    #print("KARMAX model",quote=F)
    beta1<- mqo[(p1+2):length(mqo)]
    reg <- c(mqo[1:(p1+1)], rep(0,q1), prec, beta1) # initializing the parameter values
    
    loglik <- function(z) 
    {
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      theta = z[(p1+2):(p1+q1+1)]
      prec <- z[p1+q1+2] # precision parameter
      beta <- z[(p1+q1+3):length(z)]
      
      error<-rep(0,n) # E(error)=0 
      eta<-rep(NA,n)
      
      for(i in (m+1):n)
      {
        
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i] # predictor scale
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings( log( dkum(y1,mu,prec) ) )
      sum(ll)
    } 
    
    escore <- function(z) 
    {
      alpha <- z[1]
      phi <- z[2:(p1+1)]
      theta <- z[(p1+2):(p1+q1+1)]
      prec <- z[p1+q1+2] # precision parameter
      beta <- z[(p1+q1+3):length(z)]
      
      error<- rep(0,n) # E(error)=0 
      eta<- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i]<- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%error[i-ma])
        error[i]<- ynew[i]-eta[i] 
      }
      
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m))
      {
        R[i,] <- error[i+m-ma]
      }
      
      P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
      for(i in 1:(n-m))
      {
        P[i,] <- ynew[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
      }
      
      k1<- length(beta)
      M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
      for(i in 1:(n-m))
      {
        # for(j in 1:k1)
        #   M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])  

        M[i,] <- X[i+m,]-(phi%*%X[i+m-ar,]) 
      }
      
      ###FB  recorrences
      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      deta.dbeta<- matrix(0, ncol=k1,nrow=n)
      
      for(i in (m+1):n)
      {
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
        deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
      }
      
      v <- deta.dalpha[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rR <- deta.dtheta[(m+1):n,]
      rM <- deta.dbeta[(m+1):n,]
      
      mT <- diag(mu.eta(eta[(m+1):n]))
      
      delta <- log(0.5)/log(1-mu^prec)
      c1 <- mu^(prec-1)/( (1-mu^prec)*log(1-mu^prec) )
      c2 <- delta*log(1-y1^prec)+1
      vc <- prec*c1*c2
      
      a1 <- 1/prec
      a2 <- log(y1)
      a3 <- vc*mu*log(mu)/prec
      a4 <- delta-1
      a5 <- ((y1^prec)*log(y1))/(1-y1^prec)
      a <- as.vector(a1+a2+a3-(a4*a5)) # derivada em relação a phi
      
      Ualpha <- t(v) %*% mT %*% vc
      Uphi <-   t(rP) %*% mT %*% vc
      Utheta <- t(rR) %*% mT %*% vc
      Uprec <-   sum(a)      
      Ubeta <-  t(rM) %*% mT %*% vc
      
      rval <- c(Ualpha,Uphi,Utheta,Uprec,Ubeta)
    }
    names_par <- c("alpha",names_phi,names_theta,"precision",names_beta)
    
    opt <- optim(reg, loglik, escore, 
                 method = "BFGS", 
                 control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      opt <- optim(reg, loglik, #escore, 
                   method = "BFGS", 
                   control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
   
    
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+q1+2+ncol(X))]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    prec <- coef[p1+q1+2] # precision parameter
    beta <- coef[(p1+q1+3):length(coef)]
    
    
    z$alpha <- alpha
    z$phi <- phi
    z$theta <- theta
    z$prec <- prec
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%errorhat[i-ma])
      errorhat[i]<- ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * prec) - digamma((1 - muhat) * prec)
    
    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      R[i,] <- errorhat[i+m-ma]
    }
    
    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
    for(i in 1:(n-m))
    {
      P[i,] <- ynew[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
    }
    k1<- length(beta)
    M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
    for(i in 1:(n-m))
    {
      for(j in 1:k1)
        M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])  
    }
    
    
    ###FB  recorrences
    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    deta.dbeta<- matrix(0, ncol=k1,nrow=n)
    
    for(i in (m+1):n)
    {
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
    }
    
    v <- deta.dalpha[(m+1):n]
    rP <- deta.dphi[(m+1):n,]
    rR <- deta.dtheta[(m+1):n,]
    rM <- deta.dbeta[(m+1):n,]
    
    
    mT <- diag(mu.eta(etahat[(m+1):n]))
    c <- as.vector((muhat^(prec-1))/((1-muhat^prec)*log(1-muhat^prec)) * ((log(0.5)*log(1-y1^prec) )/log(1-muhat^prec) +1))
    vc <- prec*c
    lambda2 <- (muhat^(2*prec-2))/((1-muhat^prec)^2 * (log(1-muhat^prec))^2)
    W <- diag(as.vector(- prec^2 * lambda2))
    vI <- matrix(rep(1,(n-m)),ncol=1)
    lambda1 <- (muhat^(prec-2))/((1-muhat^prec) * (log(1-muhat^prec)))
    euler <- 0.5772156649
    delta <- log(0.5)/(log(1-muhat^prec))
    d <- -prec*muhat*log(muhat)*lambda2-prec*muhat*delta*lambda1*((1-digamma(delta+1)-euler)/((delta-1)*prec))
    D <- diag(as.vector(d))
    L1 <- 1/(prec^2) 
    L2<- delta*(muhat^2)*lambda2*(log(muhat)^2)*log(1-y1^prec)
    L3<- 2*delta*(muhat^2)*log(muhat)*lambda1*(( (y1^prec) * log(y1))/(1-y1^prec))
    L4<- (delta-1)*((y1^prec)*log(y1)^2)/((1-y1^prec)^2)   
    L5<- c*muhat*(log(muhat)^2)*(lambda1*(muhat^2)+1/(1-muhat^prec))
 
    L<- -L1+L2-L3-L4+L5

    Kaa <- t(v) %*% W %*% mT^2 %*% v
    Kap <- t(v) %*% W %*% mT^2 %*% rP
    Kpa <- t(Kap)
    Kat <- t(v) %*% W %*% mT^2 %*% rR
    Kta <- t(Kat)
    Kaprec <- t(v) %*% D %*% mT %*% vI
    Kpreca <- t(Kaprec)
    
    Kpp <- t(rP) %*% W %*% mT^2 %*% rP
    Kpt <- t(rP) %*% W %*% mT^2 %*% rR
    Ktp <- t(Kpt)
    Kpprec <- t(rP) %*% D %*% mT %*% vI 
    Kprecp <- t(Kpprec)
    
    Ktt <- t(rR) %*% W %*% mT^2 %*% rR
    Ktprec <- t(rR) %*% D %*% mT %*% vI
    Kprect <- t(Ktprec)
    
    Kprecprec <- sum(L)
    
    Kab <- t(v) %*% W %*% mT^2 %*% rM
    Kba <- t(Kab)
    Kbb <- t(rM) %*% W %*% mT^2 %*% rM
    Kpb <- t(rP) %*% W %*% mT^2 %*% rM
    Kbp <- t(Kpb)
    Ktb <- t(rR) %*% W %*% mT^2 %*% rM
    Kbt <- t(Ktb)
    Kbprec <- t(rM) %*% D %*% mT %*% vI
    Kprecb <- t(Kbprec)
    
    K <- -rbind(
      cbind(Kaa,Kap,Kat,Kaprec,Kab),
      cbind(Kpa,Kpp,Kpt,Kpprec,Kpb),
      cbind(Kta,Ktp,Ktt,Ktprec,Ktb),
      cbind(Kpreca,Kprecp,Kprect,Kprecprec,Kprecb),
      cbind(Kba,Kbp,Kbt,Kbprec,Kbb)
    )
    
    z$K <- K
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    X_prev<- rbind(X,X_hat)
    
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta) )) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 
    }
  }  

  
  
  ############# KARX model
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==F))
  { 
    q1<-0
    #print("KARX model",quote=F)
    beta1<- mqo[(p1+2):length(mqo)]
    reg <- c(mqo[1:(p1+1)], prec, beta1) # initializing the parameter values
    
    loglik <- function(z) 
    {
      alpha <- z[1]
      phi = z[2:(p1+1)] 
      prec <- z[p1+2] # precision parameter
      beta <- z[(p1+3):length(z)]
      
      error<-rep(0,n) # E(error)=0 
      eta<-rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) 
        error[i] <- ynew[i]-eta[i] 
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings( log( dkum(y1,mu,prec) ) )
      sum(ll)
    } 
    
    escore <- function(z) 
    {
      
      alpha <- z[1]
      phi <- z[2:(p1+1)]
      prec <- z[p1+2] # precision parameter
      beta <- z[(p1+3):length(z)]
      
      error<- rep(0,n) # E(error)=0 
      eta<- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) ))
        error[i] <- ynew[i]-eta[i] 
      }
      
      ###########
      
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
      for(i in 1:(n-m))
      {
        P[i,] <- ynew[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
      }
      
      k1<- length(beta)
      M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
      for(i in 1:(n-m))
      {
        for(j in 1:k1)
          M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])  
      }
      
      v <- matrix(rep(1,(n-m)),ncol=1)
      rP <- P
      rM <- M
      
      mT <- diag(mu.eta(eta[(m+1):n]))
      
      delta <- log(0.5)/log(1-mu^prec)
      c1 <- mu^(prec-1)/( (1-mu^prec)*log(1-mu^prec) )
      c2 <- delta*log(1-y1^prec)+1
      vc <- prec*c1*c2
      
      a1 <- 1/prec
      a2 <- log(y1)
      a3 <- vc*mu*log(mu)/prec
      a4 <- delta-1
      a5 <- ((y1^prec)*log(y1))/(1-y1^prec)
      a <- as.vector(a1+a2+a3-(a4*a5)) # derivada em relação a phi
      
      Ualpha <- t(v) %*% mT %*% vc
      Uphi <-   t(rP) %*% mT %*% vc
      Uprec <-   sum(a)      
      Ubeta <-  t(rM) %*% mT %*% vc
      
      rval <- c(Ualpha,Uphi,Uprec,Ubeta)
    }
    names_par <- c("alpha",names_phi,"precision",names_beta)
    
    opt <- optim(reg, loglik, escore, 
                 method = "BFGS", 
                 control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      opt <- optim(reg, loglik, #escore, 
                   method = "BFGS", 
                   control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
    
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+2+ncol(X))]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    prec <- coef[p1+2] # precision parameter
    beta <- coef[(p1+3):length(coef)]
    
    z$alpha <- alpha
    z$phi <- phi
    z$prec <- prec
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) 
      errorhat[i] <- ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * prec) - digamma((1 - muhat) * prec)
    
    
    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
    for(i in 1:(n-m))
    {
      P[i,] <- ynew[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
    }
    
    k1<- length(beta)
    M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
    for(i in 1:(n-m))
    {
      for(j in 1:k1)
        M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])  
    }
    
    v <- matrix(rep(1,(n-m)),ncol=1)
    rP <- P
    rM <- M
    
    mT <- diag(mu.eta(etahat[(m+1):n]))
    c <- as.vector((muhat^(prec-1))/((1-muhat^prec)*log(1-muhat^prec)) * ((log(0.5)*log(1-y1^prec) )/log(1-muhat^prec) +1))
    vc <- prec*c
    lambda2 <- (muhat^(2*prec-2))/((1-muhat^prec)^2 * (log(1-muhat^prec))^2)
    W <- diag(as.vector(- prec^2 * lambda2))
    vI <- matrix(rep(1,(n-m)),ncol=1)
    lambda1 <- (muhat^(prec-2))/((1-muhat^prec) * (log(1-muhat^prec)))
    euler <- 0.5772156649
    delta <- log(0.5)/(log(1-muhat^prec))
    d <- -prec*muhat*log(muhat)*lambda2-prec*muhat*delta*lambda1*((1-digamma(delta+1)-euler)/((delta-1)*prec))
    D <- diag(as.vector(d))
    L1 <- 1/(prec^2) 
    L2<- delta*(muhat^2)*lambda2*(log(muhat)^2)*log(1-y1^prec)
    L3<- 2*delta*(muhat^2)*log(muhat)*lambda1*(( (y1^prec) * log(y1))/(1-y1^prec))
    L4<- (delta-1)*((y1^prec)*log(y1)^2)/((1-y1^prec)^2)   
    L5<- c*muhat*(log(muhat)^2)*(lambda1*(muhat^2)+1/(1-muhat^prec))
    
    L<- -L1+L2-L3-L4+L5
    
    Kaa <- t(v) %*% W %*% mT^2 %*% v
    Kap <- t(v) %*% W %*% mT^2 %*% rP
    Kpa <- t(Kap)
    Kaprec <- t(v) %*% D %*% mT %*% vI
    Kpreca <- t(Kaprec)
    
    Kpp <- t(rP) %*% W %*% mT^2 %*% rP
    Kpprec <- t(rP) %*% D %*% mT %*% vI 
    Kprecp <- t(Kpprec)
    
    Kprecprec <- sum(L)
    
    Kab <- t(v) %*% W %*% mT^2 %*% rM
    Kba <- t(Kab)
    Kbb <- t(rM) %*% W %*% mT^2 %*% rM
    Kpb <- t(rP) %*% W %*% mT^2 %*% rM
    Kbp <- t(Kpb)
    Kbprec <- t(rM) %*% D %*% mT %*% vI
    Kprecb <- t(Kbprec)
    
    K <- -rbind(
      cbind(Kaa,Kap,Kaprec,Kab),
      cbind(Kpa,Kpp,Kpprec,Kpb),
      cbind(Kpreca,Kprecp,Kprecprec,Kprecb),
      cbind(Kba,Kbp,Kbprec,Kbb)
    )
    
    z$K <- K
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    X_prev<- rbind(X,X_hat)
    
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta) ))
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 
    }
  }
  
  
  #   ######### KMAX model
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==F))
  {
    p1<-0
    #print("KMAX model",quote=F)
    beta1<- mqo[(2):length(mqo)]
    reg <- c(mqo[1], rep(0,q1), prec, beta1) # initializing the parameter values
    
    loglik <- function(z)
    {
      alpha <- z[1]
      theta = z[(2):(q1+1)]
      prec <- z[q1+2] # precision parameter
      beta <- z[(q1+3):length(z)]
      
      error<-rep(0,n) # E(error)=0
      eta<-rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i]
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings( log( dkum(y1,mu,prec) ) )
      sum(ll)
    }
    
    escore <- function(z)
    {
      alpha <- z[1]
      theta <- z[(2):(q1+1)]
      prec <- z[q1+2] # precision parameter
      beta <- z[(q1+3):length(z)]
      
      error<- rep(0,n) # E(error)=0
      eta<- rep(NA,n)
      
      for(i in (m+1):n)
      {
        eta[i]<- alpha + X[i,]%*%as.matrix(beta) + (theta%*%error[i-ma])
        error[i] <- ynew[i]-eta[i]
      }
      
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m))
      {
        R[i,] <- error[i+m-ma]
      }
      k1<- length(beta)
      M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
      for(i in 1:(n-m))
      {
        for(j in 1:k1)
          M[i,j] <- X[i+m,j] 
      }
      
      ###FB  recorrences
      deta.dalpha <- rep(0,n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      deta.dbeta<- matrix(0, ncol=k1,nrow=n)
      
      for(i in (m+1):n)
      {
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
        deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
      }
      
      v <- deta.dalpha[(m+1):n]
      rR <- deta.dtheta[(m+1):n,]
      rM <- deta.dbeta[(m+1):n,]
      
      mT <- diag(mu.eta(eta[(m+1):n]))
      
      delta <- log(0.5)/log(1-mu^prec)
      c1 <- mu^(prec-1)/( (1-mu^prec)*log(1-mu^prec) )
      c2 <- delta*log(1-y1^prec)+1
      vc <- prec*c1*c2
      
      a1 <- 1/prec
      a2 <- log(y1)
      a3 <- vc*mu*log(mu)/prec
      a4 <- delta-1
      a5 <- ((y1^prec)*log(y1))/(1-y1^prec)
      a <- as.vector(a1+a2+a3-(a4*a5)) # derivada em relação a phi
      
      Ualpha <- t(v) %*% mT %*% vc
      Utheta <- t(rR) %*% mT %*% vc
      Uprec <-   sum(a)      
      Ubeta <-  t(rM) %*% mT %*% vc
      
      rval <- c(Ualpha,Utheta,Uprec,Ubeta)
    }
    names_par <- c("alpha",names_theta,"precision",names_beta)
    
    opt <- optim(reg, loglik, escore, 
                 method = "BFGS", 
                 control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
      opt <- optim(reg, loglik, #escore, 
                   method = "BFGS", 
                   control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))
      if (opt$conv != 0)
      {
        warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
      }else{
        warning("IT WORKS WITH NUMERICAL GRADIENT!")
      }
    }
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(q1+2+ncol(X) )]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    theta <- coef[(2):(q1+1)]
    prec <- coef[q1+2] # precision parameter
    beta <- coef[(q1+3):length(coef)]
    
    z$alpha <- alpha
    z$theta <- theta
    z$prec <- prec
    
    errorhat<-rep(0,n) # E(error)=0
    etahat<-rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + X[i,]%*%as.matrix(beta) + (theta%*%errorhat[i-ma])
      errorhat[i] <- ynew[i]-etahat[i] # predictor scale
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * prec) - digamma((1 - muhat) * prec)
    
    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      R[i,] <- errorhat[i+m-ma]
    }
    k1<- length(beta)
    M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
    for(i in 1:(n-m))
    {
      for(j in 1:k1)
        M[i,j] <- X[i+m,j] 
    }
    
    
    ###FB  recorrences
    deta.dalpha <- rep(0,n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    deta.dbeta<- matrix(0, ncol=k1,nrow=n)
    
    for(i in (m+1):n)
    {
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
    }
    
    v <- deta.dalpha[(m+1):n]
    rR <- deta.dtheta[(m+1):n,]
    rM <- deta.dbeta[(m+1):n,]
    
    
    mT <- diag(mu.eta(etahat[(m+1):n]))
    c <- as.vector((muhat^(prec-1))/((1-muhat^prec)*log(1-muhat^prec)) * ((log(0.5)*log(1-y1^prec) )/log(1-muhat^prec) +1))
    vc <- prec*c
    lambda2 <- (muhat^(2*prec-2))/((1-muhat^prec)^2 * (log(1-muhat^prec))^2)
    W <- diag(as.vector(- prec^2 * lambda2))
    vI <- matrix(rep(1,(n-m)),ncol=1)
    lambda1 <- (muhat^(prec-2))/((1-muhat^prec) * (log(1-muhat^prec)))
    euler <- 0.5772156649
    delta <- log(0.5)/(log(1-muhat^prec))
    d <- -prec*muhat*log(muhat)*lambda2-prec*muhat*delta*lambda1*((1-digamma(delta+1)-euler)/((delta-1)*prec))
    D <- diag(as.vector(d))
    L1 <- 1/(prec^2) 
    L2<- delta*(muhat^2)*lambda2*(log(muhat)^2)*log(1-y1^prec)
    L3<- 2*delta*(muhat^2)*log(muhat)*lambda1*(( (y1^prec) * log(y1))/(1-y1^prec))
    L4<- (delta-1)*((y1^prec)*log(y1)^2)/((1-y1^prec)^2)   
    L5<- c*muhat*(log(muhat)^2)*(lambda1*(muhat^2)+1/(1-muhat^prec))
    
    L<- -L1+L2-L3-L4+L5
    
    Kaa <- t(v) %*% W %*% mT^2 %*% v
    Kat <- t(v) %*% W %*% mT^2 %*% rR
    Kta <- t(Kat)
    Kaprec <- t(v) %*% D %*% mT %*% vI
    Kpreca <- t(Kaprec)
    
    Ktt <- t(rR) %*% W %*% mT^2 %*% rR
    Ktprec <- t(rR) %*% D %*% mT %*% vI
    Kprect <- t(Ktprec)
    
    Kprecprec <- sum(L)
    
    Kab <- t(v) %*% W %*% mT^2 %*% rM
    Kba <- t(Kab)
    Kbb <- t(rM) %*% W %*% mT^2 %*% rM
    Ktb <- t(rR) %*% W %*% mT^2 %*% rM
    Kbt <- t(Ktb)
    Kbprec <- t(rM) %*% D %*% mT %*% vI
    Kprecb <- t(Kbprec)
    
    
    K <- -rbind(
      cbind(Kaa,Kat,Kaprec,Kab),
      cbind(Kta,Ktt,Ktprec,Ktb),
      cbind(Kpreca,Kprect,Kprecprec,Kprecb),
      cbind(Kba,Kbt,Kbprec,Kbb)
    )
    
    z$K <- K
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    X_prev<- rbind(X,X_hat)#X[i - ar, ]
    
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 
    }
  } 
  
  z$serie <- y
  z$barma <- names_par
  z$forecast <- y_prev[(n+1):(n+h1)]
  
  ##############################################
  
  # #library("rootSolve")
  # print("AQUI1")
  # escores<-rbind(escore(z$coef),gradient(loglik, z$coef))
  # print(escores)
  # print(gradient(escore, z$coef))
  # print(-z$K)
  
  # residuals
  res1 <- y-z$fitted
  vary <- (b-a)^2 * ( (log(0.5)/log(1-z$fitted^prec)) * beta(1+2/prec, log(0.5)/log(1-z$fitted^prec))
                      - (log(0.5)/log(1-z$fitted^prec) * beta(1+1/prec, log(0.5)/log(1-z$fitted^prec)))^2 )  
  
  z$resid1 <- (res1/sqrt(vary))[(m+1):n]
  
  l_tilde <- log(dkum(y,y,z$prec,a=a,b=b))
  l_hat <- log(dkum(y,z$fitted,z$prec,a=a,b=b))
  
  dt <- (l_tilde-l_hat)[(m+1):n]
  dt[which(dt<0)]<-0
  
  z$l_hat <- l_hat
  
  z$resid2 <- sign(y[(m+1):n]-z$fitted[(m+1):n])*sqrt(2*(dt))
  
  z$resid3 <- as.vector(qnorm(pkum(y[(m+1):n],z$fitted[(m+1):n],z$prec,a=a,b=b)))
  
  if(resid==1) residc <- z$resid1
  if(resid==2) residc <- z$resid2
  if(resid==3) residc <- z$resid3
  
  ############################################
  
  #vcov <- try(solve(K))
  Kchol<- tryCatch(chol(z$K), error = function(e) return("error"), warning = function(o) return("error"))
  
  if(Kchol[1] == "error")
  { 
    z$vcov <- try(solve(K))
    warning("We have problems with information matrix inversion!")

  }else{
    vcov <- try(chol2inv(Kchol))
    z$vcov <- vcov
  }
  
  #z$vcov <- vcov
  
  stderror <- sqrt(diag(z$vcov))
  z$stderror <- stderror
  
  z$zstat <- abs(z$coef/stderror)
  z$pvalues <- 2*(1 - pnorm(z$zstat) )
  
  z$loglik <- opt$value
  #z$loglik <- opt$value*(n/(n-m))
  z$counts <- as.numeric(opt$counts[1])
  
  if(any(is.na(X)==F))
  {
    z$k<- (p1+q1+2+length(beta))
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }else{
    z$k<- (p1+q1+2)
    z$aic <- -2*(z$loglik*(n/(n-m)))+2*(z$k)
    z$bic <- -2*(z$loglik*(n/(n-m)))+log(n)*(z$k)
    z$hq <- -2*(z$loglik*(n/(n-m)))+log(log(n))*(z$k)
  }
  
  model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
  
  z$model <- model_presentation
  z$link <- link
  
  
  ###################################################
  ######### GRAPHICS ################################
  
  if(diag>0)
  {
    
    print(model_presentation)
    print(" ",quote=F)
    print(c("Log-likelihood:",round(z$loglik,4)),quote=F)
    print(c("Number of iterations in BFGS optim:",z$counts),quote=F)
    print(c("AIC:",round(z$aic,4)," SIC:",round(z$bic,4)," HQ:",round(z$hq,4)),quote=F)
    
    print("Residuals:",quote=F)
    print(summary(residc))
    
    t<-seq(-5,n+6,by=1)
    
    par(mfrow=c(1,1))
    par(mar=c(2.8, 2.7, 1.2, 1)) 
    par(mgp=c(1.7, 0.45, 0))
    plot(residc,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
    lines(t,rep(-3,n+12),lty=2,col=1)
    lines(t,rep(3,n+12),lty=2,col=1)
    lines(t,rep(-2,n+12),lty=3,col=1)
    lines(t,rep(2,n+12),lty=3,col=1)
    
    
    max_y<- max(c(z$fitted,y),na.rm=T)
    min_y<- min(c(z$fitted,y),na.rm=T)
    plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
         xlab="Fitted values",ylab="Observed data",
         xlim=c(0.95*min_y,max_y*1.05),
         ylim=c(0.95*min_y,max_y*1.05))
    lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
    
    plot(as.vector(z$fitted[(m+1):n]),as.vector(residc), main=" ", pch = "+",
         xlab="Fitted values",ylab="Residuals")
    #lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
    
    
    densidade<-density(residc)
    plot(densidade,ylab="density",main=" ")
    lines(densidade$x,dnorm(densidade$x),lty=2)
    legend("topleft",c("Exact distribution of residuals","Normal approximation"),#pch=vpch,
           pt.bg="white", lty=c(1,2), bty="n")
    
    acf(residc,ylab="ACF",xlab="Lag") 
    
    pacf(residc,ylab="PACF",xlab="Lag") 
    
    max_r<- max(residc,na.rm=T)
    min_r<- min(residc,na.rm=T)
    qqnorm(residc, pch = "+",
           xlim=c(0.95*min_r,max_r*1.05),
           ylim=c(0.95*min_r,max_r*1.05),
           main="",xlab="Normal quantiles",ylab="Empirical quantiles")
    lines(c(-10,10),c(-10,10),lty=2)
    
    par(mfrow=c(1,1))
    plot(y,type="l",ylab="serie",xlab="tempo")
    lines(z$fitted,col="red")
    
    fim<-end(y)[1]+end(y)[2]/12
    
    y_prev <- ts(y_prev, start=start(y), frequency=frequency(y))
    par(mfrow=c(1,1))
    plot(y_prev,type="l",col="red", ylim=c(min(y),max(y)),ylab="Serie",xlab="Time")
    abline(v=fim,lty=2)
    lines(y)
    
    w1<-5
    h1<-4
    
    if(diag>1)
    {
      postscript(file = "resid_v_ind.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1))
        par(mgp=c(1.7, 0.45, 0))
        plot(residc,main=" ",xlab="Index",ylab="Residuals", pch = "+",ylim=c(-4,4))
        lines(t,rep(-3,n+12),lty=2,col=1)
        lines(t,rep(3,n+12),lty=2,col=1)
        lines(t,rep(-2,n+12),lty=3,col=1)
        lines(t,rep(2,n+12),lty=3,col=1)
      }
      dev.off()
      
      postscript(file = "resid_v_fitted.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        plot(as.vector(z$fitted[(m+1):n]),as.vector(residc), main=" ", pch = "+",
             xlab="Fitted values",ylab="Residuals",ylim=c(-4,4))
        lines(t,rep(-3,n+12),lty=2,col=1)
        lines(t,rep(3,n+12),lty=2,col=1)
        lines(t,rep(-2,n+12),lty=3,col=1)
        lines(t,rep(2,n+12),lty=3,col=1)
      }
      dev.off()
      
      postscript(file = "obs_v_fit.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
             xlab="Fitted values",ylab="Observed data",
             xlim=c(0.95*min_y,max_y*1.05),
             ylim=c(0.95*min_y,max_y*1.05))
        lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
      }
      dev.off()
      
      postscript(file = "resid_density.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(1.5, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        
        plot(densidade,ylab="Density",main=" ",xlab=" ",ylim=c(0,1.15*max(densidade$y)))
        lines(densidade$x,dnorm(densidade$x),lty=2)
        legend("topleft",c("Exact distribution of residuals","Normal approximation"),#pch=vpch,
               pt.bg="white", lty=c(1,2), bty="n")
      }
      dev.off()
      
      postscript(file = "resid_FAC.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        acf(residc,ylab="ACF",xlab="Lag") 
      }
      dev.off()
      
      postscript(file = "resid_FACP.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        pacf(residc,ylab="PACF",xlab="Lag")
      }
      dev.off()
      
      postscript(file = "qq_plot.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {  
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) 
        par(mgp=c(1.7, 0.45, 0))
        qqnorm(residc, pch = "+",
               xlim=c(0.95*min_r,max_r*1.05),
               ylim=c(0.95*min_r,max_r*1.05),
               main="",xlab="Normal quantiles",ylab="Empirical quantiles")
        lines(c(-10,10),c(-10,10),lty=2)
      }
      dev.off()
      
      postscript(file = "adjusted.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        plot(y,type="l",ylab="Serie",xlab="Time")
        lines(z$fitted,col="red")
      }
      dev.off()
      

      postscript(file = "forecast.eps",horizontal=F,paper="special",width = 6, height = 4.7,family = "Times")
      {
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        plot(y_prev,type="l",lty=2,col="red", ylim=c(min(y),max(y)),ylab="RH",xlab="Times")
        abline(v=fim,lty=2)
        lines(y)
        legend("bottomleft",c("Observed data","Fitted and forecast values"),#pch=vpch,
               pt.bg="white", lty=c(1,2), bty="n",col=c(1,"red"))
      }
      dev.off()
      
    }    
  }  

  return(z)
}