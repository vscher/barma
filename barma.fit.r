# Implementado por Fabio M Bayer (bayer@ufsm.br) em 15/10/2015

barma.fit<- function (y, ar, ma, link, names_phi,names_theta,names_beta,diag,h1,X,X_hat,resid)
{
  maxit1<-1000
  
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
  
  # inicializacao dos parametros alpha e phi (beta)
  if(any(is.na(ar)==F)) # se nao tem componente ar
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
  
  if(any(is.na(X)==T)) # nao tem regressores
  {
    x <- as.matrix(Z)
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
    prec = 1/n1 * sum(mean * (1 - mean)/sigma2 - 1)
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
    prec = 1/n1 * sum(mean * (1 - mean)/sigma2 - 1)
    
  }
  
  
  
  
  ############
  
  ######### ARMA model
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==T))
  { 
    #print("BARMA model",quote=F)
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
        #error[i] <- y[i]-linkinv(eta[i])
        #mui<-linkinv(eta[i])
        #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE))
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
        #error[i] <- y[i]-linkinv(eta[i])
        #mui<-linkinv(eta[i])
        #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
      }
      
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      ystar <- log(y1/(1-y1))
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
      
      mT <- diag(mu.eta(mu))
      
      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m))
      {
        R[i,] <- error[i+m-ma]
      }
      
      ###FB  recorrencias
      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)

      # print(dim(P))
      # print(dim(deta.dphi))
      
      for(i in (m+1):n)
      {
        #print(P[(i-m),(i-m)])
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
        
      }
      
      a <- deta.dalpha[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rR <- deta.dtheta[(m+1):n,]

      Ualpha <- prec * a %*% mT %*% (ystar-mustar)
      Uphi <-   prec * t(rP) %*% mT %*% (ystar-mustar)
      Utheta <- prec * t(rR) %*% mT %*% (ystar-mustar)
      Uprec <-  sum(mu * (ystar-mustar) + log(1-y1) 
                    - digamma((1 - mu) * prec) + digamma(prec) )
      
      rval <- c(Ualpha,Uphi,Utheta,Uprec)
    }
    names_par <- c("alpha",names_phi,names_theta,"precision")
    
    opt <- optim(reg, loglik, escore, method = "BFGS", 
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
    
    z <- c()
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
      #errorhat[i] <- y[i]-linkinv(etahat[i]) # original scale
      #mui<-linkinv(etahat[i])
      #errorhat[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
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
    
    vI <- as.vector(rep(1,n-m)) 
    
    
    ###FB  recorrencias
    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    
    # print(dim(P))
    # print(dim(deta.dphi))
    
    for(i in (m+1):n)
    {
      #print(P[(i-m),(i-m)])
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      
    }
    
    a <- deta.dalpha[(m+1):n]
    rP <- deta.dphi[(m+1):n,]
    rR <- deta.dtheta[(m+1):n,]
    
    psi1 = trigamma(muhat * prec)
    psi2 = trigamma((1 - muhat) * prec)
    mT <- diag(mu.eta(muhat))
    W = diag(c(prec * (psi1 + psi2))) %*% mT^2
    vc = prec * (psi1 * muhat - psi2 * (1 - muhat))
    D = diag(c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(prec)))
    
  
    Kaa <- prec* t(a) %*% W %*% a # r ok
    Kpa <- prec* t(rP) %*% W %*% a # r ok
    Kap <- t(Kpa) # ok
    Kta <- prec * t(rR) %*% W %*% a # ok
    Kat <- t(Kta) # ok 
    Kaprec <- t(a) %*% mT %*% vc # ok
    Kpreca <- Kaprec #ok 
    Kpp <- prec * t(rP) %*% W %*% rP
    Kpt <- prec * t(rP) %*% W %*% rR
    Ktp <- t(Kpt)
    Kpprec <- t(rP) %*% mT %*% vc
    Kprecp <- t(Kpprec)
    Ktt <- prec * t(rR) %*% W %*% rR
    Ktprec <- t(rR) %*% mT %*% vc
    Kprect <- t(Ktprec)
    Kprecprec <- sum(diag(D))
    
    K <- rbind(
      cbind(Kaa,Kap,Kat,Kaprec),
      cbind(Kpa,Kpp,Kpt,Kpprec),
      cbind(Kta,Ktp,Ktt,Ktprec),
      cbind(Kpreca,Kprecp,Kprect,Kprecprec)
    )
    
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar]) + (theta%*%errorhat[n+i-ma])
      #errorhat[i]<-ynew[i]-etahat[i] # predictor scale
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 # original scale
    }
  }  
  
  
  ##### only AR model
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==T))
  {
    #print("BAR model",quote=F)
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
      
      ll <- suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE))
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
      ystar <- log(y1/(1-y1))
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
      
      mT <- diag(mu.eta(mu))
      
      Ualpha <- prec * sum((mu.eta(mu)) * (ystar-mustar))
      Uphi <-   prec * t(P) %*% mT %*% (ystar-mustar)
      Uprec <-  sum(mu * (ystar-mustar) + log(1-y1) 
                    - digamma((1 - mu) * prec) + digamma(prec) )
      
      rval <- c(Ualpha,Uphi,Uprec)
    }
    names_par <- c("alpha",names_phi,"precision")
    
    opt <- optim(reg, loglik, escore, method = "BFGS", 
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
    
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
      #errorhat[i] <- y[i]-linkinv(etahat[i]) # original scale
      #mui<-linkinv(etahat[i])
      #errorhat[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * prec) - digamma((1 - muhat) * prec)
    
    #z$resid1 <- errorhat
    
    vI <- as.vector(rep(1,n-m)) 
    
    psi1 = trigamma(muhat * prec)
    psi2 = trigamma((1 - muhat) * prec)
    mT <- diag(mu.eta(muhat))
    W = diag(c(prec * (psi1 + psi2))) %*% mT^2
    vc = prec * (psi1 * muhat - psi2 * (1 - muhat))
    D = diag(c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(prec)))
    
    Kaa <- as.matrix(prec * sum(diag(W)))
    Kpa <- prec* t(P) %*% W %*% vI
    Kap <- t(Kpa)
    Kaprec <- vI %*% mT %*% vc
    Kpreca <- Kaprec
    Kpp <- prec * t(P) %*% W %*% P
    Kpprec <- t(P) %*% mT %*% vc
    Kprecp <- t(Kpprec)
    Kprecprec <- sum(diag(D))
    
    K <- rbind(
      cbind(Kaa,Kap,Kaprec),
      cbind(Kpa,Kpp,Kpprec),
      cbind(Kpreca,Kprecp,Kprecprec)
    )
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar]) 
      #errorhat[i]<-ynew[i]-etahat[i] # predictor scale
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
    }
    
  }
  
  ######### MA model
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==T))
  { 
    p1<-0
    #print("BMA model",quote=F)
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
        #error[i] <- y[i]-linkinv(eta[i])
        #mui<-linkinv(eta[i])
        #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE))
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
        #error[i] <- y[i]-linkinv(eta[i])
        #mui<-linkinv(eta[i])
        #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
      }
      
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      ystar <- log(y1/(1-y1))
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
      
      mT <- diag(mu.eta(mu))
      
      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      
      for(i in 1:(n-m))
      {
        R[i,] <- error[i+m-ma]
      }
      
      ###FB  recorrencias
      deta.dalpha <- rep(0,n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      
      # print(dim(P))
      # print(dim(deta.dphi))
      
      for(i in (m+1):n)
      {
        #print(P[(i-m),(i-m)])
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
        
      }
      
      a <- deta.dalpha[(m+1):n]
      rR <- deta.dtheta[(m+1):n,]
      
      Ualpha <- prec * a %*% mT %*% (ystar-mustar)
      Utheta <- prec * t(rR) %*% mT %*% (ystar-mustar)
      Uprec <-  sum(mu * (ystar-mustar) + log(1-y1) 
                    - digamma((1 - mu) * prec) + digamma(prec) )
      
      
      rval <- c(Ualpha,Utheta,Uprec)
    }
    names_par <- c("alpha",names_theta,"precision")
    
    opt <- optim(reg, loglik, escore, method = "BFGS", 
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
    
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
      #errorhat[i] <- y[i]-linkinv(etahat[i]) # original scale
      #mui<-linkinv(etahat[i])
      #errorhat[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
    }
    muhat <- linkinv(etahat[(m+1):n])
    y1<-y[(m+1):n]
    
    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    z$mustarhat <- digamma(muhat * prec) - digamma((1 - muhat) * prec)
    #z$resid1 <- errorhat
    
    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m))
    {
      R[i,] <- errorhat[i+m-ma]
    }
    
    vI <- as.vector(rep(1,n-m)) 
    
    ###FB  recorrencias
    deta.dalpha <- rep(0,n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    
    # print(dim(P))
    # print(dim(deta.dphi))
    
    for(i in (m+1):n)
    {
      #print(P[(i-m),(i-m)])
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      
    }
    
    a <- deta.dalpha[(m+1):n]
    rR <- deta.dtheta[(m+1):n,]
    
    psi1 = trigamma(muhat * prec)
    psi2 = trigamma((1 - muhat) * prec)
    mT <- diag(mu.eta(muhat))
    W = diag(c(prec * (psi1 + psi2))) %*% mT^2
    vc = prec * (psi1 * muhat - psi2 * (1 - muhat))
    D = diag(c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(prec)))
    
    Kaa <- prec* t(a) %*% W %*% a # r ok
    Kta <- prec * t(rR) %*% W %*% a # ok
    Kat <- t(Kta) # ok 
    Kaprec <- t(a) %*% mT %*% vc # ok
    Kpreca <- Kaprec #ok 
    Ktt <- prec * t(rR) %*% W %*% rR
    Ktprec <- t(rR) %*% mT %*% vc
    Kprect <- t(Ktprec)
    Kprecprec <- sum(diag(D))
    

    K <- rbind(
      cbind(Kaa,Kat,Kaprec),
      cbind(Kta,Ktt,Ktprec),
      cbind(Kpreca,Kprect,Kprecprec)
    )
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + (theta%*%errorhat[n+i-ma])
      #errorhat[i]<-ynew[i]-etahat[i] # predictor scale
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 # original scale
    }
  }
  
  

  
  
  ######### BARMAX model
  if(any(is.na(ar)==F) && any(is.na(ma)==F) && any(is.na(X)==F))
  { 
    print("BARMAX model",quote=F)
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
        #error[i] <- y[i]-linkinv(eta[i])
        #mui<-linkinv(eta[i])
        #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE))
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
        #error[i] <- y[i]-linkinv(eta[i])
        #mui<-linkinv(eta[i])
        #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
      }
      
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      ystar <- log(y1/(1-y1))
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
      
      mT <- diag(mu.eta(mu))
      
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
        for(j in 1:k1)
          M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])  
      }
      
      
      ###FB  recorrencias
      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      deta.dbeta<- matrix(0, ncol=k1,nrow=n)
      
      # print(dim(P))
      # print(dim(deta.dphi))
      
      for(i in (m+1):n)
      {
        #print(P[(i-m),(i-m)])
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
        deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
        
      }
      
      a <- deta.dalpha[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rR <- deta.dtheta[(m+1):n,]
      rM <- deta.dbeta[(m+1):n,]
      
      Ualpha <- prec * a %*% mT %*% (ystar-mustar)
      Uphi <-   prec * t(rP) %*% mT %*% (ystar-mustar)
      Utheta <- prec * t(rR) %*% mT %*% (ystar-mustar)
      Uprec <-  sum(mu * (ystar-mustar) + log(1-y1) 
                    - digamma((1 - mu) * prec) + digamma(prec) )
      Ubeta <-   prec * t(rM) %*% mT %*% (ystar-mustar)
      
      rval <- c(Ualpha,Uphi,Utheta,Uprec,Ubeta)
    }
    names_par <- c("alpha",names_phi,names_theta,"precision",names_beta)
    
    opt <- optim(reg, loglik, escore, method = "BFGS", 
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
    
    z <- c()
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
      #errorhat[i] <- y[i]-linkinv(etahat[i]) # original scale
      #mui<-linkinv(etahat[i])
      #errorhat[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
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
    
    k1<-length(beta)
    M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
    for(i in 1:(n-m))
    {
      for(j in 1:k1)
        M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])  
    }
    
    
    
    ###FB  recorrencias
    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    deta.dbeta<- matrix(0, ncol=k1,nrow=n)
    
    # print(dim(P))
    # print(dim(deta.dphi))
    
    for(i in (m+1):n)
    {
      #print(P[(i-m),(i-m)])
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
      
    }
    
    a <- deta.dalpha[(m+1):n]
    rP <- deta.dphi[(m+1):n,]
    rR <- deta.dtheta[(m+1):n,]
    rM <- deta.dbeta[(m+1):n,]
    
#    vI <- as.vector(rep(1,n-m)) 
    
    psi1 = trigamma(muhat * prec)
    psi2 = trigamma((1 - muhat) * prec)
    mT <- diag(mu.eta(muhat))
    W = diag(c(prec * (psi1 + psi2))) %*% mT^2
    vc = prec * (psi1 * muhat - psi2 * (1 - muhat))
    D = diag(c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(prec)))
    
    
    Kaa <- prec* t(a) %*% W %*% a # r ok
    Kpa <- prec* t(rP) %*% W %*% a # r ok
    Kap <- t(Kpa) # ok
    Kta <- prec * t(rR) %*% W %*% a # ok
    Kat <- t(Kta) # ok 
    Kaprec <- t(a) %*% mT %*% vc # ok
    Kpreca <- Kaprec #ok 
    Kpp <- prec * t(rP) %*% W %*% rP
    Kpt <- prec * t(rP) %*% W %*% rR
    Ktp <- t(Kpt)
    Kpprec <- t(rP) %*% mT %*% vc
    Kprecp <- t(Kpprec)
    Ktt <- prec * t(rR) %*% W %*% rR
    Ktprec <- t(rR) %*% mT %*% vc
    Kprect <- t(Ktprec)
    Kprecprec <- sum(diag(D))
    
    Kba <- prec * t(rM)%*%W %*% a
    Kbb <- prec * t(rM) %*% W %*% rM
    Kbprec <- t(rM) %*% mT %*% vc 
    Kbp <- prec * t(rM) %*% W %*% rP
    Kbt <- prec * t(rM) %*% W %*% rR
    
    Kab <- t(Kba)
    Kprecb <- t(Kbprec)
    Kpb <- t(Kbp)
    Ktb <- t(Kbt)
    
    K <- rbind(
      cbind(Kaa,Kap,Kat,Kaprec,Kab),
      cbind(Kpa,Kpp,Kpt,Kpprec,Kpb),
      cbind(Kta,Ktp,Ktt,Ktprec,Ktb),
      cbind(Kpreca,Kprecp,Kprect,Kprecprec,Kprecb),
      cbind(Kba,Kbp,Kbt,Kbprec,Kbb)
    )
    
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    X_prev<- rbind(X,X_hat)
    
    #     (phi%*%(ynew[i-ar]-X[i-ar,]%*%beta ))
    #     (phi%*%(ynew_prev[n+i-ar]-X_prev[i-ar,]%*%beta ))
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta) )) + (theta%*%errorhat[n+i-ma])
      #errorhat[i]<-ynew[i]-etahat[i] # predictor scale
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 # original scale
    }
  }  # fim BARMAX
  
  
  
  
  
  
  
  ######### BARX model
  if(any(is.na(ar)==F) && any(is.na(ma)==T) && any(is.na(X)==F))
  { 
    q1<-0
    print("BARX model",quote=F)
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
        #error[i] <- y[i]-linkinv(eta[i])
        #mui<-linkinv(eta[i])
        #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE))
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
        #error[i] <- y[i]-linkinv(eta[i])
        #mui<-linkinv(eta[i])
        #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
      }
      
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      ystar <- log(y1/(1-y1))
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
      
      mT <- diag(mu.eta(mu))
      
      P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
      for(i in 1:(n-m))
      {
        P[i,] <- ynew[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
      }
      
      M <- matrix(rep(NA,(n-m)*length(beta)),ncol=length(beta))
      for(i in 1:(n-m))
      {
        for(j in 1:length(beta))
          M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])  
      }
      
      Ualpha <- prec * sum((mu.eta(mu)) * (ystar-mustar))
      Uphi <-   prec * t(P) %*% mT %*% (ystar-mustar)
      Uprec <-  sum(mu * (ystar-mustar) + log(1-y1) 
                    - digamma((1 - mu) * prec) + digamma(prec) )
      Ubeta <-   prec * t(M) %*% mT %*% (ystar-mustar)
      
      rval <- c(Ualpha,Uphi,Uprec,Ubeta)
    }
    names_par <- c("alpha",names_phi,"precision",names_beta)
    
    opt <- optim(reg, loglik, escore, method = "BFGS", 
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    
    print("AQUI")
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
    
    z <- c()
    z$conv <- opt$conv
    coef <- (opt$par)[1:(p1+2+ncol(X))]
    names(coef)<-names_par
    z$coeff <- coef
    
    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    prec <- coef[p1+2] # precision parameter
    beta <- coef[(p1+3):length(coef)]
    
    print(coef)
    print(beta)
    
    z$alpha <- alpha
    z$phi <- phi
    z$prec <- prec
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta) )) 
      errorhat[i] <- ynew[i]-etahat[i] # predictor scale
      #errorhat[i] <- y[i]-linkinv(etahat[i]) # original scale
      #mui<-linkinv(etahat[i])
      #errorhat[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
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
    
    M <- matrix(rep(NA,(n-m)*length(beta)),ncol=length(beta))
    for(i in 1:(n-m))
    {
      for(j in 1:length(beta))
        M[i,j] <- X[i+m,j]-sum(phi*X[i+m-ar,j])  
    }
    
    
    vI <- as.vector(rep(1,n-m)) 
    
    psi1 = trigamma(muhat * prec)
    psi2 = trigamma((1 - muhat) * prec)
    mT <- diag(mu.eta(muhat))
    W = diag(c(prec * (psi1 + psi2))) %*% mT^2
    vc = prec * (psi1 * muhat - psi2 * (1 - muhat))
    D = diag(c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(prec)))
    
    Kaa <- as.matrix(prec * sum(diag(W)))
    Kpa <- prec* t(P) %*% W %*% vI
    Kap <- t(Kpa)
    Kaprec <- vI %*% mT %*% vc
    Kpreca <- Kaprec
    Kpp <- prec * t(P) %*% W %*% P
    Kpprec <- t(P) %*% mT %*% vc
    Kprecp <- t(Kpprec)
    Kprecprec <- sum(diag(D))
    
    Kba <- prec * t(M)%*%W %*% vI
    Kbb <- prec * t(M) %*% W %*% M
    Kbprec <- t(M) %*% mT %*% vc 
    Kbp <- prec * t(M) %*% W %*% P
    
    Kab <- t(Kba)
    Kprecb <- t(Kbprec)
    Kpb <- t(Kbp)
    
    K <- rbind(
      cbind(Kaa,Kap,Kaprec,Kab),
      cbind(Kpa,Kpp,Kpprec,Kpb),
      cbind(Kpreca,Kprecp,Kprecprec,Kprecb),
      cbind(Kba,Kbp,Kbprec,Kbb)
    )
    
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    X_prev<- rbind(X,X_hat)
    
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta) ))
      #errorhat[i]<-ynew[i]-etahat[i] # predictor scale
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 # original scale
    }
  }  # fim BARX
  
  
  
  
  
  
  
  
  ######### BMAX model
  if(any(is.na(ar)==T) && any(is.na(ma)==F) && any(is.na(X)==F))
  { 
    p1<-0
    print("BMAX model",quote=F)
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
        #error[i] <- y[i]-linkinv(eta[i])
        #mui<-linkinv(eta[i])
        #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
      }
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      
      ll <- suppressWarnings(dbeta(y1, mu * prec, (1 - mu) * prec, log = TRUE))
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
        #error[i] <- y[i]-linkinv(eta[i])
        #mui<-linkinv(eta[i])
        #error[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
      }
      
      mu <- linkinv(eta[(m+1):n])
      y1<-y[(m+1):n]
      ystar <- log(y1/(1-y1))
      mustar <- digamma(mu * prec) - digamma((1 - mu) * prec)
      
      mT <- diag(mu.eta(mu))
      
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
      
      
      ###FB  recorrencias
      deta.dalpha <- rep(0,n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      deta.dbeta<- matrix(0, ncol=k1,nrow=n)
      
      # print(dim(P))
      # print(dim(deta.dphi))
      
      for(i in (m+1):n)
      {
        #print(P[(i-m),(i-m)])
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
        deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
        
      }
      
      a <- deta.dalpha[(m+1):n]
      rR <- deta.dtheta[(m+1):n,]
      rM <- deta.dbeta[(m+1):n,]
      
      Ualpha <- prec * a %*% mT %*% (ystar-mustar)
      Utheta <- prec * t(rR) %*% mT %*% (ystar-mustar)
      Uprec <-  sum(mu * (ystar-mustar) + log(1-y1) 
                    - digamma((1 - mu) * prec) + digamma(prec) )
      Ubeta <-   prec * t(rM) %*% mT %*% (ystar-mustar)
      
      rval <- c(Ualpha,Utheta,Uprec,Ubeta)
    }
    names_par <- c("alpha",names_theta,"precision",names_beta)
    
    opt <- optim(reg, loglik, escore, method = "BFGS", 
                 control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    
    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")
    
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
      #errorhat[i] <- y[i]-linkinv(etahat[i]) # original scale
      #mui<-linkinv(etahat[i])
      #errorhat[i] <- (y[i]-mui)/sqrt((mui*(1-mui))/(1+prec)) # standardized
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
    
    k1<-length(beta)
    M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
    for(i in 1:(n-m))
    {
      for(j in 1:k1)
        M[i,j] <- X[i+m,j]
    }
    
    
    
    ###FB  recorrencias
    deta.dalpha <- rep(0,n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    deta.dbeta<- matrix(0, ncol=k1,nrow=n)
    
    # print(dim(P))
    # print(dim(deta.dphi))
    
    for(i in (m+1):n)
    {
      #print(P[(i-m),(i-m)])
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
      
    }
    
    a <- deta.dalpha[(m+1):n]
    rR <- deta.dtheta[(m+1):n,]
    rM <- deta.dbeta[(m+1):n,]
    
    #    vI <- as.vector(rep(1,n-m)) 
    
    psi1 = trigamma(muhat * prec)
    psi2 = trigamma((1 - muhat) * prec)
    mT <- diag(mu.eta(muhat))
    W = diag(c(prec * (psi1 + psi2))) %*% mT^2
    vc = prec * (psi1 * muhat - psi2 * (1 - muhat))
    D = diag(c(psi1 * (muhat^2) + psi2 * (1 - muhat)^2 - trigamma(prec)))
    
    
    Kaa <- prec* t(a) %*% W %*% a # r ok
    Kta <- prec * t(rR) %*% W %*% a # ok
    Kat <- t(Kta) # ok 
    Kaprec <- t(a) %*% mT %*% vc # ok
    Kpreca <- Kaprec #ok 
    Ktt <- prec * t(rR) %*% W %*% rR
    Ktprec <- t(rR) %*% mT %*% vc
    Kprect <- t(Ktprec)
    Kprecprec <- sum(diag(D))
    
    Kba <- prec * t(rM)%*%W %*% a
    Kbb <- prec * t(rM) %*% W %*% rM
    Kbprec <- t(rM) %*% mT %*% vc 
    Kbt <- prec * t(rM) %*% W %*% rR
    
    Kab <- t(Kba)
    Kprecb <- t(Kbprec)
    Ktb <- t(Kbt)
    
    K <- rbind(
      cbind(Kaa,Kat,Kaprec,Kab),
      cbind(Kta,Ktt,Ktprec,Ktb),
      cbind(Kpreca,Kprect,Kprecprec,Kprecb),
      cbind(Kba,Kbt,Kbprec,Kbb)
    )
    
    
    #### Forecasting
    ynew_prev <- c(ynew,rep(NA,h1))
    y_prev[1:n] <- z$fitted
    
    X_prev<- rbind(X,X_hat)#X[i - ar, ]
    
    #     (phi%*%(ynew[i-ar]-X[i-ar,]%*%beta ))
    #     (phi%*%(ynew_prev[n+i-ar]-X_prev[i-ar,]%*%beta ))
    
    for(i in 1:h1)
    {
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (theta%*%errorhat[n+i-ma])
      #errorhat[i]<-ynew[i]-etahat[i] # predictor scale
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 # original scale
    }
  }  # fim BMAX
  
  
  z$serie <- y
  z$barma <- names_par
  z$forecast <- y_prev[(n+1):(n+h1)]
  
  ##############################################
  
  # residuals
  res1 <- y-z$fitted
  res2 <- ynew-z$etahat
  res3 <- as.vector(ystar)[(m+1):n]-z$mustarhat
  z$resid1 <- (res1[(m+1):n])/sqrt(z$fitted[(m+1):n]*(1-z$fitted[(m+1):n])/(1+z$prec))
  z$resid2 <- (res2[(m+1):n])/sqrt((mu.eta(z$fitted[(m+1):n])^{-2})*(z$fitted[(m+1):n]*(1-z$fitted[(m+1):n])/(1+z$prec)) )
  z$resid3 <- (res3/sqrt( trigamma(z$fitted[(m+1):n]*z$prec)+trigamma((1-z$fitted[(m+1):n])*z$prec) ) )
  
  l_tilde <- (dbeta(y[(m+1):n], y[(m+1):n] * z$prec, (1 - y[(m+1):n]) * z$prec, log = TRUE))
  l_hat <- (dbeta(y[(m+1):n], z$fitted[(m+1):n] * z$prec, (1 - z$fitted[(m+1):n]) * z$prec, log = TRUE))
  
  dt <- (l_tilde-l_hat)
  dt[which(dt<0)]<-0
  
  z$resid4 <- sign(y[(m+1):n]-z$fitted[(m+1):n])*sqrt(2*(dt))
  
  if(resid==1) residc <- z$resid1
  if(resid==2) residc <- z$resid2
  if(resid==3) residc <- z$resid3
  if(resid==4) residc <- z$resid4 #<- sign(y[(m+1):n]-z$fitted[(m+1):n])*sqrt(2*(dt))
  
  
  ############################################
  
  #vcov <- solve(K)
  vcov <- chol2inv(chol(K))
  z$vcov <- vcov
  
  stderror <- sqrt(diag(vcov))
  z$stderror <- stderror
  
  z$zstat <- abs(z$coef/stderror)
  z$pvalues <- 2*(1 - pnorm(z$zstat) )
  
  z$loglik <- opt$value
  #z$loglik <- opt$value*(n/(n-m))
  z$counts <- as.numeric(opt$counts[1])
  #   z$aic <- -2*z$loglik+2*(p1+q1+2)
  #   z$bic <- -2*z$loglik+log(n)*(p1+q1+2)
  
  if(any(is.na(X)==F))
  {
    z$aic <- -2*z$loglik+2*(p1+q1+2+length(beta))
    z$bic <- -2*z$loglik+log(n)*(p1+q1+2+length(beta))
  }else{
    z$aic <- -2*z$loglik+2*(p1+q1+2)
    z$bic <- -2*z$loglik+log(n)*(p1+q1+2)
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
    print(c("AIC:",round(z$aic,4)," BIC:",round(z$bic,4)),quote=F)
    
    print("Residuals:",quote=F)
    print(summary(residc))
    
    t<-seq(-5,n+6,by=1)
    
    par(mfrow=c(1,1))
    par(mar=c(2.8, 2.7, 1.2, 1)) # margens c(baixo,esq,cima,direia)
    par(mgp=c(1.7, 0.45, 0))
    plot(residc,main=" ",xlab="índices",ylab="resíduos", pch = "+",ylim=c(-4,4))
    lines(t,rep(-3,n+12),lty=2,col=1)
    lines(t,rep(3,n+12),lty=2,col=1)
    lines(t,rep(-2,n+12),lty=3,col=1)
    lines(t,rep(2,n+12),lty=3,col=1)
    
    
    max_y<- max(c(z$fitted,y),na.rm=T)
    min_y<- min(c(z$fitted,y),na.rm=T)
    plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
         xlab="valores ajustados",ylab="valores observados",
         xlim=c(0.95*min_y,max_y*1.05),
         ylim=c(0.95*min_y,max_y*1.05))
    lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
    
    densidade<-density(residc)
    plot(densidade,ylab="density",main=" ")
    lines(densidade$x,dnorm(densidade$x),lty=2)
    legend("topleft",c("Densidade estimada","Normal padrão"),#pch=vpch,
           pt.bg="white", lty=c(1,2), bty="n")
    
    acf(residc,ylab="FAC",xlab="defasagem") # função de autocorrelação
    
    pacf(residc,ylab="FACP",xlab="defasagem") # função de autocorrelação parcial
    
    max_r<- max(residc,na.rm=T)
    min_r<- min(residc,na.rm=T)
    qqnorm(residc, pch = "+",
           xlim=c(0.95*min_r,max_r*1.05),
           ylim=c(0.95*min_r,max_r*1.05),
           main="",xlab="quantis normais",ylab="quantis empíricos")
    lines(c(-10,10),c(-10,10),lty=2)
    
    par(mfrow=c(1,1))
    plot(y,type="l",ylab="serie",xlab="tempo")
    lines(z$fitted,col="red")
    
    fim<-end(y)[1]+end(y)[2]/12
    
    y_prev <- ts(y_prev, start=start(y), frequency=frequency(y))
    par(mfrow=c(1,1))
    plot(y_prev,type="l",col="red", ylim=c(min(y),max(y)),ylab="serie",xlab="tempo")
    abline(v=fim,lty=2)
    lines(y)
    
    w1<-3
    h1<-3
    
    if(diag>1)
    {
      pdf(file = "resid_v_ind.pdf",width = w1, height = h1,family = "Times")
{
        par(mfrow=c(1,1))
        par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
        par(mgp=c(1.7, 0.45, 0))
        plot(residc,main=" ",xlab="índices",ylab="resíduos", pch = "+",ylim=c(-4,4))
        lines(t,rep(-3,n+12),lty=2,col=1)
        lines(t,rep(3,n+12),lty=2,col=1)
        lines(t,rep(-2,n+12),lty=3,col=1)
        lines(t,rep(2,n+12),lty=3,col=1)
      }
dev.off()

pdf(file = "obs_v_fit.pdf",width = w1, height = h1,family = "Times")
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  plot(as.vector(z$fitted), as.vector(y), main=" ", pch = "+",
       xlab="valores ajustados",ylab="valores observados",
       xlim=c(0.95*min_y,max_y*1.05),
       ylim=c(0.95*min_y,max_y*1.05))
  lines(c(-0.2,1.2),c(-0.2,1.2),lty=2)
}
dev.off()

pdf(file = "resid_density.pdf",width = w1, height = h1,family = "Times")
{
  par(mfrow=c(1,1))
  par(mar=c(1.5, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  
  plot(densidade,ylab="densidade",main=" ",xlab=" ",ylim=c(0,1.15*max(densidade$y)))
  lines(densidade$x,dnorm(densidade$x),lty=2)
  legend("topleft",c("Densidade estimada","Normal padrão"),#pch=vpch,
         pt.bg="white", lty=c(1,2), bty="n")
}
dev.off()

pdf(file = "resid_FAC.pdf",width = w1, height = h1,family = "Times")
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  acf(residc,ylab="FAC",xlab="defasagem") # função de autocorrelação
}
dev.off()

pdf(file = "resid_FACP.pdf",width = w1, height = h1,family = "Times")
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  pacf(residc,ylab="FACP",xlab="defasagem") # função de autocorrelação parcial
}
dev.off()

pdf(file = "qq_plot.pdf",width = w1, height = h1,family = "Times")
{  
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  qqnorm(residc, pch = "+",
         xlim=c(0.95*min_r,max_r*1.05),
         ylim=c(0.95*min_r,max_r*1.05),
         main="",xlab="quantis normais",ylab="quantis empíricos")
  lines(c(-10,10),c(-10,10),lty=2)
}
dev.off()

pdf(file = "adjusted.pdf",width = 6, height = 4,family = "Times")
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  plot(y,type="l",ylab="serie",xlab="tempo")
  lines(z$fitted,col="red")
}
dev.off()

pdf(file = "forecast.pdf",width = 6, height = 4,family = "Times")
{
  par(mfrow=c(1,1))
  par(mar=c(2.8, 2.7, 1, 1)) # margens c(baixo,esq,cima,direia)
  par(mgp=c(1.7, 0.45, 0))
  plot(y_prev,type="l",col="red", ylim=c(min(y),max(y)),ylab="serie",xlab="tempo")
  abline(v=fim,lty=2)
  lines(y)
}
dev.off()
    }    
  }  
return(z)
}
