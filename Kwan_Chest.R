#Simulações para dissertação de Mestrado UFPE
#Programa que realiza os 4 testes de Kwan-Sim com e sem autocorrelação parcial 
#Alterações em (05/06/2016) by Vinícius Scher (vinitscher@gmail.com)
#Versão (1.01)

#Informações:
#x: vetor de residuos
#lag: valor de m (defasagens da função de autocorrelação)
#fitfd: graus de liberdade que são subtraidos de x (se x é uma série de resíduos)
#type: tipo de correlação parcial="partial" ou normal="correlation"(default)
#test: qual tipo de teste de Kwan-Sim será realizado 1, 2, 3, 4.

####################################################################################################  
Kwan.sim.chest<-function (x, lag = 1, fitdf = 0, type="correlation", test=1) 
{
  if (NCOL(x) > 1) 
    stop("x não é um vetor ou uma serie temporal unidimensional ")
  if (lag < 1) 
    stop("O valor do lag deve ser positivo")
  if (fitdf < 0) 
    stop("Fitdf não pode ser negativo")
  if (fitdf >= lag) 
    stop("Lag deve ser maior que os graus de liberade de Fitdf")
  if (type == "partial")
  {
    METHOD <- "Teste de kwan and Sim com autocorrelação parcial"
    cor <- acf(x, lag.max = lag, type = "partial", plot = FALSE, na.action = na.pass)
    obs <- cor$acf[1:lag]
    name<-"X-squared on Residuals(ACF - Partial) for fitted ARMA process"
  }
   if(type == "correlation") 
     {
     METHOD <- "Teste de kwan and Sim com autocorrelação normal"
     cor <- acf(x, lag.max = lag, type = "correlation", plot = FALSE, na.action = na.pass)
     obs <- cor$acf[2:(lag + 1)]
     name<-"X-squared on Residuals(ACF - Normal) for fitted ARMA process "
   } 
  
  if(test<0 || test>4 )
    stop("O valor para o teste a ser realizado deve ser 1, 2, 3 ou 4.")
  DNAME <- deparse(substitute(x))
  n <- sum(!is.na(x))
  
  
  
  if(test == 1)
  {
    z1<-0.5*(log((1+obs)/(1-obs)))
    Q<-sum((seq.int(n - 4,n - lag - 3) * z1^2))
    result<-0
    for(k in 1:lag)
    {
      aux1<-0
      aux2<-0
      N<-0
      D<-0
      aux1<-(n-k)/(n*(n+2))
      aux2<-((3*((n^2)-((2*k-6)*n)+(k-10)))/(n*(n+2)*(n+4)*(n+6)))
      N<-(n-k-3)
      a1<-1
      result<-result+(N*((a1^2)*aux1+((2/3)*(a1)*aux2)))+(1/n^2)
    }
    EQ<-result
    METHOD <- "Kwan and Sim(1) test"
  }
  
  if(test == 2)
  {
    z1<-0.5*(log((1+obs)/(1-obs)))
    z2<-z1-(0.25*((seq.int(n-1,n-lag)^-1)*(3*z1+obs)))
    Q<-sum((seq.int(n - 2,n - lag - 1) * z2^2))
    result<-0
    for(k in 1:lag)
    {
      aux1<-0
      aux2<-0
      N<-0
      D<-0
      aux1<-(n-k)/(n*(n+2))
      aux2<-((3*((n^2)-((2*k-6)*n)+(k-10)))/(n*(n+2)*(n+4)*(n+6)))
      N<-(n-k-1)
      a2<-1-(n-k)^-1
      result<-result+(N*((a2^2)*aux1+((2/3)*(a2)*aux2)))+(1/n^2)
    }
    EQ<-result
    METHOD <- "Kwan and Sim(2) test"
  }
  
  if(test == 3)
  {
    z1<-0.5*(log((1+obs)/(1-obs)))
    z2<-z1-(0.25*((seq.int(n-1,n-lag)^-1)*(3*z1+obs)))
    z3<-z2-((1/96)*(seq.int(n-1,n-lag)^-2)*(23*z1+33*obs-5*(obs^3)))
    Q<-sum((seq.int(n - 2,n - lag - 1) * z3^2))
    
    result<-0
    for(k in 1:lag)
    {
      aux1<-0
      aux2<-0
      N<-0
      D<-0
      aux1<-(n-k)/(n*(n+2))
      aux2<-((3*((n^2)-((2*k-6)*n)+(k-10)))/(n*(n+2)*(n+4)*(n+6)))
      N<-(n-k-1)
      a2<-1-(n-k)^-1
      a3<-a2-((7/12)*(n-k)^-2)
      result<-result+(N*((a3^2)*aux1+((2/3)*(a3)*aux2)))+(1/n^2)
    }
    EQ<-result
    METHOD <- "Kwan and Sim(3) test"
  }
    
    if(test ==4)
    {
      z4<-asin(obs)
      Q<-sum(((seq.int(n-1,n-lag)^2)/(seq.int(n - 2,n - lag -1  ))) * z4^2) #Q*
      
      result<-0
      for(k in 1:lag)
      {
        aux1<-0
        aux2<-0
        N<-0
        D<-0
        aux1<-(n-k)/(n*(n+2))
        aux2<-((3*((n^2)-((2*k-6)*n)+(k-10)))/(n*(n+2)*(n+4)*(n+6)))
        N<-(n-k)^2
        D<-(n-k-1)
        result<-result+((N/D)*(aux1+(aux2/3)))
      }
      EQ<-result
      METHOD <- "Kwan and Sim(4) test"
    }
  
  STATISTIC <- Q
  names(STATISTIC) <- name
  mydf <- EQ - fitdf
  PARAMETER <- c(mydf)
  names(PARAMETER) <- c("df")
  PVAL <- 1 - pchisq(STATISTIC, df = mydf)
  names(PVAL) <- "p-value"
  structure(list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, method = METHOD,  data.name = DNAME), class = "htest")
}




