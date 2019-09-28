#Simulações para dissertação de Mestrado UFPE
#Prgrama que realiza teste de Dufour and Roy robusto(utilizando os postos dos resíduos )
#Alterações em (22/05/2015) by Vinícius Scher (vinitscher@gmail.com)
#Versão (1.1)

#Informações:
#x:vetor de residuos
#lag: valor de m (defasagens da função de autocorrelação)
#fitfd: graus de liberdade que são subtraidos de x (se x é uma série de resíduos)

####################################################################################################  
Dufour.test<-function (x, lag = 1, fitdf = 0) 
{
  if (NCOL(x) > 1) 
    stop("x não é um vetor ou uma serie temporal unidimensional ")
  if (lag < 1) 
    stop("O valor do lag deve ser positivo")
  if (fitdf < 0) 
    stop("Fitdf não pode ser negativo")
  if (fitdf >= lag) 
    stop("Lag deve ser maior que os graus de liberade de Fitdf")
  METHOD <- "Dufour and Roy Robust test"
  DNAME <- deparse(substitute(x))
  aux<-rank(x) #recebe os postos das observações 
  cor <- acf(aux, lag.max = lag, type = "correlation", plot = FALSE, na.action = na.pass) #calcula autocorrelação dos postos das observações.
  obs <- cor$acf[2:(lag + 1)] #recebe as m primeiras autocorrelações excluindo a primeira que assume valor 1. 
  n <- sum(!is.na(x)) #tamanho da amostra.
  l<-lag #número de defasagens
  vec.var<-matrix(rep(0,l),nrow=l,ncol=l) #inicia a matriz de variâncias
  mu<-matrix(rep(NA,l),nrow=l,ncol=1) #inicia a matriz de médias
  
  for(k in 1:l)
  {
  mu[k,1]<-(-(n-k)/(n*(n-1))) #calcula as m médias  
  var1<-((5*n^4)-((5*k+9)*n^3)+(9*n^2*(k-2))+(2*k*n*(5*k+8))+(16*k^2))
  var2<-(5*n^2*(n-1)^2*(n+1))
  vec.var[k,k]<-(var1/var2) #calcula as m variâncias 
  }
 
  A<-t((obs-t(mu)))
  QDR<-t(A) %*% solve(vec.var) %*% A #realiza o produto matricial para a estatística do teste
  STATISTIC <- as.numeric(QDR) #estatística de teste
  names(STATISTIC) <- "X-squared on Rank of Residuals for fitted ARMA process"
  mydf <- lag - fitdf #graus de liberdade
  PARAMETER <- c(mydf)
  names(PARAMETER) <- c("df")
  PVAL <- 1 - pchisq(STATISTIC, df = mydf) #p-valor
  names(PVAL) <- "p-value"
  structure(list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, method = METHOD, data.name = DNAME), class = "htest")
  
}