#Simulações para dissertação de Mestrado UFPE
#Prgrama que realiza teste de Monti
#Alterações em (20/05/2015) by Vinícius Scher (vinitscher@gmail.com)
#Versão (1.1)

#Informações:
#x:vetor de residuos
#lag: valor de m (defasagens da função de autocorrelação)
#fitfd: graus de liberdade que são subtraidos de x (se x é uma série de resíduos)

####################################################################################################  
Monti.test<-function (x, lag = 1, fitdf = 0) 
{
  if (NCOL(x) > 1) 
    stop("x não é um vetor ou uma serie temporal unidimensional ")
  if (lag < 1) 
    stop("O valor do lag deve ser positivo")
  if (fitdf < 0) 
    stop("Fitdf não pode ser negativo")
  if (fitdf >= lag) 
    stop("Lag deve ser maior que os graus de liberade de Fitdf")
  METHOD <- "Monti test"
  DNAME <- deparse(substitute(x))
  cor <- acf(x, lag.max = lag, type = "partial", plot = FALSE, na.action = na.pass)
  obs <- cor$acf[1:lag]
  n <- sum(!is.na(x))
  STATISTIC <- n * (n + 2) * sum((1/seq.int(n - 1,n - lag) * obs^2))
  names(STATISTIC) <- "X-squared on Residuals for fitted ARMA process"
  mydf <- lag - fitdf
  PARAMETER <- c(mydf)
  names(PARAMETER) <- c("df")
  PVAL <- 1 - pchisq(STATISTIC, df = mydf)
  PVAL
  names(PVAL) <- "p-value"
  structure(list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, method = METHOD, data.name = DNAME), class = "htest")
  
}