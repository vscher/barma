#Simulações para dissertação de Mestrado UFPE
#Prgrama que realiza teste de Ljung-Box
#Alterações em (19/05/2015) by Vinícius Scher (vinitscher@gmail.com)
#Versão (1.1)

#Informações:
#x:vetor de residuos
#lag: valor de m (defasagens da função de autocorrelação)
#fitfd: graus de liberdade que são subtraidos de x (se x é uma série de resíduos)

####################################################################################################  
Ljung.Box<-function (x, lag = 1, fitdf = 0) 
{
  if (NCOL(x) > 1) 
    stop("x não é um vetor ou uma serie temporal unidimensional ")
  if (lag < 1) 
    stop("O valor do lag deve ser positivo")
  if (fitdf < 0) 
    stop("Fitdf não pode ser negativo")
  if (fitdf >= lag) 
    stop("Lag deve ser maior que os graus de liberade de Fitdf")
    METHOD <- "Ljung-Box test"
    DNAME <- deparse(substitute(x))
    cor <- acf(x, lag.max = lag, type = "correlation", plot = FALSE, na.action = na.pass) #calcula a função de autocorrelação.
    obs <- cor$acf[2:(lag + 1)] #recebe as m autocorrelações,exclui a primeira autocorrelação pois assume valor 1.
    n <- sum(!is.na(x)) #tamanho da amostra.
    
    STATISTIC <- n * (n + 2) * sum((1/seq.int(n - 1,n - lag) * obs^2)) #estatística de teste.
    names(STATISTIC) <- "X-squared on Residuals for fitted ARMA process"
    mydf <- lag - fitdf #calcula os graus de liberdade
    PARAMETER <- c(mydf) #graus de liberdade
    names(PARAMETER) <- c("df")
    PVAL <- 1 - pchisq(STATISTIC, df = mydf) #calcula o p-valor
    names(PVAL) <- "p-value"
    structure(list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, method = METHOD, data.name = DNAME), class = "htest")
    
}

