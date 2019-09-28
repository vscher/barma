
# density function
dkum<-function(y,mu,phi,a=0,b=1)
{
  d<-(1/(b-a))*((phi*log(0.5))/(log(1-mu^phi)))*y^(phi-1)*(1-y^phi)^( ((log(0.5))/(log(1-mu^phi)))-1 )
  d
}

# cumulative distribution function
pkum<-function(y,mu,phi,a=0,b=1)
{
  p<- 1- (1-y^phi)^( (log(0.5))/(log(1-mu^phi)) )
  p
}

# quantile function
qkum<-function(u,mu,phi,a=0,b=1)
{
  q<- a + (b-a)*( 1-(1-u)^( (log(1-mu^phi))/(log(0.5)) ) )^(1/phi)
  q
}

# inversion method for randon generation
rkum<-function(n,mu,phi,a=0,b=1)
{
  u<- runif(n)
  y<- a + (b-a)*( 1-(1-u)^( (log(1-mu^phi))/(log(0.5)) ) )^(1/phi)
  y
}

# pkum(0.3,0.25,2)
# qkum(0.6368364,0.25,2)
# 
# x<-seq(0,1,by=0.01)
# plot(x,dkum(x,0.5,2),type="l")
# plot(x,pkum(x,0.25,2),type="l")
# plot(qkum(x,0.25,2),x,type="l")
# plot(x,qkum(x,0.25,2),type="l")
# 
# plot(x,dkum(x,0.25,2),type="l")
# hist(rkum(1000,0.25,2),xlim=c(0,1))




# 
# w1<- 4
# h1<- 3.5
# pdf("kuma-var-phi.pdf",width = w1, height = h1,family = "Times")
# {
# par(mfrow=c(1,1))
# par(mar=c(2.8, 2.7, 0.5,0.5)) # margens c(baixo,esq,cima,direia)
# par(mgp=c(1.2, 0.5, 0))
# 
# x<-seq(0,1,by=0.01)
# plot(x,dkum(x,0.25,4),type="l",ylab=expression(f[mu](y)),xlab="y")
# lines(x,dkum(x,0.25,2))
# lines(x,dkum(x,0.25,0.8))
# lines(x,dkum(x,0.25,0.25))
# 
# text(0.29,5.5,4)
# text(0.254,3.08,2)
# text(0.233,1.66,0.8)
# text(0.22,1,0.25)
# }
# dev.off()
# 
# 
# 
# 
# w1<- 4
# h1<- 3.5
# pdf("kuma-var-mu.pdf",width = w1, height = h1,family = "Times")
# {
# par(mfrow=c(1,1))
# par(mar=c(2.8, 2.7, 0.5,0.5)) # margens c(baixo,esq,cima,direia)
# par(mgp=c(1.2, 0.5, 0))
# 
# x<-seq(0,1,by=0.01)
# plot(x,dkum(x,0.1,2),type="l",ylab=expression(f[mu](y)),xlab="y")
# lines(x,dkum(x,0.25,2))
# lines(x,dkum(x,0.5,2))
# lines(x,dkum(x,0.75,2))
# 
# text(0.14,7,0.1)
# text(0.24,3.15,0.25)
# text(0.55,1.85,0.5)
# text(0.84,2.05,0.75)
# }
# dev.off()
