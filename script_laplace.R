xrange<-seq(2,4,0.001)


approx <- function(x) sqrt((2*pi)/x)*(x/exp(1))^x
# ap2 <- function(x) exp(x*log(x))*exp(-x)*sqrt((2*pi)/x)

# plot(ap2(xrange))


par(mfrow=c(1,1))
  par(mar=c(2.6, 3.3, 0.5, 0.5)) # margens c(baixo,esq,cima,direia)
par(mgp=c(1.5, 0.45, 0))


plot(xrange,gamma(xrange),type="l",cex.lab=1.3,lty=2,lwd=1.5,ylab=expression(paste(sqrt(2*pi/x)," ",(x/exp)^x)),xlab=expression(x))
lines(xrange,approx(xrange))
legend("topleft",c("Exata","Aprox. Laplace"),lty=c(1,2),bty="n",lwd=1.8,cex=1.4)





integrate e^(-t) t^(x-1) dt  from 0 to inf