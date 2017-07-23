setwd('~/Working/Draft/gmd/')
BX    <- 1
muX   <- 5
varX  <- 1
BY    <- .2
muY   <- 4
varY  <- .8
 muZ_e  <- (BX*muX+BY*muY)/(BX+BY)
varZ_e1 <- (BX*varX+BY*varY)/(BX+BY)
varZ_e2 <- (BX**2*varX+BY**2*varY)/(BX**2+BY**2)
#A function that generates the concentrations within each size class
mu=seq(0,10,length.out=200)
dB <- function(B, mean, var, mu=seq(1,20,length.out=100),N=10000){
   M  <- length(mu)
   X  <- rnorm(N, mean=mean, sd = sqrt(var))
   NewX <- numeric(M)
   for (i in 1:M){
       #Count the number of X within one interval:
       NewX[i] <- length(which(X >= mu[i] & X < mu[i+1]))
   }
   NewX <- NewX/N*B
   return(NewX)
}
X    <-  dB(BX, muX, varX, mu) 
Y    <-  dB(BY, muY, varY, mu) 
#Calculate the resultant mean and var:
Z   <- X+Y   #Concentration of Z over the size spectrum
Z_e1 <- dB(BX+BY, muZ_e, varZ_e1, mu) 
Z_e2 <- dB(BX+BY, muZ_e, varZ_e2, mu) 

pdf('moment_transport.pdf',height=4,width=4,paper='a4')
op <- par( font.lab = 1,
             family ="serif",
             mar    = c(3.5,4,2,2),
             mgp    = c(2.3,1,0),
             mfrow  = c(1,1),
             cex.lab= 1.2)

  plot(mu,X,pch=16, type='l',ylim=c(0,.06),
       xlab=expression(paste('Ln mean size ('*µm^3*')')),
       ylab='Biomass')
points(mu,Y,pch=16, type='l',col=2)

points(mu,Z,pch=16, col=3, type='l')
#points(mu,Z_e1,col=4,type='l',lty=3)
points(mu,Z_e2,col=4,type='l',lty=3)
txt=c('Community A','Community B','A + B', 'Approximation of A+B')
legend('topright', txt, col=1:4, lty=c(1,1,1,3) )
dev.off()
#Use a normal distribution curve to fit Z:
NORM <- function(B,mean,VAR,x)  {#B: total conc.
     B*pnorm(mean,mean,VAR)*2/sqrt(2*pi*VAR)*exp(-(x-mean)**2/2/VAR)
}



