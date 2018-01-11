setwd('~/Working/Draft/gmd/')
source('~/Working/FlexEFT1D/Rscripts/mu_EFT.R')
source('~/Working/FlexEFT1D/Rscripts/plot_1D.R')
BX    <- 1
muX   <- 5
varX  <- 1

BY    <- 1
muY   <- 1
varY  <- 1

ADD_mom <- function(BX, muX, varX, BY, muY, varY, legend=F){
    muZ  <- (BX*muX+BY*muY)/(BX+BY)
   varZ  <- (BX*(varX+muX^2)+BY*(varY+muY^2))/(BX+BY)-muZ^2
   mu    <-  seq(0,10,length.out=200)
   X     <-  dB(BX, muX, varX, mu) 
   Y     <-  dB(BY, muY, varY, mu) 
   #Calculate the resultant mean and var:
   Z   <- X+Y   #Concentration of Z over the size spectrum
   Z1  <- dB(BX+BY, muZ, varZ, mu) 
   plot(mu,X,pch=16, type='l',ylim=c(0,.04),
       xlab=expression(paste('Ln size ('*µm^3*')')),
       ylab='Biomass')
   points(mu,Y,pch=16, type='l',col=2)
   points(mu,Z,pch=16, col=3, type='l')
   points(mu,Z1,col=4,type='l',lty=3)
   text(7,.02,  bquote(italic(P)[A] ~ ' = '~.(BX)),  pos=4)
   text(7,.017, bquote(bar(italic(l))[A] ~ ' = '~.(muX)),  pos=4)
   text(7,.014, bquote(italic(v)[A] ~ ' = '~.(varX)),pos=4)
   text(7,.011, bquote(italic(P)[B] ~ ' = '~.(BY)),  pos=4)
   text(7,.008, bquote(bar(italic(l))[B] ~ ' = '~.(muY)),  pos=4)
   text(7,.005, bquote(italic(v)[B] ~ ' = '~.(varY)),pos=4)
   if (legend) {
      txt=c('Community A','Community B','A + B',
            'Approximation of A+B')
      legend('topright', txt, col=1:4, lty=c(1,1,1,3) )
   }
}

model <- 'NPZDcont'
Stn   <- 'K2'
filedir  <- paste0('~/working/FlexEFT1D/DRAM/',model,'/BOTH_TD/')
PHY      <- get_Final_Data(filedir, Stn, 'PHY_T')
Y        <- as.numeric(PHY$data[60,])

PMU      <- get_Final_Data(filedir, Stn, 'R_PMU')
P1       <- as.numeric(PMU$data[60,])

VAR      <- get_Final_Data(filedir, Stn, 'R_VAR')
V1       <- as.numeric(VAR$data[60,])

Aks      <- get_Final_Data(filedir, Stn, 'Aks')
A1       <- as.numeric(Aks$data[60,])


pdf('moment_transport.pdf',height=6,width=6,paper='a4')
op <- par( font.lab = 1,
             family ="serif",
             mar    = c(3.5,4,1,.2),
             mgp    = c(2.3,1,0),
             mfrow  = c(2,2),
             cex.lab= 1.2)

 #Mixing between two communities with equal biomass but very different mean size
  ADD_mom(BX=.6,muX=5,varX=1,BY=.6,muY=3,varY=1,legend=T)
  ADD_mom(BX=.6,muX=5,varX=.0001,BY=.6,muY=3,varY=1,legend=T)
  mtext('a)',adj=0, outer=F)

  ADD_mom(BX=.6,muX=5,varX=1,BY=.6,muY=2,varY=1)
  mtext('b)',adj=0, outer=F)
  ADD_mom(BX=.6,muX=5,varX=1,BY=.1,muY=3,varY=1)
  mtext('c)',adj=0, outer=F)
  ADD_mom(BX=.6,muX=5,varX=1,BY=.6,muY=3,varY=5)
  mtext('d)',adj=0, outer=F)

 #Plot a vertical profile of B, mu and Var at the end of Feb. from K2:
  #Varname  <- expression(paste("PHY, "*bar(italic(l))*', or '*italic(v)^2*' (Relative units)'))
  #plot(Y, PHY$depth, type='b', xlab=Varname, ylab='Depth (m)')
  #points(P1/10,  PHY$depth, type='b', col=2)
  #points(V1/100, PHY$depth, type='b', col=3)
  #legend('bottomright',legend=c('PHY',expression(bar(italic(l))),expression(italic(v)^2)), 
  #       lty=1,pch=1,col=1:3)
  #mtext('c)',adj=0, outer=F)
dev.off()

#Use a normal distribution curve to fit Z:
NORM <- function(B,mean,VAR,x)  {#B: total conc.
     B*pnorm(mean,mean,VAR)*2/sqrt(2*pi*VAR)*exp(-(x-mean)**2/2/VAR)
}



