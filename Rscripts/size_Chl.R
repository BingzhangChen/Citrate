#Check how different for the Chl distribution against phytoplankton biomass distribution:
N <- 200 #Number of phytoplankton size classes
setwd('~/Working/FlexEFT1D/Rscripts')
source('~/Working/FlexEFT1D/Rscripts/mu_EFT.R')

#Generate a phytoplankton histogram:
logV  <- function(ESD)log(pi/6*ESD**3)
PMU   <- logV(5)
VAR   <- .1
mu    <-  seq(-3,15,length.out=200)
P     <-  dB(B = 1, PMU, VAR, mu)

plot(mu,P,pch=16, type='l',ylim=c(0,.06),
     xlab=expression(paste('Ln mean size ('*µm^3*')')),
     ylab='Biomass')

#Calculate Chl = P*theta/QN
ESD <- c(0.5, 2,  5,  20, 100)
PMU_<- log(pi/6*ESD**3)

#Assume environmental condition:
plot_theta = function(NO3, PAR, Temp=15, aI0=0.05,alphaI=0.08,mu0=2,alphamu=.2,
                      betamu=-.017,KN0=.2,alphaK=.27){
   #Calculate theta:
   Theta  =numeric(length(mu))
   QN     =Theta
   for (i in 1:length(mu)){
       L   =mu[i]
       rmax=mu0*exp(L*alphamu+L**2*betamu)
       KN  =ScaleTrait(L,KN0,alphaK)
       aI  =ScaleTrait(L,aI0,alphaI)
      Theta[i]=mu_fix(NO3=NO3,par_=PAR, mu0=rmax,KN0=KN,
                     aI0C=aI, thetamax=.5)$Theta
      QN[i]=mu_fix(NO3=NO3,par_=PAR, mu0=rmax,KN0=KN,
                     aI0C=aI, thetamax=.5)$QN
   
   }
   XLAB <- 'ESD (µm)'
   plot(mu, Theta,      xaxt='n',
        xlab=XLAB,
        ylab=expression(paste('Chl:C (µg '*mol^-1*')'))) 
   axis(1, at=PMU_, labels=as.character(ESD))
   txt=paste0('NO3 = ',NO3,' PAR = ',PAR)
   mtext(txt)
   LM=lm(Theta~mu)
   abline(LM) 
   text(2,min(Theta)+.002,paste0('Slope = ',round(as.numeric(coef(LM)[2]),3)))
}

pdf('Scaling_Theta.pdf', width = 6, height = 8)
op <- par(font.lab = 1,
            family ="serif", cex.axis=1.2, cex.lab=1.2,
            mar    = c(4,3.8,.5,.5),
            mgp    = c(2.3,1,0),
            mfcol  = c(2,1),
            oma    = c(4,4,1,0)) 

#Surface:
plot_theta(NO3=.1, PAR=300)

#Deep
plot_theta(NO3=10, PAR=10)


dev.off()
plot(mu, QN,    type='l')
