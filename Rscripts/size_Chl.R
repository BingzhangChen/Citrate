#Check how different for the Chl distribution against phytoplankton biomass distribution:
N <- 200 #Number of phytoplankton size classes

source('~/Working/FlexEFT1D/Rscripts/mu_EFT.R')

#Generate a phytoplankton histogram:
logV  <- function(ESD)log(pi/6*ESD**3)
PMU   <- logV(5)
VAR   <- .1
mu    <-  seq(0,10,length.out=200)
P     <-  dB(B = 1, PMU, VAR, mu)

plot(mu,P,pch=16, type='l',ylim=c(0,.06),
     xlab=expression(paste('Ln mean size ('*µm^3*')')),
     ylab='Biomass')

#Calculate Chl = P*theta/QN

#Assume environmental condition:
NO3 <- .2
PAR <- 100
Temp <- 15

#Calculate theta:



