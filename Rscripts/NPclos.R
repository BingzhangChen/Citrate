#Validate derivative functions for NP closure model
SUBROUTINE PHY_NPCLOSURE(NO3,PAR_,Temp_,PHY,VPHY,VNO3, COVNP, muNet,SI,theta,QN, Snp, Chl, NPP, SVPHY, SVNO3, SCOVNP)
USE bio_MOD, ONLY : TEMPBOL, params, mu_Edwards2015 
USE bio_MOD, ONLY : iIopt, imu0, iaI0_C, iKN, Ep
USE bio_MOD, ONLY : thetamax, thetamin

NPCLOS <- function(NO3, PAR_, Temp_,PHY,VPHY,VNO3, COVNP, Snp, Chl, NPP, SVPHY, SVNO3, SCOVNP){

  Qmin = 0.06; Qmax=0.18
  alphaI =  0.055

#The temperature and light dependent component
 mu0hat  = 1.0
SI       = mu_Edwards2015(PAR_, params(iIopt),mu0hat, alphaI) 
mu0hatSI = mu0hat*SI

KN = exp(params(iKN))
fN = NO3/(NO3 + Kn)  #Nitrogen limitation index

# Phytoplankton growth rate at the mean nitrogen:
muNet = mu0hatSI*fN

# Snp: ensemble mean production (nitrogen based)
Snp   = PHY*muNet + mu0hatSI*(Kn*COVNP/(Kn+NO3)**2 - Kn*PHY*VNO3/(Kn+NO3)**3)

#N:C ratio at <N>
cff1  = 1d0-Qmin/Qmax
cff   = 1d0-cff1*fN
QN    = Qmin/cff
Q     = 1./QN
#Chl:C ratio at <N>
cff   = (thetamax - thetamin)/PAR_
theta = thetamin+muNet/alphaI*cff   #Unit: gChl/molC

#Chl:N ratio at <N>
eta   = theta/QN
dYdN  = Kn/(NO3 + Kn)**2
cff1  = 1./Qmax - 1./Qmin
dEta_dN  = dYdN*(theta*cff1 + Q*mu0hatSI/alphaI*cff)
d2ChldN2 = -2.*PHY*Kn/(NO3 + Kn)**3 * (thetamin*cff1 + Q*mu0hatSI/alphaI*cff)

#Ensemble mean Chl
Chl = PHY*eta + .5*(2.*COVNP*dEta_dN + VNO3*d2ChldN2)

dmuQ_dN = dYdN*(muNet*cff1 + Q*mu0hatSI)
d2NPPdN2= 2.*PHY*mu0hatSI*Kn/(NO3+Kn)**3 * (cff1*(Kn-2.*NO3)/(NO3+Kn) + 1./Qmin)
# NPP: ensemble mean carbon based primary production
NPP = PHY*muNet*Q + .5*(2.*COVNP*dmuQ_dN + VNO3*d2NPPdN2)

# Calculate sources and sinks of variances of N, P, and covariance of NP
SVPHY = 2.*mu0hatSI*(fN*VPHY         + Kn*PHY*       COVNP/(Kn+NO3)**2)
SVNO3 =-2.*mu0hatSI*(fN*COVNP        + Kn*PHY*        VNO3/(Kn+NO3)**2) 
SCOVNP=    mu0hatSI*(fN*(COVNP-VPHY) + Kn*PHY*(VNO3-COVNP)/(Kn+NO3)**2) 
return
END SUBROUTINE PHY_NPCLOSURE

