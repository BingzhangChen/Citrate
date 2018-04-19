#Validate derivative functions for NP closure model
SUBROUTINE PHY_NPCLOSURE(NO3,PAR_,Temp_,PHY,VPHY,VNO3, COVNP, muNet,SI,theta,QN, Snp, Chl, NPP, SVPHY, SVNO3, SCOVNP)
USE bio_MOD, ONLY : TEMPBOL, params, mu_Edwards2015 

 mu_Edwards2015 <- function(I_, lnIopt, mumax, alpha) {
  return(1./( I_/(alpha*exp(2.*lnIopt)) + 1./mumax -2./(alpha*exp(lnIopt))  + 1./(I_*alpha)  ) 
}

NPCLOS <- function(NO3, PAR_, Temp_,PHY,VPHY,VNO3, COVNP,
                   Iopt = log(1e3),mu0hat = 1.0,
                 alphaI =  0.055, KN = 0.5,
                   Snp, Chl, NPP, SVPHY, SVNO3, SCOVNP){

  Qmin   = 0.06; Qmax=0.18
  thetamin = 0.02
  thetamax = 0.63
  #The temperature and light dependent component
  SI       = mu_Edwards2015(PAR_, Iopt,mu0hat, alphaI) 
  mu0hatSI = mu0hat*SI
  fN = NO3/(NO3 + Kn)  #Nitrogen limitation index

  # Phytoplankton growth rate at the mean nitrogen:
  muNet = mu0hatSI*fN

  # Snp: ensemble mean production (nitrogen based)
  Snp   = PHY*muNet +
      mu0hatSI*(Kn*COVNP/(Kn+NO3)**2 - Kn*PHY*VNO3/(Kn+NO3)**3)

  #N:C ratio at <N>
  cff1  = 1-Qmin/Qmax
  cff   = 1-cff1*fN
  QN    = Qmin/cff
  Q     = 1./QN
  #Chl:C ratio at <N>
  cff   = (thetamax - thetamin)/PAR_
  theta = thetamin+muNet/alphaI*cff   #Unit: gChl/molC

  #Chl:N ratio at <N>
  eta   = theta/QN
  dYdN  = Kn/(NO3 + Kn)**2
  d2YdN2= -2.*dYdN/(NO3 + Kn)
  cff1  = 1./Qmax - 1./Qmin
  dEta_dN  = dYdN*(theta*cff1 + Q*mu0hatSI/alphaI*cff)
  d2ChldN2 = 2.*PHY*(dYdN**2*(cff1*cff*mu0hatSI/alphaI)-dEta_dN/(NO3+Kn))

 #Ensemble mean Chl
 Chl = PHY*eta + .5*(2.*COVNP*dEta_dN + VNO3*d2ChldN2)
 
 dmuQ_dN = dYdN*(muNet*cff1 + Q*mu0hatSI)
 d2NPPdN2= PHY*(d2YdN2*dmuQ_dN/dYdN + 2.*dYdN**2*mu0hatSI *cff1)
 # NPP: ensemble mean carbon based primary production
 NPP = PHY*muNet*Q + .5*(2.*COVNP*dmuQ_dN + VNO3*d2NPPdN2)

  # Calculate sources and sinks of variances of N, P, and covariance of NP
  SVPHY = 2.*mu0hatSI*(fN*VPHY         + Kn*PHY*       COVNP/(Kn+NO3)**2)
  SVNO3 =-2.*mu0hatSI*(fN*COVNP        + Kn*PHY*        VNO3/(Kn+NO3)**2) 
  SCOVNP=    mu0hatSI*(fN*(COVNP-VPHY) + Kn*PHY*(VNO3-COVNP)/(Kn+NO3)**2) 
}

#Derive derivatives using R

eta = expression((thetamin+mu0*N/(N+Kn)/(alpha*PAR)*(thetamax - thetamin)) * (1-(1-Qmin/Qmax)*N/(N+Kn))/Qmin)
detadN = D(eta, 'N')
d2etadN2 = D(detadN, 'N')
eval(eta, list(thetamin = 0.02, thetamax = 0.63, mu0 = 1, Kn=0.5, N = 0.1, alpha=0.055, PAR=1000, Qmin=0.06, Qmax=0.18))
eval(detadN, list(thetamin = 0.02, thetamax = 0.63, mu0 = 1, Kn=0.5, N = 0.1, alpha=0.055, PAR=1000, Qmin=0.06, Qmax=0.18))
eval(d2etadN2, list(thetamin = 0.02, thetamax = 0.63, mu0 = 1, Kn=0.5, N = 0.1, alpha=0.055, PAR=1000, Qmin=0.06, Qmax=0.18))


muQ = expression( mu0*N/(N+Kn)*(1-(1-Qmin/Qmax)*N/(N+Kn))/Qmin)
dmuQdN = D(muQ, 'N')
d2muQdN2 = D(dmuQdN, 'N')
eval(muQ, list(mu0 = 1, Kn=0.5, N = 0.1, Qmin=0.06, Qmax=0.18))
eval(dmuQdN, list(mu0 = 1, Kn=0.5, N = 0.1, Qmin=0.06, Qmax=0.18))
eval(d2muQdN2, list(mu0 = 1, Kn=0.5, N = 0.1, Qmin=0.06, Qmax=0.18))
