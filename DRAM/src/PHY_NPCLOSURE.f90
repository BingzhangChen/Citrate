! The subroutine for phytoplankton closure model
SUBROUTINE PHY_NPCLOSURE(NO3,PAR_,Temp_,PHY,VPHY,VNO3, COVNP,mu0hatSI, muNet,SI,theta,QN, Snp, Chl, NPP, SVPHY, SVNO3, SCOVNP)
USE PARAM_MOD
Implicit none
real, intent(in)  :: NO3, PAR_,Temp_,PHY, VPHY, VNO3, COVNP 
real, intent(out) :: muNet, theta, QN, SI, Snp, Chl, NPP, SVPHY, SVNO3, &
SCOVNP,mu0hatSI

! muNet: mean growth rate at <N>
real :: mu0hat
real :: alphaI, cff,cff1, eta, dYdN,d2YdN2, dEta_dN, d2ChldN2, Q
real :: KN, tf, fN, dmuQ_dN, d2NPPdN2
real, parameter :: Qmin = 0.06, Qmax=0.18

alphaI   = exp(params(iaI0_C))
tf       = TEMPBOL(Ep,Temp_)   !Temperature effect

!The temperature and light dependent component
mu0hat   = tf*exp(params(imu0))
mu0hatSI = mu_Edwards2015(PAR_, params(iIopt),mu0hat, alphaI) 
SI       = mu0hatSI/mu0hat
KN       = exp(params(iKN))
fN       = NO3/(NO3 + Kn)  !Nitrogen limitation index

! Phytoplankton growth rate at the mean nitrogen:
muNet = mu0hatSI*fN

! Snp: ensemble mean production (nitrogen based)
Snp   = PHY*muNet + mu0hatSI*(Kn*COVNP/(Kn+NO3)**2 - Kn*PHY*VNO3/(Kn+NO3)**3)

!N:C ratio at <N>
cff1  = 1.-Qmin/Qmax
cff   = 1.-cff1*fN
QN    = Qmin/cff
Q     = 1./QN

!Chl:C ratio at <N>
cff   = (thetamax - thetamin)/PAR_
theta = thetamin + muNet/alphaI*cff   !Unit: gChl/molC

!Chl:N ratio at <N>
eta      = theta*Q
dYdN     = Kn/(NO3 + Kn)**2
d2YdN2   = -dYdN * 2. /(NO3 + Kn) 
cff1     = 1./Qmax - 1./Qmin
dEta_dN  = dYdN*(theta*cff1 + Q*mu0hatSI/alphaI*cff)
d2ChldN2 = 2.*PHY*(dYdN**2 * (cff1*cff*mu0hatSI/alphaI) - dEta_dN/(NO3+Kn))

!Ensemble mean Chl
Chl      = PHY*eta + .5*(2.*COVNP*dEta_dN + VNO3*d2ChldN2)
Chl      = max(Chl, eps)

dmuQ_dN  = dYdN*(muNet*cff1 + Q*mu0hatSI)
d2NPPdN2 = PHY*(d2YdN2*dmuQ_dN/dYdN + 2.*dYdN**2*mu0hatSI*cff1)

! NPP: ensemble mean carbon based primary production
NPP      = PHY*muNet*Q + .5*(2.*COVNP*dmuQ_dN + VNO3*d2NPPdN2)
NPP      = max(NPP, 0.)

! Calculate sources and sinks of variances of N, P, and covariance of NP
SVPHY    = 2.*mu0hatSI*(fN*VPHY         + Kn*PHY*       COVNP/(Kn+NO3)**2)
SVNO3    =-2.*mu0hatSI*(fN*COVNP        + Kn*PHY*        VNO3/(Kn+NO3)**2) 
SCOVNP   =    mu0hatSI*(fN*(COVNP-VPHY) + Kn*PHY*(VNO3-COVNP)/(Kn+NO3)**2) 
return
END SUBROUTINE PHY_NPCLOSURE
