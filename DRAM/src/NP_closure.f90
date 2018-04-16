SUBROUTINE NP_CLOSURE
    ! This NP_closure model has 5 tracers: <N>, <P>, <N'>^2, <P'>^2, <N'P'>
    ! Governing functions follow Mandal et al. JPR (2016)
USE bio_MOD
USE MOD_1D, only: it, nsave
IMPLICIT NONE
integer :: k,j,i
!INPUT PARAMETERS:
real :: tC,par_
!LOCAL VARIABLES of phytoplankton:
real :: NO3,PHY, VPHY, VNO3, COVNP, SVNO3, SVPHY, SCOVNP
real :: QN  ! cell quota related variables
real :: muNet, Dp, PP_PN 
real :: SI
real :: theta, Snp
real :: Chl, NPPt
!-----------------------------------------------------------------------
DO k = nlev, 1, -1   
   ! Retrieve current (local) state variable values.
   tC      = Temp(k)
   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif
   Varout(oPAR_,k) = par_
   NO3    = Vars(iNO3,    k)
   PHY    = Vars(iPHY(1), k)
   VPHY   = Vars(iVPHY,   k)
   VNO3   = Vars(iVNO3,   k)
   COVNP  = Vars(iCOVNP,  k)

   IF (mod(it, nsave) .EQ. 1) THEN
     !Use nonmixed light to calculate NPP
     call PHY_NPCLOSURE(NO3,PAR_,tC,PHY,VPHY,VNO3, COVNP, muNet,SI,theta,QN, &
                        Snp, Chl, NPPt, SVPHY, SVNO3, SCOVNP)
     Varout(oPPt, k)  = NPPt * 12d0 ! Correct the unit to ug C/L
   ENDIF
   !Use mixed light to estimate source and sink
   call PHY_NPCLOSURE(NO3,PAR_,tC,PHY,VPHY,VNO3, COVNP, muNet,SI,theta,QN, &
                        Snp, Chl, NPPt, SVPHY, SVNO3, SCOVNP)
   ! Save some model outputs:
   Varout(oTheta(1),k)= theta! Chl:C ratio at <N>
   Varout(oQN(1)   ,k)= QN   ! N:C ratio at <N> 
   Varout(oSI(1),   k)= SI   ! Light limitation
   Varout(oCHLt,k)    = Chl  ! ensemble mean Chl
!=============================================================
!! Solve ODE functions:
!All rates have been multiplied by dtdays to get the real rate correponding to the actual time step
   Dp   = params(iDp)  ! Mortality rate of phytoplankton
   PP_PN= Snp - PHY*Dp

!Update tracers:
   NO3  = NO3   - PP_PN*dtdays
   PHY  = PHY   + PP_PN*dtdays
   VPHY = VPHY  + (SVPHY-2.*Dp*VPHY )*dtdays
   VNO3 = VNO3  + (SVNO3+2.*Dp*COVNP)*dtdays
   COVNP= COVNP + (SCOVNP+Dp*(VPHY-COVNP))*dtdays
   Varout(oNO3,k)      = NO3
   Varout(oPHY(1),k)   = PHY
   Varout(oVPHY  ,k)   = VPHY
   Varout(oVNO3  ,k)   = VNO3
   Varout(oCOVNP ,k)   = COVNP
   Varout(oPHYt,  k)   = PHY
   Varout(omuNet(1),k) = muNet               !Growth rate at <N>
ENDDO
END SUBROUTINE NP_CLOSURE 

! The subroutine only for phytoplankton
SUBROUTINE PHY_NPCLOSURE(NO3,PAR_,Temp_,PHY,VPHY,VNO3, COVNP, muNet,SI,theta,QN, Snp, Chl, NPP, SVPHY, SVNO3, SCOVNP)
USE bio_MOD, ONLY : TEMPBOL, params, mu_Edwards2015 
USE bio_MOD, ONLY : iIopt, imu0, iaI0_C, iKN, Ep
USE bio_MOD, ONLY : thetamax, thetamin
implicit none
real, intent(in)  :: NO3, PAR_,Temp_,PHY, VPHY, VNO3, COVNP 
real, intent(out) :: muNet, theta, QN, SI, Snp, Chl, NPP, SVPHY, SVNO3, SCOVNP

! muNet: mean growth rate at <N>
real :: mu0hat, mu0hatSI
real :: alphaI, cff,cff1, eta, dYdN, dEta_dN, d2ChldN2, Q
real :: KN, tf, fN, dmuQ_dN, d2NPPdN2
real, parameter :: Qmin = 0.06, Qmax=0.18

alphaI   = params(iaI0_C)
tf       = TEMPBOL(Ep,Temp_)   !Temperature effect

!The temperature and light dependent component
mu0hat   = tf*params(imu0)
SI       = mu_Edwards2015(PAR_, params(iIopt),mu0hat, alphaI) 
mu0hatSI = mu0hat*SI

KN = params(iKN)
fN = NO3/(NO3 + Kn)  !Nitrogen limitation index

! Phytoplankton growth rate at the mean nitrogen:
muNet = mu0hatSI*fN

! Snp: ensemble mean production (nitrogen based)
Snp   = PHY*muNet + mu0hatSI*(Kn*COVNP/(Kn+NO3)**2 - Kn*PHY*VNO3/(Kn+NO3)**3)

!N:C ratio at <N>
cff1  = 1d0-Qmin/Qmax
cff   = 1d0-cff1*fN
QN    = Qmin/cff
Q     = 1./QN
!Chl:C ratio at <N>
cff   = (thetamax - thetamin)/PAR_
theta = thetamin+muNet/alphaI*cff   !Unit: gChl/molC

!Chl:N ratio at <N>
eta   = theta/QN
dYdN  = Kn/(NO3 + Kn)**2
cff1  = 1./Qmax - 1./Qmin
dEta_dN  = dYdN*(theta*cff1 + Q*mu0hatSI/alphaI*cff)
d2ChldN2 = -2.*PHY*Kn/(NO3 + Kn)**3 * (thetamin*cff1 + Q*mu0hatSI/alphaI*cff)

!Ensemble mean Chl
Chl = PHY*eta + .5*(2.*COVNP*dEta_dN + VNO3*d2ChldN2)

dmuQ_dN = dYdN*(muNet*cff1 + Q*mu0hatSI)
d2NPPdN2= 2.*PHY*mu0hatSI*Kn/(NO3+Kn)**3 * (cff1*(Kn-2.*NO3)/(NO3+Kn) + 1./Qmin)
! NPP: ensemble mean carbon based primary production
NPP = PHY*muNet*Q + .5*(2.*COVNP*dmuQ_dN + VNO3*d2NPPdN2)

! Calculate sources and sinks of variances of N, P, and covariance of NP
SVPHY = 2.*mu0hatSI*(fN*VPHY         + Kn*PHY*       COVNP/(Kn+NO3)**2)
SVNO3 =-2.*mu0hatSI*(fN*COVNP        + Kn*PHY*        VNO3/(Kn+NO3)**2) 
SCOVNP=    mu0hatSI*(fN*(COVNP-VPHY) + Kn*PHY*(VNO3-COVNP)/(Kn+NO3)**2) 
return
END SUBROUTINE PHY_NPCLOSURE

