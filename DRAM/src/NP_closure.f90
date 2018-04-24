SUBROUTINE NP_CLOSURE
    ! This NP_closure model has 5 tracers: <N>, <P>, <N'>^2, <P'>^2, <N'P'>
    ! Governing functions follow Mandal et al. JPR (2016)
USE bio_MOD
USE MOD_1D, only: it, nsave
IMPLICIT NONE
integer :: k
!INPUT PARAMETERS:
real :: par_
!LOCAL VARIABLES of phytoplankton:
real :: NO3,PHY, VPHY, VNO3, COVNP, SVNO3, SVPHY, SCOVNP
real :: QN  ! cell quota related variables
real :: muNet, Dp, PP_PN, muSI 
real :: SI
real :: theta, Snp
real :: Chl, NPPt
!-----------------------------------------------------------------------
DO k = nlev, 1, -1   
   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif
   Varout(oPAR_,k) = par_
   tf_P   = TEMPBOL(Ep, Temp(k))
   NO3    = Vars(iNO3,    k)
   NO3    = max(NO3,    eps)
   PHY    = Vars(iPHY(1), k)
   VPHY   = Vars(iVPHY,   k)
   VPHY   = max(VPHY,    0.)
   VNO3   = Vars(iVNO3,   k)
   VNO3   = max(VNO3,    0.)
   COVNP  = Vars(iCOVNP,  k)
   IF (mod(it, nsave) .EQ. 1) THEN
     !Use nonmixed light to calculate NPP
     call PHY_NPCLOSURE(NO3,PAR_,Temp(k),PHY,VPHY,VNO3, COVNP,muSI, muNet,SI,theta,QN, &
                        Snp, Chl, NPPt, SVPHY, SVNO3, SCOVNP)
     Varout(oPPt, k)  = NPPt * 12d0 ! Correct the unit to ug C/L
   ENDIF
   !Use mixed light to estimate source and sink
   call PHY_NPCLOSURE(NO3,PAR_,Temp(k),PHY,VPHY,VNO3, COVNP,muSI, muNet,SI,theta,QN, &
                        Snp, Chl, NPPt, SVPHY, SVNO3, SCOVNP)
   ! Save some model outputs:
   Varout(oTheta(1),k)= theta! Chl:C ratio at <N>
   Varout(oQN(1)   ,k)= QN   ! N:C ratio at <N> 
   Varout(oSI(1),   k)= SI   ! Light limitation
   Varout(oCHLt,    k)= Chl  ! ensemble mean Chl
   Varout(oCHL(1),  k)= Chl  ! ensemble mean Chl
!=============================================================
!! Solve ODE functions:
!All rates have been multiplied by dtdays to get the real rate correponding to the actual time step
   Dp   = exp(params(iDp))  ! Mortality rate of phytoplankton
   PP_PN= Snp - PHY*Dp*tf_P
   PP_PN= min(PP_PN, (NO3-eps)/dtdays)  !To make NO3 positive
   PP_PN= max(PP_PN, (eps-PHY)/dtdays)  !To make PHY positive

!Update tracers:
   NO3  = NO3   -               PP_PN*dtdays
   PHY  = PHY   +               PP_PN*dtdays
   VPHY = VPHY  + (SVPHY-2.*Dp*VPHY*tf_P )*dtdays
   VNO3 = VNO3  + (SVNO3+2.*Dp*COVNP*tf_P)*dtdays
   COVNP= COVNP + (SCOVNP+Dp*(VPHY-COVNP)*tf_P)*dtdays
   Varout(oNO3,k)      = NO3
   Varout(oPHY(1),k)   = PHY
   Varout(oVPHY  ,k)   = max(VPHY,0d0)
   Varout(oVNO3  ,k)   = max(VNO3,0d0)
   Varout(oCOVNP ,k)   = COVNP
   Varout(oPHYt,  k)   = PHY
   Varout(omuNet(1),k) = muNet               !Growth rate at <N>
   Varout(omuAvg,   k) = Snp/PHY             !Ensemble mean growth rate 
ENDDO
END SUBROUTINE NP_CLOSURE 
