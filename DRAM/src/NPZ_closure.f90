SUBROUTINE NPZ_CLOSURE
    ! This NPZ_closure model has 5 tracers: <N>, <P>, <Z>,<N'>^2, <P'>^2,<Z'>^2, <N'P'>,<N'Z'>,<P'Z'>
    ! Governing functions follow Priyadarshi et al. JTB (2017)
USE PARAM_MOD
USE MOD_1D, only: it, nsave
IMPLICIT NONE
integer :: k
!INPUT PARAMETERS:
real :: par_
!LOCAL VARIABLES of phytoplankton:
real :: NO3,PHY, ZOO,VPHY, VNO3,VZOO, COVNP,COVPZ,COVNZ 
real :: SVNO3, SVPHY, SCOVNP
real :: QN  ! cell quota related variables
real :: muNet, Dp, PP_PN, muSI, KN, gmax
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
   NO3    = Vars(iNO3,    k)
   NO3    = max(NO3,    eps)
   PHY    = max(Vars(iPHY(1),k), eps)
   ZOO    = Vars(iZOO,    k)
   ZOO    = max(ZOO,    eps)
   VPHY   = Vars(iVPHY,   k)
   VPHY   = max(VPHY,   eps)
   VNO3   = Vars(iVNO3,   k)
   VNO3   = max(VNO3,   eps)
   VZOO   = max(Vars(iVZOO,k),eps)
   COVNP  = Vars(iCOVNP,  k)
   COVPZ  = Vars(iCOVPZ,  k)
   COVNZ  = Vars(iCOVNZ,  k)
   IF (mod(it, nsave) .EQ. 1) THEN
     !Use nonmixed light to calculate NPP
     call PHY_NPCLOSURE(NO3,PAR_,Temp(k),PHY,VPHY,VNO3, COVNP,muSI, muNet,SI,theta,QN, &
                        Snp, Chl, NPPt, SVPHY, SVNO3, SCOVNP)
     Varout(oNPP, k)  = NPPt * 12d0 ! Correct the unit to ug C/L
   ENDIF
   !Use mixed light to estimate source and sink
   call PHY_NPCLOSURE(NO3,PAR_,Temp(k),PHY,VPHY,VNO3, COVNP, muSI, muNet,SI,theta,QN, &
                        Snp, Chl, NPPt, SVPHY, SVNO3, SCOVNP)
   ! Save some model outputs:
   Varout(oTheta,k) = theta! Chl:C ratio at <N>
   Varout(oQN,   k) = QN   ! N:C ratio at <N> 
   Varout(oSI,   k) = SI   ! Light limitation
   Varout(oCHLt, k) = Chl  ! ensemble mean Chl
!=============================================================
!! Solve ODE functions:
!All rates have been multiplied by dtdays to get the real rate correponding to the actual time step
   Dp    = exp(params(iDp))*exp(params(imu0))  ! Mortality rate of phytoplankton
   tf_P  = TEMPBOL(Ep, Temp(k))
   PP_PN = Snp - PHY*Dp*tf_P
   tf_Z  = TEMPBOL(Ez, Temp(k))
   gmax  = exp(params(iGmax))*tf_Z
   Zmort = exp(params(imz))  *tf_Z   ! Mortality rate of zooplankton
   INGES = gmax*(PHY*ZOO + COVPZ)
   PP_NZ = (1.-GGE)*INGES + Zmort*ZOO
   PP_ZP = GGE*INGES

!Update tracers:
   Varout(iNO3, k) = max(NO3  + (PP_NZ - PP_PN)*dtdays, eps)
   Varout(iPHY, k) = max(PHY  + (PP_PN - INGES)*dtdays, eps)
   Varout(iZOO, k) = ZOO  + (PP_ZP - Zmort*ZOO)*dtdays
   Varout(iZOO, k) = max(Varout(iZOO,k), eps)
   Varout(iVPHY,k) = VPHY + (SVPHY - 2.*Dp*VPHY*tf_P - 2.*gmax*(ZOO*VPHY + PHY*COVPZ))*dtdays
   Varout(iVPHY,k) = max(Varout(iVPHY,k), eps)
   Varout(iVNO3,k) = VNO3 + (SVNO3 + 2.*Dp*COVNP*tf_P        &
        + 2.*(1.-GGE)*gmax*(ZOO*COVNP + PHY*COVNZ)           &
        + 2.*Zmort*COVNZ)*dtdays
   Varout(iVNO3,k) = max(Varout(iVNO3,k), eps)

   Varout(iVZOO,k) = VZOO + 2.*dtdays*(GGE*gmax*(ZOO*COVPZ+PHY*VZOO) - Zmort*VZOO)
   Varout(iVZOO,k) = max(Varout(iVZOO,k), eps)

   Varout(iCOVNP,k)= COVNP + (SCOVNP+ Dp*tf_P*(VPHY-COVNP) + Zmort*COVPZ + &
   gmax*ZOO*((1.-GGE)*VPHY - COVNP) +gmax*PHY*((1.-GGE)* COVPZ - COVNZ) )*dtdays

   KN   = exp(params(iKN))
   SI   = NO3/(NO3+KN)
   NPPt = muSI*(SI*COVPZ + PHY*COVNZ*KN/(KN+NO3)**2)

   Varout(iCOVPZ,k) = COVPZ + ( NPPt &
   - (Dp*tf_P+Zmort)*COVPZ + gmax*ZOO*(GGE*VPHY - COVPZ) +&
   gmax * PHY*(GGE*COVPZ - VZOO))*dtdays

   Varout(iCOVNZ,k) = COVNZ + (-NPPt + Dp*tf_P*COVPZ + Zmort*(VZOO-COVNZ) &
   + gmax*ZOO*(GGE*COVNP + (1.-GGE)*COVPZ) &
   + gmax*PHY*(GGE*COVNZ + (1.-GGE)*VZOO))*dtdays

   Varout(omuNet, k) = muNet               !Growth rate at <N>
ENDDO
END SUBROUTINE NPZ_CLOSURE 
