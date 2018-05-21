SUBROUTINE NPZ_CLOSURE
    ! This NPZ_closure model has 5 tracers: <N>, <P>, <Z>,<N'>^2, <P'>^2,<Z'>^2, <N'P'>,<N'Z'>,<P'Z'>
    ! Governing functions modified following Priyadarshi et al. JTB (2017)
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
real :: SI, CFF1, CFF2
real :: theta, Snp, Kp
real :: Chl, NPPt, Zmort1, Zmort2

!Zooplankton feeding threshold for phytoplankton
real, parameter :: Pt_bio = 0.1d0

!-----------------------------------------------------------------------
KN   = exp(params(iKN))
Kp   = exp(params(iKPHY))
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
     Varout(oNPP, k) = NPPt * 12d0 ! Correct the unit to ug C/L
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

   ! params(igmax) is scaled factor to imu0
   gmax  = exp(params(igmax)) * exp(params(imu0)) * tf_Z
   Zmort = exp(params(imz))   * gmax   ! Mortality rate of zooplankton

   Select case(ZOO_MORT)
   case(LINEAR)
      Zmort1 = Zmort * ZOO
      Zmort2 = Zmort
   case(QUADRATIC)
      Zmort1 = Zmort * (ZOO**2 + VZOO)
      Zmort2 = 2.*Zmort * ZOO
   case default
      stop "The correct zooplankton mortality option is NOT selected!"
   end select
   ! Define two scratch variables that are repeatly used below
   select case(grazing_formulation)
   case(LINEAR)
      CFF1  = PHY
      CFF2  = 1d0
      INGES = gmax*(CFF1 * ZOO + COVPZ)
   case(H2T) !Holling type 2 with threshold
      NPPt   = PHY - Pt_bio
      if (NPPt < 0.) then
         CFF1  = 0d0
         CFF2  = 1d0
         INGES = 0d0
      else
         CFF1  = NPPt/ (Kp + NPPt)
         CFF2  = Kp  / (Kp + NPPt)**2
         INGES = gmax*( CFF1 * ZOO                                      &
               - Kp * ZOO * VPHY /(Kp + NPPt)**3 + CFF2 * COVPZ)
      endif

   case(H3)  !Holling type 3
      CFF1  =    PHY**2 / (Kp**2 + PHY**2)
      CFF2  = 2.*Kp**2 * PHY / (Kp**2 + PHY**2)**2

      INGES = gmax*( CFF1 * ZOO                                          &
            + Kp**2 * ZOO * VPHY*(Kp**2 - 3.*PHY**2)/(Kp**2 + PHY**2)**3 &
            + CFF2 * COVPZ )
   case default
      stop "The correct option for zooplankton grazing function is NOT selected!"
   end select

   PP_NZ = (1. - GGE)*INGES + Zmort1
   PP_ZP =        GGE*INGES

!Update tracers:
   Varout(iNO3, k) = max(NO3  + (PP_NZ - PP_PN )*dtdays, eps)
   Varout(iPHY, k) = max(PHY  + (PP_PN - INGES )*dtdays, eps)
   Varout(iZOO, k) = max(ZOO  + (PP_ZP - Zmort1)*dtdays, eps)

   SI   = NO3/(NO3 + KN)
   NPPt = muSI*(SI * COVPZ + PHY * COVNZ * KN / (KN + NO3)**2 )

   Varout(iVPHY,k) = VPHY + (SVPHY - 2.*Dp*VPHY*tf_P                    &
      - 2.*gmax*(CFF2 * ZOO * VPHY + CFF1 * COVPZ) )*dtdays

   Varout(iVNO3,k) = VNO3 + (SVNO3 + 2.*Dp*COVNP*tf_P                   &
        + 2.*(1.-GGE)*gmax*(CFF2 * ZOO * COVNP + CFF1 * COVNZ)          &
        + 2.* Zmort2 * COVNZ) * dtdays

   Varout(iVZOO,k) = VZOO + 2.*dtdays*(GGE*gmax*( CFF2 * ZOO * COVPZ    &
        + CFF1 * VZOO) - 2. * Zmort2 * VZOO )

   Varout(iCOVNP,k)= COVNP + (SCOVNP + Dp*tf_P*(VPHY - COVNP)           &
       +     Zmort    * COVPZ                                           &
       +     gmax*ZOO * CFF2 * ((1.-GGE)* VPHY - COVNP)                 &
       +     gmax     * CFF1 * ((1.-GGE)*COVPZ - COVNZ) )*dtdays

   Varout(iCOVPZ,k) = COVPZ + ( NPPt                                    &
       - (Dp*tf_P +  Zmort2)* COVPZ                                     &
       + gmax * ZOO * CFF2 * (GGE * VPHY - COVPZ)                       &
       + gmax       * CFF1 * (GGE * COVPZ - VZOO))*dtdays

   Varout(iCOVNZ,k) = COVNZ + (-NPPt + Dp*tf_P*COVPZ                    & 
       + Zmort2*(VZOO - COVNZ)                                          &
       + gmax * ZOO  * CFF2 * (GGE*COVNP + (1.-GGE)*COVPZ)              &
       + gmax        * CFF1 * (GGE*COVNZ + (1.-GGE)*VZOO))*dtdays

   Varout(omuNet, k) = muNet               !Growth rate at <N>
   Varout(iVPHY,k)   = max(Varout(iVPHY,k), eps)
   Varout(iVNO3,k)   = max(Varout(iVNO3,k), eps)
   Varout(iVZOO,k)   = max(Varout(iVZOO,k), eps)
   Varout(oPON, k)   = PHY + ZOO
ENDDO
END SUBROUTINE NPZ_CLOSURE 
