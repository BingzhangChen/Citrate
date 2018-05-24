PROGRAM main
USE PARAM_MOD
IMPLICIT NONE

!constant temperature
real, parameter   :: temp_ = 15.

!constant light
real, parameter   :: PAR_  = 1000.

! Time step
real, parameter   :: dtsec   = 60.   ! Unit: seconds
real, parameter   :: d_per_s = 864d2 ! how many seconds in one day
integer, parameter:: nsave   = INT(d_per_s)/INT(dtsec) ! Save results once per day
!Total time (days) of model running
integer, parameter:: NDAYS = 360

!Total number of time steps
integer, parameter:: Nstep = NDAYS*INT(d_per_s)/INT(dtsec) 

! Total nutrient concentration:
real, parameter   :: A    = 2.

! Number of beta values
integer, parameter:: NB   = 5
real, parameter   :: beta(NB) = [0.01, 0.1,0.5, 1., 2.]

! Number of parameters to test
integer, parameter:: NPS  = 20

real    :: B 
real    :: current_sec 

! Biological tracers
real    :: NO3, PHY, ZOO, VNO3, VPHY, VZOO, COVNP, COVNZ, COVPZ

! Fluxes
real    :: PP_P_N   =0d0
real    :: PP_N_P   =0d0
real    :: PP_PP_NP =0d0
real    :: PP_NP_NN =0d0
real    :: PP_NP_PP =0d0
real    :: PP_PZ_PP =0d0
real    :: PP_Z_P   =0d0
real    :: PP_N_Z   =0d0
real    :: PP_NP_PZ =0d0
real    :: PP_NN_NP =0d0
real    :: PP_NN_NZ =0d0
real    :: PP_ZZ_PZ =0d0
real    :: PP_PZ_NZ =0d0
real    :: PP_NZ_PZ =0d0
real    :: PP_NZ_NP =0d0
real    :: PP_NZ_ZZ =0d0

! To count time
real    :: start, finish
integer :: i,j,k,p,it, current_day
character(20)       :: str
character(LEN=50)   :: outfile = ' '


call cpu_time(start) 

!Set parameters (the parameters will be back transformed to normal values in the main subroutine):
params(imu0)    = log(1.)
params(iIopt)   = log(1d3)
params(iaI0_C)  = log(0.05)
params(iDp)     = log(0.1)
params(iKPHY)   = log(1d0)

! Test the sensitivities of gmax and mz
MaxValue(igmax) = 5.
MinValue(igmax) = 0.01
MaxValue(imz)   = 0.8
MinValue(imz)   = 0.01
MaxValue(iKN)   = 3.
MinValue(iKN)   = 0.01

DO p = 1, NPS
 params(iKN) = log(A) + log((MaxValue(iKN) - MinValue(iKN))/dble(NPS-1)*dble(p-1)&
       + MinValue(iKN))

 DO j = 1, NPS
 
   params(igmax) = log((MaxValue(igmax) - MinValue(igmax))/dble(NPS-1)*dble(j-1)&
       + MinValue(igmax))
   DO i = 1, NPS
     params(imz) = log((MaxValue(imz) - MinValue(imz))/dble(NPS-1)*dble(i-1)&
       + MinValue(imz))
 
     DO k = 1, NB
         !Initialize:
         NO3   = 1.0
         PHY   = 0.6
         ZOO   = A - NO3 - PHY
 
         ! beta = B/A**2
         B     = A*A*beta(k)
         
         ! B   = VNO3 + VPHY + VZOO + 2*(COVNP + COVPZ + COVNZ)
         VNO3     =B*0.2
         VPHY     =B*0.3
         VZOO     =B*0.2
         COVNP    =B*.05
         COVPZ    =B*.05
         COVNZ    =B*.05
         PP_P_N   =0d0
         PP_N_P   =0d0
         PP_PP_NP =0d0
         PP_NP_NN =0d0
         PP_NP_PP =0d0
         PP_PZ_PP =0d0
         PP_Z_P   =0d0
         PP_N_Z   =0d0
         PP_NP_PZ =0d0
         PP_NN_NP =0d0
         PP_NN_NZ =0d0
         PP_ZZ_PZ =0d0
         PP_PZ_NZ =0d0
         PP_NZ_PZ =0d0
         PP_NZ_NP =0d0
         PP_NZ_ZZ =0d0
 
         ! Create data out file
         write(str, '(F12.2)') beta(k)
         str     = adjustl(str)
         outfile = 'beta'//trim(str)
         write(str, '(F12.2)') exp(params(igmax))
         str     = adjustl(str)
         outfile = trim(outfile)//'Gm'//trim(str)
         write(str, '(F12.2)') exp(params(imz))
         str     = adjustl(str)
         outfile = trim(outfile)//'mz'//trim(str)
         write(str, '(F12.2)') exp(params(iKN))
         str     = adjustl(str)
         outfile = trim(outfile)//'KN'//trim(str)

         open (unit=10, file = outfile, status = 'replace')
         write(10, 101) 'Timestep', 'Days', 'NO3', 'PHY', 'ZOO', 'VNO3', 'VPHY', &
           'VZOO', 'COVNP', 'COVNZ', 'COVPZ','P_N', 'N_P', 'Z_P','N_Z','PP_NP','NP_NN',& 
           'NP_PP','PZ_PP', 'NP_PZ',                                             &
           'NN_NP','NN_NZ','ZZ_PZ','PZ_NZ','NZ_PZ','NZ_NP','NZ_ZZ'
   
         write(10, 200) 0, 0, NO3,PHY,ZOO, VPHY,VNO3,VZOO,COVNP,COVNZ,COVPZ,PP_P_N,&
          PP_N_P, PP_Z_P,PP_N_Z, PP_PP_NP, PP_NP_NN,PP_NP_PP,PP_PZ_PP,PP_NP_PZ,PP_NN_NP,&
          PP_NN_NZ,PP_ZZ_PZ,PP_PZ_NZ,PP_NZ_PZ,PP_NZ_NP,PP_NZ_ZZ
         
         !the fraction of a time step in one day
         dtdays = dtsec/d_per_s
         
         do it = 1, (NSTEP+1)
         ! Calculate current timing (zero is starting time):
           current_sec = float(it-1)*dtsec
         
           call NPZ_CLOSURE
         
           IF (mod(it, nsave) .EQ. 1) THEN !Save results
               ! Calculate model time in days:
               current_day = int(current_sec/d_per_s)
               write(10, 200) it, current_day, NO3,PHY,ZOO, VPHY,VNO3,VZOO,COVNP,COVNZ,COVPZ,PP_P_N,&
                PP_N_P, PP_Z_P,PP_N_Z, PP_PP_NP, PP_NP_NN,PP_NP_PP,PP_PZ_PP,PP_NP_PZ,PP_NN_NP,&
                PP_NN_NZ,PP_ZZ_PZ,PP_PZ_NZ,PP_NZ_PZ,PP_NZ_NP,PP_NZ_ZZ
           ENDIF
         enddo
         close(10)
       ENDDO
   ENDDO
 ENDDO
ENDDO

Call cpu_time(finish)
print '("Time = ",f8.3," seconds.")', (finish-start)
101 format(100(A5,3x))
200 format(I10, 2x, I5, 100(2x, 1pe12.2))

CONTAINS 

SUBROUTINE NPZ_CLOSURE
! This NPZ_closure model has 5 tracers: <N>, <P>, <Z>,<N'>^2, <P'>^2,<Z'>^2, <N'P'>,<N'Z'>,<P'Z'>
! Governing functions modified following Priyadarshi et al. JTB (2017)
USE PARAM_MOD
IMPLICIT NONE
real :: QN  ! cell quota related variables
real :: muNet, Dp, muSI, KN, gmax
real :: SI, CFF1, CFF2
real :: theta,  Kp 
real :: Chl, NPPt, Zmort1, Zmort2

!Zooplankton feeding threshold for phytoplankton
real, parameter :: Pt_bio = 0.001d0

!-----------------------------------------------------------------------
KN   = exp(params(iKN))
Kp   = exp(params(iKPHY))

!Phytoplankton equations 
call PHY_NPCLOSURE(NO3,PAR_,temp_,PHY,VPHY,VNO3, COVNP,muSI,  &
        muNet,SI,theta,QN, PP_P_N, PP_N_P, Chl, NPPt, PP_PP_NP, PP_NP_NN)
!=============================================================
!! Solve ODE functions:
!All rates have been multiplied by dtdays to get the real rate correponding to the actual time step
Dp     = exp(params(iDp))*exp(params(imu0))  ! Mortality rate of phytoplankton
tf_P   = tempBOL(Ep, temp_)
tf_Z   = tempBOL(Ez, temp_)
PP_N_P = PP_N_P + PHY*Dp*tf_P

! params(igmax) is scaled factor to imu0
! Because temp = 15 ºC, no effect of temperature actually
gmax  = exp(params(igmax)) * exp(params(imu0)) * tf_Z
Zmort = exp(params(imz))   * exp(params(imu0)) * tf_Z   ! Mortality rate of zooplankton
Select case(ZOO_MORT)
case(LINEAR)
   Zmort1 = Zmort * ZOO
   Zmort2 = Zmort
case(QUADRATIC)
   Zmort1 = Zmort *(ZOO*ZOO + VZOO)
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

PP_N_Z = (1. - GGE)*INGES + Zmort1
PP_Z_P = INGES

!Update tracers:
NO3 = NO3  + (PP_N_Z - PP_P_N + PP_N_P)*dtdays
PHY = PHY  + (PP_P_N - PP_Z_P - PP_N_P)*dtdays
ZOO = ZOO  + (PP_Z_P - PP_N_Z)*dtdays

SI   = NO3/(NO3 + KN)

PP_NP_PP = Dp * tf_P * VPHY
PP_PZ_PP = gmax*(CFF2*ZOO*VPHY + CFF1*COVPZ)

PP_PZ_NZ = muSI*(SI * COVPZ + PHY * COVNZ * KN / (KN + NO3)**2 )
PP_ZZ_PZ = gmax * (CFF2 * ZOO * COVPZ + CFF1 * VZOO)
PP_NZ_PZ = Dp*tf_P*COVPZ
PP_NZ_NP = gmax * (CFF2 * COVNP * ZOO + CFF1 * COVNZ)
PP_NZ_ZZ = Zmort2 * VZOO + gmax * (1.-GGE)*(CFF2*ZOO*COVPZ + CFF1 *VZOO)
PP_NP_PZ = Zmort2*COVPZ + gmax*(1.-GGE)*(COVPZ*CFF1 + VPHY*ZOO*CFF2) 
PP_NN_NP = Dp*COVNP*tf_P
PP_NN_NZ = (CFF2*ZOO*COVNP + CFF1*COVNZ)*(1.-GGE)*gmax + Zmort2*COVNZ
VPHY     = VPHY + 2.*(PP_PP_NP - PP_NP_PP - PP_PZ_PP)*dtdays
VNO3     = VNO3 + 2.*(PP_NN_NP + PP_NN_NZ - PP_NP_NN)*dtdays
VZOO     = VZOO + 2.*(PP_ZZ_PZ - PP_NZ_ZZ           )*dtdays
COVNP    = COVNP+    (PP_NP_PZ + PP_NP_PP + PP_NP_NN - PP_PP_NP - PP_NN_NP - PP_NZ_NP)*dtdays
COVNZ    = COVNZ+    (PP_NZ_NP + PP_NZ_ZZ + PP_NZ_PZ - PP_NN_NZ - PP_PZ_NZ)*dtdays
COVPZ    = COVPZ+    (PP_PZ_NZ + PP_PZ_PP - PP_NZ_PZ - PP_NP_PZ - PP_ZZ_PZ)*dtdays

!VPHY = VPHY + (SVPHY - 2.*Dp*VPHY*tf_P                    &
!   - 2.*gmax*(CFF2 * ZOO * VPHY + CFF1 * COVPZ) )*dtdays

!VNO3 = VNO3 + (SVNO3 + 2.*Dp*COVNP*tf_P                   &
!     + 2.*(1.-GGE)*gmax*(CFF2 * ZOO * COVNP + CFF1 * COVNZ)          &
!     + 2.* Zmort2 * COVNZ) * dtdays

!VZOO = VZOO + 2.*dtdays*(GGE*gmax*( CFF2 * ZOO * COVPZ    &
!     + CFF1 * VZOO) - 2. * Zmort2 * VZOO )

!COVNP = COVNP + (SCOVNP + Dp*tf_P*(VPHY - COVNP)                     &
!    +     Zmort2    * COVPZ                                           &
!    +     gmax*ZOO * CFF2 * ((1.-GGE)* VPHY - COVNP)                 &
!    +     gmax     * CFF1 * ((1.-GGE)*COVPZ - COVNZ) )*dtdays

!COVPZ = COVPZ + ( NPPt                                               &
!    - (Dp*tf_P +  Zmort2)* COVPZ                                     &
!    + gmax * ZOO * CFF2 * (GGE * VPHY - COVPZ)                       &
!    + gmax       * CFF1 * (GGE * COVPZ - VZOO))*dtdays

!COVNZ = COVNZ + (-NPPt + Dp*tf_P*COVPZ                               & 
!    + Zmort2*(VZOO - COVNZ)                                          &
!    + gmax * ZOO  * CFF2 * (GGE*COVNP + (1.-GGE)*COVPZ)              &
!    + gmax        * CFF1 * (GGE*COVNZ + (1.-GGE)*VZOO))*dtdays

return
END SUBROUTINE NPZ_CLOSURE 

END PROGRAM
