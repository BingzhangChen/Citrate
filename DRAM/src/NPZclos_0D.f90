program main
USE PARAM_MOD
IMPLICIT NONE
real, parameter   :: tC   = 15.
real, parameter   :: PAR_ = 1000.
! Time step
real, parameter   :: dtsec   = 60.
real, parameter   :: d_per_s = 864d2 ! how many seconds in one day
integer, parameter:: nsave   = INT(d_per_s)/INT(dtsec) ! Timesteps to save
!Total time (days) of model running
integer, parameter:: NDAYS = 360.

!Total number of time steps
integer, parameter:: Nstep = NDAYS*INT(d_per_s)/INT(dtsec) 

! Total nutrient concentration:
real, parameter   :: A    = 2.

real    :: B, beta
real    :: current_sec 

real    :: NO3, PHY, ZOO, VNO3, VPHY, VZOO, COVNP, COVNZ, COVPZ

! To count time
real    :: start, finish
integer :: i,j,k,it, current_day
character(20)       :: str
character(LEN=20)   :: outfile = ' '


call cpu_time(start) 
!Set parameters:
params(imu0)    = log(2.)
params(igmax)   = log(1.)
params(iKN)     = log(A * 0.5)
params(iIopt)   = log(1d3)
params(iaI0_C)  = log(0.05)
params(iDp)     = log(0.4)
params(imz)     = log(0.07)
params(iKPHY)   = log(1d0)

!Initialize:
NO3   = 1.0
PHY   = 0.5
ZOO   = A - NO3 - PHY
beta  = 0.1
B     = A*A*beta

! B   = VNO3 + VPHY + VZOO + 2*(COVNP + COVPZ + COVNZ)
VNO3  = B*0.3
VPHY  = B*0.3
VZOO  = B*0.4
COVNP = B*.0
COVPZ = B*.0
COVNZ = B*.0

! Create data out file
write(str, '(F12.2)') beta
str     = adjustl(str)
outfile = 'beta'//str//'.out'

open (unit=10, file = outfile, status = 'replace')
write(10, 101) 'Timestep', 'Days', 'NO3', 'PHY', 'ZOO', 'VNO3', 'VPHY',   &
  'VZOO', 'COVNP', 'COVNZ', 'COVPZ' 
write(10, 200) 0, 0, NO3,PHY,ZOO, VPHY,VNO3,VZOO,COVNP,COVNZ,COVPZ

101 format(100(A5,3x))

!the fraction of a time step in one day
dtdays = dtsec/d_per_s
write(6,*) 'dtdays = ',dtdays

do it = 1, (NSTEP+1)
! Calculate current timing (zero is starting time):
  current_sec = float(it-1)*dtsec

  call NPZ_CLOSURE(PAR_, tC, NO3, PHY, ZOO, VPHY, VNO3,VZOO,COVNP,&
    COVPZ,COVNZ)
  IF (mod(it, nsave) .EQ. 1) THEN !Save results
      ! Calculate model time in days:
      current_day = int(current_sec/d_per_s)
      write(10, 200) it, current_day, NO3,PHY,ZOO, VPHY,VNO3,VZOO,COVNP,&
        COVNZ,COVPZ
  ENDIF
enddo
close(10)
call cpu_time(finish)
print '("Time = ",f8.3," seconds.")', (finish-start)
200 format(I10, 2x, I5, 100(2x, 1pe12.2))
end program

SUBROUTINE NPZ_CLOSURE(par_, temp_, NO3, PHY, ZOO, VPHY, VNO3,VZOO,COVNP,&
    COVPZ,COVNZ)
! This NPZ_closure model has 5 tracers: <N>, <P>, <Z>,<N'>^2, <P'>^2,<Z'>^2, <N'P'>,<N'Z'>,<P'Z'>
! Governing functions modified following Priyadarshi et al. JTB (2017)
USE PARAM_MOD
IMPLICIT NONE
!INPUT/output PARAMETERS:
real, intent(in)    :: par_, temp_
real, intent(inout) :: NO3,PHY, ZOO, VPHY, VNO3,VZOO,COVNP, COVPZ,COVNZ

real :: SVNO3, SVPHY, SCOVNP
real :: QN  ! cell quota related variables
real :: muNet, Dp, PP_PN, muSI, KN, gmax
real :: SI, CFF1, CFF2
real :: theta,  Kp, PP_P_N, PP_N_P, PP_PP_NP, PP_NP_NN, PP_NP_PP, PP_PZ_PP
real :: PP_NP_PZ
real :: PP_NN_NP
real :: PP_NN_NZ
real :: PP_ZZ_PZ
real :: PP_PZ_NZ
real :: PP_NZ_PZ
real :: PP_NZ_NP
real :: PP_NZ_ZZ
real :: Chl, NPPt, Zmort1, Zmort2

!Zooplankton feeding threshold for phytoplankton
real, parameter :: Pt_bio = 0.1d0

!-----------------------------------------------------------------------
KN   = exp(params(iKN))
Kp   = exp(params(iKPHY))

if (NO3 < 0.) then
  write(6,*) "NO3 is negative!"
  stop
endif
if (PHY < 0.) then
  write(6,*) "PHY is negative!"
  stop
endif
if (ZOO < 0.) then
  write(6,*) "ZOO is negative!"
  stop
endif
if (VPHY < 0.) then
  write(6,*) "VPHY is negative!"
  stop
endif
if (VZOO < 0.) then
  write(6,*) "VZOO is negative!"
  stop
endif
if (VNO3 < 0.) then
  write(6,*) "VNO3 is negative!"
  stop
endif

!Phytoplankton equations 
call PHY_NPCLOSURE(NO3,PAR_,Temp_,PHY,VPHY,VNO3, COVNP,muSI,  &
        muNet,SI,theta,QN, PP_P_N, PP_N_P, Chl, NPP, PP_PP_NP, PP_NP_NN)
!=============================================================
!! Solve ODE functions:
!All rates have been multiplied by dtdays to get the real rate correponding to the actual time step
Dp     = exp(params(iDp))*exp(params(imu0))  ! Mortality rate of phytoplankton
tf_P   = tempBOL(Ep, temp_)
tf_Z   = tempBOL(Ez, temp_)
PP_N_P = PP_N_P + PHY*Dp*tf_P

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

PP_N_Z = (1. - GGE)*INGES + Zmort1
PP_Z_P = INGES

!Update tracers:
NO3 = NO3  + (PP_N_Z - PP_P_N + PP_N_P)*dtdays
PHY = PHY  + (PP_P_N - PP_Z_P - PP_N_P)*dtdays
ZOO = ZOO  + (PP_Z_P - PP_N_Z)*dtdays

SI   = NO3/(NO3 + KN)
NPPt = muSI*(SI * COVPZ + PHY * COVNZ * KN / (KN + NO3)**2 )

PP_NP_PP = Dp * tf_P * VPHY
PP_PZ_PP = gmax*(CFF2*ZOO*VPHY + CFF1*COVPZ)

!VPHY = VPHY + (SVPHY - 2.*Dp*VPHY*tf_P                    &
!   - 2.*gmax*(CFF2 * ZOO * VPHY + CFF1 * COVPZ) )*dtdays

PP_ZZ_PZ = gmax * (CFF2 * ZOO * COVPZ + CFF1 * VZOO)
PP_PZ_NZ = NPPt
PP_NZ_PZ = Dp*tf_P*COVPZ
PP_NZ_NP = gmax * (CFF2 * COVNP * ZOO + CFF1 * COVNZ)
PP_NZ_ZZ = Zmort2 * VZOO + gmax * (1.-GGE)*(CFF2*ZOO*COVPZ    &
         + CFF1 * VZOO)
PP_NP_PZ = Zmort2*COVPZ + gmax*(1.-GGE)*(COVPZ*CFF1 + VPHY*ZOO*CFF2) 
PP_NN_NP = Dp*COVNP*tf_P
PP_NN_NZ = (CFF2*ZOO*COVNP + CFF1*COVNZ)*(1.-GGE)*gmax + Zmort2*COVNZ
VPHY     = VPHY + 2.*(PP_PP_NP - PP_NP_PP - PP_PZ_PP)*dtdays
VNO3     = VNO3 + 2.*(PP_NN_NP + PP_NN_NZ - PP_NP_NN)*dtdays
VZOO     = VZOO + 2.*(PP_ZZ_PZ - PP_NZ_ZZ           )*dtdays
COVNP    = COVNP+    (PP_NP_PZ + PP_NP_PP + PP_NP_NN - PP_PP_NP - PP_NN_NP - PP_NZ_NP)*dtdays
COVNZ    = COVNZ+    (PP_NZ_NP + PP_NZ_ZZ + PP_NZ_PZ - PP_NN_NZ - PP_PZ_NZ)*dtdays
COVPZ    = COVPZ+    (PP_PZ_NZ + PP_PZ_PP - PP_NZ_PZ - PP_NP_PZ - PP_ZZ_PZ)*dtdays
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
