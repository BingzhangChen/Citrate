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
integer, parameter:: NDAYS= 360.

!Total number of time steps
integer, parameter:: Nstep = NDAYS*INT(d_per_s)/INT(dtsec) 

! Total nutrient concentration:
real, parameter   :: A    = 2.

real    :: B, beta
real    :: current_sec 

real    :: NO3, PHY, ZOO, VNO3, VPHY, VZOO, COVNP, COVNZ, COVPZ
integer :: i,j,k,it, current_day
character(20)       :: str
character(LEN=20)   :: outfile = ' '

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
NO3   = 1.5
PHY   = 0.4
ZOO   = A - NO3 - PHY
beta  = 0.1
B     = A*A*beta

! B   = VNO3 + VPHY + VZOO + 2*(COVNP + COVPZ + COVNZ)
VNO3  = B*0.2
VPHY  = B*0.2
VZOO  = B*0.2
COVNP = B*.1
COVPZ = B*.05
COVNZ = B*.05

! Create data out file
write(str, '(F12.2)') beta
str     = adjustl(str)
outfile = 'beta'//str//'.out'

open (unit=10, file = outfile, status = 'replace')
write(10, 101) 'Timestep', 'Days', 'NO3', 'PHY', 'ZOO', 'VNO3', 'VPHY',   &
  'VZOO', 'COVNP', 'COVNZ', 'COVPZ' 
101 format(100(A5,3x))

!the fraction of a time step in one day
dtdays = dtsec/d_per_s

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
200 format(I10, 2x, I5, 100(2x, 1pe12.3))
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
real :: theta, Snp, Kp
real :: Chl, NPPt, Zmort1, Zmort2

!Zooplankton feeding threshold for phytoplankton
real, parameter :: Pt_bio = 0.1d0

!-----------------------------------------------------------------------
KN   = exp(params(iKN))
Kp   = exp(params(iKPHY))

if (NO3 < 0.) then
  write(6,*) "NO3 is negative!"
  return
endif
if (PHY < 0.) then
  write(6,*) "PHY is negative!"
  return
endif
if (ZOO < 0.) then
  write(6,*) "ZOO is negative!"
  return
endif
if (VPHY < 0.) then
  write(6,*) "VPHY is negative!"
  return
endif
if (VZOO < 0.) then
  write(6,*) "VZOO is negative!"
  return
endif
if (VNO3 < 0.) then
  write(6,*) "VNO3 is negative!"
  return
endif

!Phytoplankton equations 
call PHY_NPCLOSURE(NO3,par_,temp_,PHY,VPHY,VNO3, COVNP, muSI, muNet,SI,&
  theta,QN, Snp, Chl, NPPt, SVPHY, SVNO3, SCOVNP)
!=============================================================
!! Solve ODE functions:
!All rates have been multiplied by dtdays to get the real rate correponding to the actual time step
Dp    = exp(params(iDp))*exp(params(imu0))  ! Mortality rate of phytoplankton
tf_P  = tempBOL(Ep, temp_)
PP_PN = Snp - PHY*Dp*tf_P
tf_Z  = tempBOL(Ez, temp_)

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
NO3 = NO3  + (PP_NZ - PP_PN )*dtdays
PHY = PHY  + (PP_PN - INGES )*dtdays
ZOO = ZOO  + (PP_ZP - Zmort1)*dtdays

SI   = NO3/(NO3 + KN)
NPPt = muSI*(SI * COVPZ + PHY * COVNZ * KN / (KN + NO3)**2 )

VPHY = VPHY + (SVPHY - 2.*Dp*VPHY*tf_P                    &
   - 2.*gmax*(CFF2 * ZOO * VPHY + CFF1 * COVPZ) )*dtdays

VNO3 = VNO3 + (SVNO3 + 2.*Dp*COVNP*tf_P                   &
     + 2.*(1.-GGE)*gmax*(CFF2 * ZOO * COVNP + CFF1 * COVNZ)          &
     + 2.* Zmort2 * COVNZ) * dtdays

VZOO = VZOO + 2.*dtdays*(GGE*gmax*( CFF2 * ZOO * COVPZ    &
     + CFF1 * VZOO) - 2. * Zmort2 * VZOO )

COVNP = COVNP + (SCOVNP + Dp*tf_P*(VPHY - COVNP)                     &
    +     Zmort    * COVPZ                                           &
    +     gmax*ZOO * CFF2 * ((1.-GGE)* VPHY - COVNP)                 &
    +     gmax     * CFF1 * ((1.-GGE)*COVPZ - COVNZ) )*dtdays

COVPZ = COVPZ + ( NPPt                                               &
    - (Dp*tf_P +  Zmort2)* COVPZ                                     &
    + gmax * ZOO * CFF2 * (GGE * VPHY - COVPZ)                       &
    + gmax       * CFF1 * (GGE * COVPZ - VZOO))*dtdays

COVNZ = COVNZ + (-NPPt + Dp*tf_P*COVPZ                               & 
    + Zmort2*(VZOO - COVNZ)                                          &
    + gmax * ZOO  * CFF2 * (GGE*COVNP + (1.-GGE)*COVPZ)              &
    + gmax        * CFF1 * (GGE*COVNZ + (1.-GGE)*VZOO))*dtdays
return
END SUBROUTINE NPZ_CLOSURE 
