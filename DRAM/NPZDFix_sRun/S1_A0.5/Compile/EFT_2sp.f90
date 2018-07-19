SUBROUTINE FLEXEFT_2SP
! Two phytoplankton species model, one with higher muMax, lower KN and aI0,
! the other with lower  muMax, higher KN and aI0
use bio_MOD
implicit none
integer :: i,k
real :: dx,par_
real :: NO3, ZOO, DET
real :: PHYtot = 0d0
real :: PHYtot2= 0d0
real :: ppC_PN = 0d0
real :: pp_PN  = 0d0
real :: CHLtot = 0d0
real :: alphaG,Kp,gmax,RDN,mz,PHYtot1,muNet

! The parameter that describes how feeding preference decreases with cell size
real,    parameter :: alphaGmax = -0.05
real :: RL(NPHY)  ! Feeding preference on each size class
real :: mu0(2) = [5.0, 1.0]
real ::  KN(2) = [0.1, 1.0]
real :: A0N(2) = [1D4, 0.001]
real :: aI0(2) = [0.05, 0.4]

mu0(1)=params(imu0)
mu0(2)=mu0(1)*params(imu0B)
A0N(1)=params(iA0N)
A0N(2)=A0N(1)*10**(params(iA0N2))
aI0(1)=params(iaI0)
aI0(2)=aI0(1)*10**(params(iaI0B))

Ez   = 0.6
Ep   = 0.5
Ez   = 0.6
Kp   = 0.5
gmax = params(igmax)
RDN  = 0.1
mz   = params(imz)

!The feeding preference of the generic zooplankton for prey i
RL(1)=1d0
RL(2)=RL(1)*params(iRL2)

alphaG=1d0+10**params(ialphaG)

DO k = nlev,1,-1
   PHYtot  = 0d0  ! Calculate total PHY biomass
   PHYtot1 = 0d0  ! Total palatable prey
   PHYtot2 = 0d0  ! total P**alphaG
   pp_PN   = 0d0  ! total primary production (mmol N per d per m3)
   ppC_PN  = 0d0  ! total primary production (mmol C per d per m3)
   CHLtot  = 0d0  ! Total CHL A
   DET     = Vars(iDET,k)
   ZOO     = Vars(iZOO,k)  !Zooplankton biomass
   NO3     = Vars(iNO3,k)

   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif
   Varout(oPAR_,k)=par_

!! Phytoplankton section:
   do i=1,NPHY
    ! Calculate phytoplankton growth rate, theta, and QN based on environmental conditions
     call EFT_phygrowth(mu0(i),KN(i),A0N(i), aI0(i), NO3,Temp(k),PAR(k),&
                     muNet, Varout(oQN(i),k),                           &
                     Varout(oTheta(i),k),                               &
                     Varout(oSI(i),k), Varout(oLno3(i),k))


    ! Calculate phytoplankton growth rate, theta, and QN based on environmental conditions
     call EFT_phygrowth(mu0(i),KN(i),A0N(i), aI0(i), NO3,Temp(k), par_, &
                     Varout(omuNet(i),k), Varout(oQN(i),k),             &
                     Varout(oTheta(i),k),                               &
                     Varout(oSI(i),k), Varout(oLno3(i),k))

     ! Total phytoplankton biomass
     PHYtot=PHYtot + Vars(iPHY(i),k)

     ! Total palatable prey biomass
     PHYtot1=PHYtot1 + Vars(iPHY(i),k)*RL(i)

     ! The denominator of the second term in Eq. 6 (for prey switching)
     PHYtot2=PHYtot2+Vars(iPHY(i),k)**alphaG*RL(i)

     pp_PN  =Vars(iPHY(i),k)*Varout(omuNet(i),k)+pp_PN

     ! Chl in this size class
     dx     =Vars(iPHY(i),k)/Varout(oQN(i),k)   *Varout(oTheta(i),k)
     ppC_PN =Vars(iPHY(i),k)*muNet/Varout(oQN(i),k)+ppC_PN
     CHLtot =CHLtot+dx
     Varout(oCHL(i),k)=dx 
   enddo

   ! save total phytoplankton biomass
   Varout(oPHYt,k) = PHYtot
   Varout(oCHLt,k) = CHLtot
   Varout(oPON, k) = PHYtot + ZOO + DET

   ! save phytoplankton average growth rate 
   Varout(omuAvg,k)= pp_PN/PHYtot/dtdays

   ! save total NPP (carbon-based, unit: µg C L-1 d-1)
   Varout(oPPt,k)  = ppC_PN/dtdays*12d0
!---------------------------------------------------------------
!! ZOOplankton section:
   tf_z = TEMPBOL(Ez,Temp(k))

 ! The grazing dependence on total palatable prey (dimensionless; Eq. 7)
   gbar = grazing(grazing_formulation,kp,PHYtot1)

 !Zooplankton total ingestion rate
   INGES = tf_z*dtdays*gmax*gbar

 !Zooplankton excretion rate (-> DOM)
   RES  = INGES*(1d0-GGE-unass)

 !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass
    
! Grazing rate on PHY each size class (specific to N-based Phy biomass, unit: d-1) (Eq. 12)

   ! Calculate the specific grazing rate for each size class
   do i = 1, NPHY
     ! Eq. 9 in Smith & Sergio
     Varout(oGraz(i),k)=(INGES*ZOO/PHYtot2)*Vars(iPHY(i),k)**(alphaG-1d0)*RL(i)
     Varout(oPHY(i),k) =Vars(iPHY(i),k)*(1d0+Varout(omuNet(i),k)-Varout(oGraz(i),k))
   enddo

!!End of zooplankton section
!=============================================================
!! Solve ODE functions:
  Zmort = ZOO*ZOO*dtdays* mz *tf_z  !Mortality term for ZOO
 
  ! For production/destruction matrix:
  pp_ND = dtdays* RDN *DET*tf_z   
  pp_NZ = ZOO*RES        
  pp_DZ = ZOO*EGES+Zmort 
  pp_ZP = ZOO*INGES      
 
  Varout(oDET,k) = DET+pp_DZ-pp_ND
  Varout(oNO3,k) = NO3+pp_ND+pp_NZ-pp_PN
  Varout(oZOO,k) = ZOO+pp_ZP-pp_NZ-pp_DZ
  Varout(oZ2N,k) = pp_NZ/dtdays
  Varout(oD2N,k) = pp_ND/dtdays

! Correct the unit for muNet:
  do i = 1, NPHY
     Varout(omuNet(i),k)=Varout(omuNet(i),k)/dtdays
     Varout(oGraz(i),k) =Varout(oGraz(i), k)/dtdays
  enddo
ENDDO
return
END SUBROUTINE
