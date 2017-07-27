SUBROUTINE FLEXEFT_DISC
use bio_MOD
implicit none
integer :: i,j,k
real :: dx,par_
real :: NO3, ZOO, DET
real :: PHYtot = 0d0
real :: PHYtot2= 0d0
real :: ppC_PN = 0d0
real :: pp_PN  = 0d0
real :: CHLtot = 0d0
real :: CHLs(4)= 0d0  ! Size fractionated Chl
real :: alphaG,Kp,gmax,RDN,mz,PHYtot1
real,    parameter :: PMU0      = log(10.)

! The parameter that describes how feeding preference decreases with cell size
real,    parameter :: alphaGmax = 0d0

! The feeding preference at 1 µm3
real,    parameter :: PREF0     = 1d0
real               :: RL(NPHY)  ! Feeding preference on each size class

Kp   = 0.5
gmax = params(igmax)
RDN  = 0.1
mz   = params(imz)
alphaG=1.+params(ialphaG)

DO k = nlev,1,-1
   PHYtot  =0d0  ! Calculate total PHY biomass
   PHYtot1 =0d0  ! Total palatable prey
   PHYtot2 =0d0  ! total P**alphaG
   pp_PN   =0d0  ! total primary production (mmol N per d per m3)
   ppC_PN  =0d0  ! total primary production (mmol C per d per m3)
   CHLtot  =0d0  ! Total CHL A
   CHLs(:) =0d0  ! Size fractionated CHL
  
   DET  = Vars(iDET,k)
   ZOO  = Vars(iZOO,k)  !Zooplankton biomass
   NO3  = Vars(iNO3,k)

   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif
   Varout(oPAR_,k)=par_

!! Phytoplankton section:
   do i=1,NPHY
     ! PMU_ is log(10*pi/6*ESD**3)
     PMU_(i) = PMU_(i) - PMU0 ! Restore to original value

     call EFT_size(PMU_(i),Vars(iNO3,k),Temp(k), par_,                 &
          Varout(omuNet(i),k), Varout(oQN(i),k), Varout(oTheta(i),k),  &
          Varout(oSI(i),k), Varout(oLno3(i), k)  )

     if (Varout(oTheta(i),k) .le. 0.01 .and. k .lt. nlev) then
        Varout(oTheta(i), k) = Varout(oTheta(i), k+1)
     endif

     PMU_(i)=PMU_(i)+PMU0 ! Restore to positive values

    !The feeding preference of the generic zooplankton for prey i
     RL(i)=ScaleTrait(PMU_(i), PREF0, alphaGmax)

     ! Total phytoplankton biomass
     PHYtot=PHYtot + Vars(iPHY(i),k)

     ! Total palatable prey biomass
     PHYtot1=PHYtot1 + Vars(iPHY(i),k)*RL(i)

     ! The denominator of the second term in Eq. 6 (for prey switching)
     PHYtot2=PHYtot2+Vars(iPHY(i),k)**alphaG*RL(i)

     pp_PN  = Vars(iPHY(i),k)*Varout(omuNet(i),k)+pp_PN

     ! Chl in this size class
     dx     =Vars(iPHY(i),k)/Varout(oQN(i),k)   *Varout(oTheta(i),k)
     ppC_PN =Vars(iPHY(i),k)*Varout(omuNet(i),k)/Varout(oQN(i),k)+ppC_PN
     CHLtot =CHLtot+ dx

     Varout(oCHL(i),k)=dx 

     do j = 1, 4
        CHLs(j) = CHLs(j) + dx*wtCHL(j,i)
     enddo
   enddo

   ! save total phytoplankton biomass
   Varout(oPHYt,k) = PHYtot
   Varout(oPON, k) = PHYtot + ZOO + DET
   Varout(oCHLt,k) = CHLtot

   ! save total NPP (carbon-based, unit: mg C L-1 d-1)
   Varout(oPPt,k)  = ppC_PN/dtdays*12.

   ! save size fractionated Chl
   do j = 1,4
      Varout(oCHLs(j),k) = CHLs(j)/CHLtot  ! Percentage of one size class
   enddo
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
  !  DET1  = (DET+pp_DZ)/(1d0 + dtdays * rDN*tf_z)
 
  Varout(oDET,k) =DET+pp_DZ-pp_ND
  Varout(oNO3,k) =NO3+pp_ND+pp_NZ-pp_PN
  Varout(oZOO,k) =ZOO+pp_ZP-pp_NZ-pp_DZ
    
  Varout(oZ2N,k) = pp_NZ/dtdays
  Varout(oD2N,k) = pp_ND/dtdays

! Correct the unit for muNet:
  
  do i = 1, NPHY
     Varout(omuNet(i),k)=Varout(omuNet(i),k)/dtdays
     Varout(oGraz(i),k) =Varout(oGraz(i), k)/dtdays
  enddo
ENDDO
return
END SUBROUTINE FLEXEFT_DISC
