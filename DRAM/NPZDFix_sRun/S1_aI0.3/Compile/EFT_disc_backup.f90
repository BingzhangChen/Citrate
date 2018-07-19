SUBROUTINE FLEXEFT_DISC
use bio_MOD
implicit none
integer :: i,j,k
real :: dx,par_
real :: NO3, ZOO, DET, DET1
real :: PHYtot = 0d0
real :: PHYtot2= 0d0
real :: ppC_PN = 0d0
real :: pp_PN  = 0d0
real :: CHLtot = 0d0
real :: CHLs(4)= 0d0  ! Size fractionated Chl
real :: alphaG,Kp,gmax,RDN,mz
real, parameter :: PMU0=log(10.)

Ez  = 0.6
Ep  = 0.5
Ez  = 0.6
Kp  = 0.5
gmax= 1.0
RDN = 0.1
mz  = 0.15

if (kill_the_winner) then
   alphaG=1.1
else
   alphaG=1d0
endif

DO k = nlev,1,-1
   PHYtot  =0d0  ! Calculate total PHY biomass
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
     PHYtot =PHYtot +Vars(iPHY(i),k)
     PHYtot2=PHYtot2+Vars(iPHY(i),k)**alphaG
     pp_PN  =Vars(iPHY(i),k)*Varout(omuNet(i),k)+pp_PN

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
   Varout(oCHLt,k) = CHLtot
   ! save total NPP (carbon-based)
   Varout(oPPt,k)  = ppC_PN/dtdays

   ! save size fractionated Chl
   do j = 1,4
      Varout(oCHLs(j),k) = CHLs(j)/CHLtot  ! Percentage of one size class
   enddo
!---------------------------------------------------------------
!! ZOOplankton section:
   tf_z = TEMPBOL(Ez,Temp(k))
 ! The grazing dependence on total prey (dimensionless)
   gbar = grazing(grazing_formulation,kp,PHYtot)
 !Zooplankton total ingestion rate
   INGES= tf_z*dtdays*gmax*gbar

 !Zooplankton excretion rate (-> DOM)
   RES  = INGES*(1d0-GGE-unass)

 !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass
    
! Grazing rate on PHY each size class (specific to N-based Phy biomass, unit: d-1) (Eq. 12)

   ! Calculate the specific grazing rate for each size class
   do i=1,NPHY
    
     ! Eq. 10 in Smith & Sergio
     Varout(oGraz(i),k)=(INGES*ZOO/PHYtot2)*Vars(iPHY(i),k)**(alphaG-1d0)
     
     Varout(oPHY(i),k) = max(Vars(iPHY(i),k)*(1d0+Varout(omuNet(i),k))  &
       /(1d0 + Varout(oGraz(i),k)), 1D-10)
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
  DET1  = (DET+pp_DZ)/(1d0 + dtdays * rDN*tf_z)
 
  Varout(oDET,k) = DET1
  Varout(oNO3,k) = max( (NO3+pp_ND+pp_NZ)/(1d0+pp_PN/NO3), 1D-10)

  Varout(oZOO,k) = (ZOO+pp_ZP)/(1d0   &
                 + EGES+ZOO*dtdays*mz*tf_z+RES)
    
  Varout(oZ2N,k) = pp_NZ/dtdays
  Varout(oD2N,k) = pp_ND/dtdays
ENDDO
return
END SUBROUTINE FLEXEFT_DISC
