subroutine NPZD_DISC
use bio_MOD
implicit none
integer :: k, i, j, itheta=1
real    :: par_,PHYtot, PHYtot2, pp_PN 
real    :: DET, ZOO, NO3, DET1
real    :: dx, alphaG
real    :: CHLs(4)= 0d0  ! Size fractionated Chl

!Ep = params(iEp)
Ep=0.5
Ez=0.6
!Ez = params(iEz)

if (kill_the_winner) then
!   alphaG=params(ialphaG)
   alphaG=1.1
else
   alphaG=1d0
endif
DO k = nlev,1,-1

   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif

   DET  = Vars(iDET,k)
   ZOO  = Vars(iZOO,k)  !Zooplankton biomass
   NO3  = Vars(iNO3,k)

   PHYtot  =0d0  ! Calculate total PHY biomass
   PHYtot2 =0d0  ! total P**alphaG
   pp_PN   =0d0  ! N-based NPP
   CHLs(:) =0d0  ! Size fractionated CHL
  
!! Phytoplankton section:
   do i = 1, NPHY

     call NPZDPhy_size(PMU_(i),Vars(iNO3,k), Temp(k), par_, Varout(omuNet(i),k), &
             Varout(oSI(i),k), Varout(oLno3(i),k))

     PHYtot  = PHYtot +Vars(iPHY(i),k)
     PHYtot2 = PHYtot2+Vars(iPHY(i),k)**alphaG
     pp_PN   = Vars(iPHY(i),k)*Varout(omuNet(i),k) + pp_PN
    
     ! Chl in this size class
     dx      =Vars(iPHY(i),k)/ params(iQ0N) * params(itheta)
     Varout(oCHL(i),k)=dx
     Varout(otheta(i),k)=params(itheta)
     Varout(oQN(i), k)=params(iQ0N)
     do j = 1, 4
        CHLs(j) = CHLs(j) + dx*wtCHL(j,i)
     enddo
  
   enddo


   ! save total phytoplankton biomass
   Varout(oPHYt,k) = PHYtot
   Varout(oCHLt,k) = PHYtot / params(iQ0N) * params(itheta)
   
   ! save total NPP (carbon-based)
   Varout(oPPt,k)  = pp_PN/dtdays/params(iQ0N)

   ! save size fractionated Chl
   do j = 1,4
      Varout(oCHLs(j),k) = CHLs(j)/Varout(oCHLt,k)  ! Percentage of one size class
   enddo
!---------------------------------------------------------------
!! ZOOplankton section:
   tf_z = TEMPBOL(Ez,Temp(k))

 ! The grazing dependence on total prey (dimensionless)
   gbar = grazing(grazing_formulation,params(ikp),PHYtot)
 !Zooplankton total ingestion rate
   INGES = tf_z*dtdays*params(igmax)*gbar

 !Zooplankton excretion rate (-> DOM)
   RES  = INGES*(1d0-GGE-unass)

 !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass
    
! Grazing rate on PHY each size class (specific to N-based Phy biomass, unit: d-1) (Eq. 12)

   ! Calculate the specific grazing rate for each size class
   do i = 1, NPHY
    
     ! Eq. 10 in Smith & Sergio
     Varout(oGraz(i),k) = (INGES*ZOO/PHYtot2)                       &
       *Vars(iPHY(i),k)**(alphaG - 1d0)
     
     Varout(oPHY(i),k) = max( Vars(iPHY(i),k)*(1d0+Varout(omuNet(i),k)) &
       /(1d0 + Varout(oGraz(i),k)), 1D-20)

   enddo

!!End of zooplankton section
!=============================================================
!! Solve ODE functions:
  Zmort = ZOO*ZOO*dtdays* params(imz) *tf_z  !Mortality term for ZOO
 
  ! For production/destruction matrix:
  pp_ND = dtdays* params(irDN) *DET*tf_z   
  pp_NZ = ZOO*RES        
  pp_DZ = ZOO*EGES+Zmort 
  pp_ZP = ZOO*INGES      
  DET1  = (DET+pp_DZ)/(1d0 + dtdays * params(irDN)*tf_z)
 
  Varout(oDET,k) = DET1
  Varout(oNO3,k) = max( (NO3+pp_ND+pp_NZ)/(1d0+pp_PN/NO3), 1D-20)

  Varout(oZOO,k) = (ZOO+pp_ZP)/(1d0   &
                 + EGES+ZOO*dtdays*params(imz)*tf_z+RES)
    
  Varout(oZ2N,k) = pp_NZ/dtdays
  Varout(oD2N,k) = pp_ND/dtdays
ENDDO

end subroutine NPZD_DISC
