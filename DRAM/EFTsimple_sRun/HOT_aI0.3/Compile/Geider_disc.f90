subroutine Geider_DISC
use bio_MOD
implicit none
integer :: k, i, j
real    :: CHLtot, theta
real    :: par_, PHYtot, PHYtot2, pp_PN 
real    :: DET, ZOO, NO3, DET1, CHL
real    :: rhochl
real    :: CHLs(4)= 0d0  ! Size fractionated Chl

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
   CHLtot  =0d0  ! Total CHL A
   CHLs(:) =0d0  ! Size fractionated CHL
  
!! Phytoplankton section:
   do i = 1, NPHY

     ! Calculate phytoplankton Chl-to-carbon ratio (gChl/molC)
     theta   = Vars(iCHL(i),k)/( Vars(iPHY(i),k) / params(iQ0N))
     Varout(oTheta(i),k) = theta

     ! PMU_ is log(10*pi/6*ESD**3)
     PMU_(i) = PMU_(i) - log(1d1) ! Restore to original value

     call GeiderPhy_size(PMU_(i),Vars(iNO3,k), Temp(k), par_, theta,    &
             Varout(omuNet(i),k), Varout(oSI(i),k), Varout(oLno3(i),k))

     PMU_(i) = PMU_(i) + log(1d1) ! Restore to positive values
     PHYtot  = PHYtot +Vars(iPHY(i),k)
     PHYtot2 = PHYtot2+Vars(iPHY(i),k)**params(ialphaG)
     pp_PN   = Vars(iPHY(i),k)*Varout(omuNet(i),k) + pp_PN
     CHLtot  = CHLtot + Vars(iCHL(i),k)
 
     do j = 1, 4
        CHLs(j) = CHLs(j) + Vars(iCHL(i),k) * wtCHL(j,i)
     enddo
   enddo


   ! save total phytoplankton biomass
   Varout(oPHYt,k) = PHYtot
   Varout(oCHLt,k) = CHLtot
   
   ! save total NPP (carbon-based)
   Varout(oPPt,k)  = pp_PN/dtdays/params(iQ0N)

   ! save size fractionated Chl
   do j = 1,4
      Varout(oCHLs(j),k) = CHLs(j)/Varout(oCHLt,k)  ! Percentage of one size class
   enddo
!---------------------------------------------------------------
!! ZOOplankton section:
   tf_z = TEMPBOL(params(iEz),Temp(k))

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
       *Vars(iPHY(i),k)**(params(ialphaG)-1d0)
     
     Varout(oPHY(i),k) = max( Vars(iPHY(i),k)*(1d0+Varout(omuNet(i),k)) &
       /(1d0 + Varout(oGraz(i),k)), 1D-20)

   ! define rhochl (gChl/molC;the fraction of phytoplankton carbon production that is devoted to Chl synthesis)
     rhochl = thetm*Varout(omuNet(i),k)/(params(iaI0)*dtdays*par_*Varout(otheta(i),k))

     ! Calculate Chl:
     Varout(oCHL(i),k) = (Vars(iCHL(i),k) + rhochl*Varout(omuNet(i),k)  &
       * Vars(iPHY(i),k)/params(iQ0N) ) / (1d0 + pp_ZP/Vars(iPHY(i),k))

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

end subroutine Geider_DISC
