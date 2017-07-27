SUBROUTINE FlexEFT_simple
use bio_MOD
implicit none
integer :: k
real    :: QN, mu0
real    :: par_
real    :: Lno3
real    :: NO3, PHY, ZOO, DET
real    :: pp_PN
real    :: KN = 1d0
real    :: muNet
real    :: Kp, gmax, RDN, mz
real,    external  :: WAPR

Kp  =0.5
RDN =0.1
mz  =params(imz)
gmax=params(igmax)
mu0 =params(imu0)
DO k = nlev,1,-1    ! k from surface to bottom, cannot be reversed
   DET  = Vars(iDET,k)
   NO3  = Vars(iNO3,k)
   PHY  = Vars(iPHY(1),k)
   ZOO  = Vars(iZOO,k)  !Zooplankton biomass

   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif
   Varout(oPAR_,k) = par_
! Phytoplankton section:
! Calculate phytoplankton growth rate, theta, and QN based on environmental conditions
   call EFT_phygrowth(mu0,KN,params(iA0N), params(iaI0), NO3,Temp(k),PAR(k),&
                   muNet, QN,                                         &
                   Varout(oTheta(1),k),                               &
                   Varout(oSI(1),k), Varout(oLno3(1),k))

   Varout(oPPt,k) = PHY*muNet/dtdays/Varout(oQN(1),k)*12d0

   !Call again with par_
   call EFT_phygrowth(mu0,KN,params(iA0N), params(iaI0), NO3,Temp(k),par_,&
                   muNet, Varout(oQN(1),k),                               &
                   Varout(oTheta(1),k),                               &
                   Varout(oSI(1),k), Varout(oLno3(1),k))
   
   !Save net growth rate
   Varout(omuNet(1),k) = muNet/dtdays
    
   !Phytoplankton sinking rate at the average size
   ! Varout(ow_p(1),  k) = abs(w_p0)  !Positive
!---------------------------------------------------------------
!! ZOOplankton section:
   tf_z = TEMPBOL(Ez,Temp(k))

 ! The grazing dependence on total prey (dimensionless)
   gbar = grazing(grazing_formulation,Kp,PHY)

 !Zooplankton total ingestion rate
   INGES = tf_z*dtdays*gmax*gbar

 !Zooplankton excretion rate (-> DOM)
   RES  = INGES*(1d0-GGE-unass)

 !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass
    
! Grazing rate on PHY(specific to N-based Phy biomass, unit: d-1)
   ! Calculate the specific grazing rate for PHY
   Varout(oGraz(1),k) = INGES*ZOO/PHY/dtdays
   
   Varout(oPHY(1),k)  = PHY*(1d0+muNet)-INGES*ZOO
!!End of zooplankton section
!=============================================================
!! Solve ODE functions:
  Zmort = ZOO*ZOO*dtdays*mz*tf_z  !Mortality term for ZOO
 
  ! For production/destruction matrix:
  pp_PN=PHY*muNet
  pp_ND=dtdays* RDN*DET*tf_z   
  pp_NZ=ZOO*RES        
  pp_DZ=ZOO*EGES+Zmort 
  pp_ZP=ZOO*INGES      
  
  Varout(oDET,k) = (DET+pp_DZ)-pp_ND
  Varout(oNO3,k) = (NO3+pp_ND+pp_NZ)-pp_PN
  Varout(oZOO,k) = (ZOO+pp_ZP)-pp_DZ-pp_NZ
  Varout(oZ2N,k) = pp_NZ/dtdays
  Varout(oD2N,k) = pp_ND/dtdays
  Varout(oCHLt,k)= Vars(iPHY(1),k)/Varout(oQN(1),k)*Varout(oTheta(1),k)
  Varout(oPON, k)= PHY+ZOO+DET
  Varout(oCHL(1),k)=Varout(oCHLt,k)
ENDDO
End subroutine FlexEFT_simple

