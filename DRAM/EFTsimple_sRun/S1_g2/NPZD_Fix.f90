subroutine NPZD_Fix
use bio_MOD
implicit none
integer :: k
real    :: NO3, PHY, ZOO, DET, DET1
real    :: par_, Ep, Ez, rmax_T, SI, muNet

Ep = params(iEp)
Ez = params(iEz)

DO k = 1, nlev
   NO3 = Vars(iNO3,   k)
   PHY = Vars(iPHY(1),k)
   ZOO = Vars(iZOO,   k)
   DET = Vars(iDET,   k)

   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif
  ! Loss rate of phytoplankton to detritus equals to zero
  ! The maximal growth rate (rmax_T) under temperature tC 
   rmax_T = params(imu0)* TEMPBOL(Ep,Temp(k)) * dtdays

   !The light limitation index (fpar)
   ! Unit of aI0: (W m-2)-1, different from other models
   SI = 1d0 - exp(-params(iaI0) * par_ )

   ! Growth rate (d-1) is the product function of temperature, light, and nutrient   
   Varout(oLno3(1),k) = NO3/(NO3 + params(iKn))
   muNet = rmax_T * Varout(oLno3(1),k) * SI

   ! The total amount of phytoplankton grazed by zooplankton (molN;gmax is the maximal specific ingestion rate!)
   tf_z   = TEMPBOL(Ez,Temp(k))

   gbar   = grazing(grazing_formulation,params(ikp),PHY)

   INGES  = ZOO*params(igmax)*tf_z * dtdays * gbar

   Zmort = ZOO*ZOO*dtdays* params(imz) *tf_z  !Mortality term for ZOO
 
   !Zooplankton excretion rate (-> DOM)
   RES = INGES*(1d0-GGE-unass)

   !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass

  ! For production/destruction matrix:
  pp_ND = dtdays* params(irDN) *DET*tf_z   
  pp_NZ = ZOO*RES        
  pp_DZ = ZOO*EGES+Zmort 
  pp_ZP = ZOO*INGES      
  
  DET1 = (DET + pp_DZ)/(1d0 + dtdays * params(irDN)*tf_z)
  Varout(oDET,k)    = DET1
  Varout(oSI(1),k)  = SI
  Varout(oNO3,k)    = (NO3+pp_ND+pp_NZ)/(1d0+ PHY*muNet/NO3)
  Varout(oPHY(1),k) = PHY*(1d0 + muNet)/(1d0 + pp_ZP/PHY)

  Varout(omuNet(1), k) = muNet/dtdays
  Varout(oPPt,k)       = PHY*muNet/params(iQ0N)/dtdays
  Varout(oGraz(1), k)  = pp_ZP/PHY/dtdays
  Varout(oZOO,k)       = (ZOO+pp_ZP)/(1d0+EGES+ZOO*dtdays*params(imz)*tf_z+RES)
  Varout(oZ2N,k)       = pp_NZ/dtdays
  Varout(oD2N,k)       = pp_ND/dtdays
  Varout(oCHLt,k)      = PHY/params(iQ0N) * params(itheta)
Enddo
return
End subroutine NPZD_Fix
