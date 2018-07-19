SUBROUTINE Geider_simple
use bio_MOD
implicit none
integer :: k
real    :: rhochl, theta, CHL
real    :: rmax_T, SI, muNet
real    :: par_
real    :: NO3, PHY, ZOO, DET, DET1
real    :: mz,RDN,gmax,Kp

!Ep = params(iEp)
!Ez = params(iEz)
Ep=0.5
Ez=0.6
mz=.15
RDN=.1
gmax=1.0
Kp=.5

if (do_IRON) KFe = params(iKFe)

DO k = 1, nlev
   NO3 = Vars(iNO3,   k)
   PHY = Vars(iPHY(1),k)
   ZOO = Vars(iZOO,   k)
   CHL = Vars(iCHL(1),k)
   DET = Vars(iDET,   k)

   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif
   ! Loss rate of phytoplankton to detritus equals to zero

   ! Define phytoplankton Chl-to-carbon ratio (gChl/molC)
   theta = CHL/(PHY/params(iQ0N))

   ! The maximal growth rate (rmax_T) under temperature tC 
   rmax_T = params(imu0)* TEMPBOL(Ep,Temp(k)) * dtdays

   !The light limitation index (fpar)
   SI = 1d0 - exp(-params(iaI0) * dtdays *par_ * theta/rmax_T)

   if (do_IRON) rmax_T = rmax_T *  DFe(k) /(DFe(k) + KFe)

   ! Growth rate (d-1) is the product function of temperature, light, and nutrient   
   Varout(oLno3(1),k) = NO3/(NO3 + params(iKn))
   muNet = rmax_T * Varout(oLno3(1),k)*SI

   ! define rhochl (gChl/molC;the fraction of phytoplankton carbon production that is devoted to Chl synthesis)
   rhochl = thetm*muNet/(params(iaI0)*dtdays*par_*theta)

   ! The total amount of phytoplankton grazed by zooplankton (molN;gmax is the maximal specific ingestion rate!)
   tf_z   = TEMPBOL(Ez,Temp(k))
   gbar   = grazing(grazing_formulation,kp,PHY)
   INGES  = ZOO*gmax*tf_z * dtdays * gbar

   Zmort = ZOO*ZOO*dtdays* imz *tf_z  !Mortality term for ZOO
 
   !Zooplankton excretion rate (-> DOM)
   RES = INGES*(1d0-GGE-unass)

   !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass

  ! For production/destruction matrix:
  pp_ND = dtdays* RDN *DET*tf_z   
  pp_NZ = ZOO*RES        
  pp_DZ = ZOO*EGES+Zmort 
  pp_ZP = ZOO*INGES      
  
  DET1 = (DET + pp_DZ)/(1d0 + dtdays * RDN *tf_z)
  Varout(oDET,k)    = DET1
  Varout(oSI(1),k)  = SI
  Varout(oNO3,k)    = (NO3+pp_ND+pp_NZ)/(1d0+ PHY*muNet/NO3)
  Varout(oPHY(1),k) = PHY*(1d0 + muNet)/(1d0 + pp_ZP/PHY)

  Varout(omuNet(1), k) = muNet/dtdays
  Varout(oPPt,k)       = PHY*muNet/params(iQ0N)/dtdays
  Varout(oGraz(1), k)  = pp_ZP/PHY/dtdays
  Varout(oZOO,k) = (ZOO+pp_ZP)/(1d0+EGES+ZOO*dtdays*mz*tf_z+RES)
  Varout(oZ2N,k) = pp_NZ/dtdays
  Varout(oD2N,k) = pp_ND/dtdays

  Varout(oCHLt,k) = (CHL + rhochl*muNet*PHY/params(iQ0N)) &
                 / (1d0 + pp_ZP/PHY)

  Varout(oCHL(1),k)   = Varout(oCHLt,k)
  Varout(oQN(1),k)    = params(iQ0N)
  Varout(oTheta(1),k) = Varout(oCHLt,k)/(Varout(oPHY(1),k)/params(iQ0N))

  !Varout(ow_p(1),k) = w_p0   !Should not be converted by dtdays
  
enddo
End subroutine Geider_simple
