SUBROUTINE Geider_Droop
    !Add phyto. Carbon into existing Geider model
use bio_MOD
implicit none
integer :: k

! Half saturation constant for Zoo grazing
real, parameter :: Kp = 0.5  

! Regenration constant from detritus to DIN (d-1)
real, parameter :: RDN= 0.1  

real    :: rhochl, theta, CHL
real    :: rmax_T, SI, muNet
real    :: par_
real    :: NO3, PHY, ZOO, DET, DET1
real    :: mz,gmax

mz   = params(imz)
gmax = params(igmax)

if (do_IRON) KFe = params(iKFe)

DO k = 1, nlev
   tf_p= TEMPBOL(Ep, Temp(k))  
   NO3 = Vars(iNO3,    k)
   PHYC= Vars(iPHYC(1),k)
   PHY = Vars(iPHY(1), k)
   ZOO = Vars(iZOO,    k)
   CHL = Vars(iCHL(1), k)
   DET = Vars(iDET,    k)
   QN  = PHY/PHYC !Cellular N:C ratio (mol N: mol C)
   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif
   ! Loss rate of phytoplankton to detritus equals to zero

   ! Define phytoplankton Chl-to-carbon ratio (gChl/molC)
   theta = CHL/PHYC

   ! DIN uptake by phytoplankton (Ward 2017)
   rhoI = rho_m *NO3/(NO3 + KN)*(QNmax-QN)/dQN*tf_p*dtdays

   ! N limitation
   Varout(oLno3(1),k) = (QN-QNmin)/dQN
   !Maximal photosynthesis rate (regulated by QN)
   PCmax = PCref * tf_P* Varout(oLno3(1),k)

   !The light limitation index (fpar)
   SI = 1d0 - exp(-params(iaI0)*par_*theta/PCmax)

   ! Photosynthesis rate
   PC = PCmax*SI

   ! define rhochl (gChl/molC;the fraction of phytoplankton carbon production that is devoted to Chl synthesis)
   rhochl = thetamax*PC/(params(iaI0)*par_*theta)

   ! The total amount of phytoplankton grazed by zooplankton (molN;gmax is the maximal specific ingestion rate!)
   tf_z   = TEMPBOL(Ez,Temp(k))
   gbar   = grazing(grazing_formulation,kp,PHY)
   INGES  = ZOO*gmax*tf_z * dtdays * gbar
   Zmort  = ZOO*ZOO*dtdays* imz *tf_z  !Mortality term for ZOO
 
   !Zooplankton excretion rate (-> DOM)
   RES = INGES*(1d0-GGE-unass)

   !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass

   ! For production/destruction matrix:
   pp_ND = dtdays* RDN *DET*tf_z   
   pp_NZ = ZOO*RES        
   pp_DZ = ZOO*EGES+Zmort 
   pp_ZP = ZOO*INGES      
  
   DET1 = DET + pp_DZ - dtdays * RDN * tf_z
   Varout(oDET,k)    = DET1
   Varout(oSI(1),k)  = SI
   Varout(oNO3,k)    = NO3+pp_ND+pp_NZ - PHYC*rhoI

   ! Equation for PHY. nitrogen
   Varout(oPHY(1),k) = PHYC*rhoI - pp_ZP

   ! Equation for PHY. carbon
   Varout(oPHYC(1),k)= PHYC*PC*dtdays - pp_ZP/QN

   Varout(omuNet(1), k) = PC
   Varout(oPPt,k)       = PHY*muNet/params(iQ0N)/dtdays
   Varout(oGraz(1), k)  = pp_ZP/PHY/dtdays
   Varout(oZOO,k) = ZOO+pp_ZP-(EGES+ZOO*dtdays*mz*tf_z+RES)
   Varout(oZ2N,k) = pp_NZ/dtdays
   Varout(oD2N,k) = pp_ND/dtdays

   Varout(oCHLt,k) = CHL + (rhochl*PC*PHYC*dtdays - pp_ZP/QN)*theta

   Varout(oCHL(1),k)   = Varout(oCHLt,k)
   Varout(oQN(1),k)    = QN
   Varout(oTheta(1),k) = theta

enddo
End subroutine Geider_Droop

subroutine PHY_GeiderDroop(Temp, PAR, NO3, QN, theta, Lno3, SI, muC)
use bio_MOD, only: params, irhom
implicit none

! Minimal and maximal N:C ratio
real, parameter :: QNmin = 0.06, QNmax = 0.18, dQN = 0.12  
  
   ! DIN uptake by phytoplankton (Ward 2017)
   rhoI = rho_m *NO3/(NO3 + KN)*(QNmax-QN)/dQN*tf_p*dtdays

   ! N limitation
   Varout(oLno3(1),k) = (QN-QNmin)/dQN
   !Maximal photosynthesis rate (regulated by QN)
   PCmax = PCref * tf_P* Varout(oLno3(1),k)

   !The light limitation index (fpar)
   SI = 1d0 - exp(-params(iaI0)*par_*theta/PCmax)

   ! Photosynthesis rate
   PC = PCmax*SI

   ! define rhochl (gChl/molC;the fraction of phytoplankton carbon production that is devoted to Chl synthesis)
   rhochl = thetamax*PC/(params(iaI0)*par_*theta)

end subroutine PHY_GeiderDroop
