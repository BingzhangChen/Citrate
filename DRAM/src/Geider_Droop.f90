SUBROUTINE Geider_Droop
    !Add phyto. Carbon into existing Geider model
USE bio_MOD
implicit none
integer :: k

! Half saturation constant for Zoo grazing
real, parameter :: Kp = 0.5  

! Regenration constant from detritus to DIN (d-1)
real, parameter :: RDN= 0.1  

real    :: rhochl, theta,QN, CHL
real    :: SI,Lno3, PC
real    :: par_, 
real    :: NO3, PHY, PHYC, ZOO, DET
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

   ! Save data:
   Varout(oQN(1),k)    = QN
   Varout(oTheta(1),k) = theta

   ! Main phyto. subroutine
   call PHY_GeiderDroop(Temp(k), par_, NO3, QN, theta, &
                        Lno3, SI, PC, rhoI, rhochl)

   ! N limitation
   Varout(oLno3(1),k) = Lno3
   Varout(oSI(1),  k) = SI

   ! The total amount of phytoplankton grazed by zooplankton (molN;gmax is the maximal specific ingestion rate!)
   tf_z   = TEMPBOL(Ez,Temp(k))
   gbar   = grazing(grazing_formulation,Kp,PHY)
   INGES  = ZOO*gmax*tf_z * gbar
   Zmort  = ZOO*ZOO*mz*tf_z  !Mortality term for ZOO
 
   !Zooplankton excretion rate (-> DOM)
   RES = INGES*(1d0-GGE-unass)

   !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass

   ! For production/destruction matrix:
   pp_ND = RDN*DET*tf_z   
   pp_NZ = ZOO*RES        
   pp_DZ = ZOO*EGES+Zmort 
   pp_ZP = ZOO*INGES      
  
   Varout(oDET,k)   = DET+dtdays*(pp_DZ - pp_ND)
   Varout(oNO3,k)   = NO3+dtdays*(pp_ND+pp_NZ-PHYC*rhoI)

   ! Equation for PHY. nitrogen (decouple nutrient uptake and photosynthesis)
   Varout(oPHY(1),k)= PHY+dtdays*(PHYC*rhoI-pp_ZP)
   Varout(oZOO,k)   = ZOO+dtdays*(pp_ZP-pp_DZ-pp_NZ)

   ! Equation for PHY. carbon
   Varout(oPHYC(1),k)=PHYC+dtdays*(PHYC*PC - pp_ZP/QN)

   Varout(omuNet(1), k) = PC
   Varout(oGraz(1), k)  = pp_ZP/PHY
   Varout(oZ2N,k) = pp_NZ
   Varout(oD2N,k) = pp_ND

   Varout(oCHLt,k)  =CHL + dtdays*(rhochl*PC*PHYC - pp_ZP/QN)*theta

   Varout(oCHL(1),k)=Varout(oCHLt,k)

   ! Calculate NPP based on incubation light
   call PHY_GeiderDroop(Temp(k), PAR(k), NO3, QN, theta, &
                        Lno3, SI, PC, rhoI, rhochl)
   Varout(oPPt,k)=PHYC*PC
enddo
End subroutine Geider_Droop

subroutine PHY_GeiderDroop(Temp, PAR, NO3, QN, theta, Lno3, SI, muC, rhoI, rhochl)
use bio_MOD, only: params, irhom, TEMPBOL, Ep, iaI0, thetamax

implicit none

real, intent(in)  :: Temp, PAR, NO3, QN, theta

! Minimal and maximal N:C ratio
real, parameter   :: QNmin = 0.06, QNmax = 0.18, dQN = 0.12  

! Nitrogen uptake rate
real, intent(out) :: rhoI 

! Carbon based photosynthetical rate
real, intent(out) :: muC

! Indices for light and nutrient limitation
real, intent(out) :: Lno3, SI, rhochl

real :: rho_m  !Maximal (reference) nitrogen uptake rate

!Half-saturation constant for N uptake
real, parameter :: KN = 0.2    

rho_m= params(irhom)

!Temperature coefficient
tf_p = TEMPBOL(Ep, Temp(k))  

! DIN uptake by phytoplankton (Ward 2017)
rhoI = rho_m *NO3/(NO3 + KN)*(QNmax-QN)/dQN*tf_p

! N limitation
Lno3 = (QN-QNmin)/dQN

!Maximal photosynthesis rate (regulated by QN)
PCmax = PCref * tf_P*Lno3

!The light limitation index (fpar)
SI = 1d0 - exp(-params(iaI0)*par_*theta/PCmax)

! Photosynthesis rate (d-1)
muC = PCmax*SI

! define rhochl (gChl/molC;the fraction of phytoplankton carbon production that is devoted to Chl synthesis)
rhochl = thetamax*muC/(params(iaI0)*par_*theta)

end subroutine PHY_GeiderDroop
