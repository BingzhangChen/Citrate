subroutine Biology
USE PARAM_MOD
IMPLICIT NONE
integer :: k
real    :: NO3, PHY, ZOO, DET
real    :: par_,  muNet
real    :: Kp, gmax, RDN,mz
real    :: QN, bI0=0d0, KFe_, QP
real, parameter :: Qpmin = 0.05/16.
real, parameter :: Qnmin = 0.05
real, parameter :: KPO4  = 0.5/16.

Kp  = 0.5
gmax= exp(params(igmax)) * exp(params(imu0))
RDN = 0.1
mz  = exp(params(imz))* exp(params(imu0))

if (do_IRON) KFe = params(iKFe)

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

   ! Calculate phytoplankton growth rate if not considering mixing
   CALL MONOD(Temp(k), PAR(k), NO3, 1.,exp(params(imu0)),Qnmin, Qpmin,  &
              exp(params(iaI0_C)), &
              bI0, exp(params(iKN)),KPO4, DFe(k), KFe_,                        &
              muNet, QN, QP,Varout(oTheta,k),                          &
              Varout(oSI,k), Varout(oLno3,k))

   Varout(oPPt,k) = Vars(iPHY(1),k)*muNet/QN*12d0   !NPP Unit: ug C L-1

   if (k .gt. N_MLD) then  !Within surface mixed layer
      ! Calculate phytoplankton growth rate, theta, and QN based on environmental conditions
      CALL MONOD(Temp(k), PAR_, NO3, 1.,exp(params(imu0)),Qnmin, Qpmin,  &
                 exp(params(iaI0_C)), &
                 bI0, exp(params(iKN)),KPO4, DFe(k), KFe_,                        &
                 muNet, Varout(oQN,k),QP,Varout(oTheta,k),             &
                 Varout(oSI,k), Varout(oLno3,k))
   else
      Varout(oQN,k) = QN
   endif

   ! Correct the unit of growth rate
   muNet = muNet*dtdays

   ! The total amount of phytoplankton grazed by zooplankton (molN;gmax is the maximal specific ingestion rate!)
   tf_z  = TEMPBOL(Ez,Temp(k))  ! Temperature effect of zooplankton
   gbar  = grazing(grazing_formulation,Kp,PHY) !Functional response of zooplankton
   INGES = gmax*tf_z*dtdays * gbar   !Specific ingestion rate of zooplankton
   Zmort = ZOO*ZOO*dtdays* mz *tf_z  !Mortality term for ZOO
 
   !Zooplankton excretion rate (-> DOM)
   RES = INGES*(1d0-GGE-unass)

   !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass

  ! For production/destruction matrix:
  pp_ND = dtdays* RDN *DET*tf_z   
  pp_NZ = ZOO*RES        
  pp_DZ = ZOO*EGES+Zmort 
  pp_ZP = ZOO*INGES      
  Varout(oDET,k)   = (DET + pp_DZ)-pp_ND
  Varout(oNO3,k)   = (NO3+pp_ND+pp_NZ)-PHY*muNet
  Varout(oPHY(1),k)= PHY*(1d0 + muNet)-pp_ZP
  Varout(omuNet, k)= muNet/dtdays
  Varout(oGraz,  k)= pp_ZP/PHY/dtdays
  Varout(oZOO,k)   = (ZOO+pp_ZP)-pp_DZ-pp_NZ
  Varout(oZ2N,k)   = pp_NZ/dtdays
  Varout(oD2N,k)   = pp_ND/dtdays
  Varout(oCHLt,k)  = PHY/Varout(oQN,k)*Varout(otheta,k)  !Chl unit: ugchl/L

  ! Total particulate organic nitrogen (PON)
  Varout(oPON,k)=Varout(oZOO,k)+Varout(oPHY(1),k)+Varout(oDET,k)   !PON unit: umol N/L
Enddo
RETURN
END SUBROUTINE BIOLOGY
