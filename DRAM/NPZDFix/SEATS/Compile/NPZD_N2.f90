subroutine NPZD_N2  ! Add N2 fixation, use P as the limiting element 
use bio_MOD
implicit none
integer :: k
real    :: NO3,PO4,PHY, ZOO, DETn,DETp,DIA
real    :: par_, muNet,muDIA
real    :: Kp, gmax, mz,KPO4
real    :: QN, Qp, bI0, KFe_, N2P_phy, theta_DIA
real    :: pp_NDn,pp_PDp,Res_DIA,Mort_DIA
real, parameter :: N2P_dia = 45. ! N:P ratio of diazotrophs
real, parameter :: N2C_dia = 0.2  ! N:C ratio of diazotrophs
real, parameter :: P2C_dia = N2C_dia/N2P_dia ! P:C ratio of diazotrophs
real, parameter :: Q0p     = 3D-3 ! Subsistence P quota of ordinary phyto.
real, parameter :: Lnifdet = 0.05 ! Diazotroph mortality   rate
real, parameter :: RDN     = 0.1  ! Regeneration rate of detritus to DIN
real            :: Lnifdin = 0.17*16 

Lnifdin=params(iLnifp)  ! Diazotroph respiration rate
Kp  = params(iKPHY)
gmax= params(igmax)
mz  = params(imz)
!bI0 = 10**params(ibI0B)
bI0 = 0.
KPO4= params(iKN)/16d0
if (do_IRON) KFe = params(iKFe)

DO k = 1, nlev
   tf_p= TEMPBOL(Ep,Temp(k))
   NO3 = Vars(iNO3,   k)
   PO4 = Vars(iPO4,   k)
   PHY = Vars(iPHY(1),k)
   ZOO = Vars(iZOO,   k)
   DETn= Vars(iDET,   k)  ! Detritus N
   DETp= Vars(iDETp,  k)  ! Detritus P

   ! Add diazotrophs:
   DIA = Vars(iDIA,   k) 

   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif
   Varout(oPAR_,k)=par_
   ! Calculate nondiazotroph growth rate if not considering mixing
   CALL MONOD(Temp(k), PAR(k), NO3,PO4,params(imu0),params(iQ0N),   &
               Q0p, params(iaI0_C),                                 &
               bI0, params(iKN), KPO4, DFe(k), KFe_,         &
               muNet, QN, QP, Varout(oTheta(1),k),                  &
               Varout(oSI(1),k), Varout(oLno3(1),k))

   ! Calculate diazotroph growth rate if not considering mixing
   call nif_growth(PO4,PAR(k),Temp(k),wstr0(1),muDIA,theta_DIA)

   ! Unit of PHY: mM P
   Varout(oPPt,k)=(PHY*muNet/QP + DIA*muDIA/P2C_dia)*12d0

   ! Calculate nondiazo growth rate, theta, and QN based on environmental conditions
   CALL MONOD(Temp(k), par_, NO3, PO4, params(imu0),params(iQ0N),  &
        Q0p,  params(iaI0_C), &
        bI0,  params(iKN), KPO4, DFe(k), KFe_,              &
        muNet,Varout(oQN(1),k),Varout(oQP(1),k),Varout(oTheta(1),k), &
        Varout(oSI(1),k), Varout(oLno3(1),k))

   N2P_phy = Varout(oQN(1),k)/Varout(oQP(1),k)  !N:P of normal phytoplankton

   ! Calculate diazotroph growth rate if not considering mixing
   call nif_growth(PO4,par_,Temp(k),wstr0(1),muDIA,theta_DIA)

   ! Correct the unit of growth rate
   muNet = muNet*dtdays

   Varout(oDIAu,k)=muDIA
   muDIA = muDIA*dtdays

   ! The total amount of phytoplankton grazed by zooplankton (molN;gmax is the maximal specific ingestion rate!)
   tf_z  = TEMPBOL(Ez,Temp(k))  ! Temperature effect of zooplankton
   gbar  = grazing(grazing_formulation,Kp,PHY) !Functional response of zooplankton
   INGES = gmax*tf_z*dtdays * gbar   !Specific ingestion rate of zooplankton
   Zmort = ZOO*ZOO*dtdays*mz*tf_z  !Mortality term for ZOO
 
   !Zooplankton excretion rate (-> DOM)
   RES = INGES*(1d0-GGE-unass)

   !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass

  ! For production/destruction matrix:
  pp_NDn = dtdays*RDN*DETn*tf_z   
  pp_PDp = dtdays*params(iRDN_P)*DETp*tf_z
  pp_NZ  = ZOO*RES        
  pp_DZ  = ZOO*EGES+Zmort 
  pp_ZP  = ZOO*INGES      
  
  ! Respiration of diazotrophs (DIA -> DIP):
  Res_DIA  = Lnifdin*tf_p*dtdays*DIA**2

  ! Mortality of diazotrophs (DIA -> POP):
  Mort_DIA = Lnifdet*tf_p*dtdays*DIA

  ! N as the unit:
  Varout(oDET,k)  = (DETn + pp_DZ*N2P_phy) - pp_NDn + Mort_DIA*N2P_dia

  ! P as the unit
  Varout(oDETp,k) = (DETp + pp_DZ)- pp_PDp + Mort_DIA

  Varout(oNO3,k)  = (NO3+pp_NDn+pp_NZ*N2P_phy+Res_DIA*N2P_dia) &
                  - PHY*muNet*N2P_phy   ! N as the unit

  Varout(oPO4,k)  = (PO4+pp_PDp+pp_NZ        +Res_DIA) &
                  - PHY*muNet-DIA*muDIA ! P as the unit

  Varout(oPHY(1),k)    = PHY*(1d0 + muNet)-pp_ZP
  Varout(oZOO,k)       = (ZOO+pp_ZP)-pp_DZ-pp_NZ
  Varout(oDIA,k)       = DIA*(1d0 + muDIA)-Mort_DIA-Res_DIA

  Varout(omuNet(1), k) = muNet/dtdays
  Varout(oGraz(1), k)  = pp_ZP/PHY/dtdays
  Varout(oZ2N,k)   = pp_NZ/dtdays
  Varout(oD2N,k)   = pp_ND/dtdays
  Varout(oCHLt,k)  = PHY/Varout(oQP(1),k)*Varout(otheta(1),k) &
                   + DIA/P2C_dia*theta_DIA
  Varout(oCHL(1),k)= Varout(oCHLt,k)

  ! Total particulate organic nitrogen (PON)
  Varout(oPON,k)=(Varout(oZOO,k)+Varout(oPHY(1),k))*N2P_phy &
            +     Varout(oDET,k)+Varout(oDIA,k)*N2P_dia

  ! Total POP:
  Varout(oPOP,k)=Varout(oZOO,k)+Varout(oPHY(1),k) + Varout(oDETp,k) + Varout(oDIA,k)
Enddo
return
End subroutine NPZD_N2
