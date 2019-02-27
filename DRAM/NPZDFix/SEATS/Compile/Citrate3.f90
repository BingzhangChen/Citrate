SUBROUTINE Citrate3_MOD
use bio_MOD
use MOD_1D, only: it, nsave
implicit none

integer :: k,j,i,k1
!INPUT PARAMETERS:
real :: tC,par_,I_min,I_max,T_min,T_max
!LOCAL VARIABLES of phytoplankton:
real :: PHY1
real :: NO3, Fe, PHY, MIC, MES, DET, DETFe, DET1
real :: PMU,VAR,PMU1,VAR1,PMUPHY,VARPHY
real :: VIo,VIo1,VTo,VTo1
real :: MToPHY,VToPHY,MIoPHY,VIoPHY
real :: mu2,dmudZ2,d2mudX22,dmudX2
real :: QN  ! cell quota related variables
real :: muNet,dmuNetdl,d2muNetdl2, dmudT,d2mudL22
real :: dmudI1,d2mudI21,dmudT1,d2mudT21,dthetadT1,d2muNet_QNdx2
real :: dmudL2,d2mudZ2,theta2,dmudI,d2muDI2,d2mudT2
real :: dmuNetdl1,d2muNetdl21,d2muNet_QNdZ2
real :: SI,Lno3,Fescav
real :: muNet1, SI1, Lno31, QN1,QN2, MTo,MTo1,MIo,MIo1
real :: theta,theta1,dthetadl,d2thetadl2,dQNdL,d2QNdL2
real :: dthetadT,d2thetadT2,dthetadI,d2thetadI2
real :: L_min,L_max,dthetadZ2,d2thetadZ22,dthetadX2,d2thetadx22
real :: dthetadl1,d2thetadl21,dQNdL1,d2QNdL21
real :: dthetadl2,d2thetadl22,dQNdL2,d2QNdL22
real :: d2thetadT21,dthetadI1,d2thetadI21
real :: VTR_L,VTR_I,VTR_T,d2muNet_QNdl2

!Declarations of zooplankton:
real :: dgdlbar1,d2gdl2bar1,dgdlbar2,d2gdl2bar2
real :: INGES1,INGES2,RES1,RES2,EGES1,EGES2,Zmort2,mz2
real :: Kp1,Cf, cff, cff1  !Zooplankton variables
real :: pp_PN, PP_MIC_P, PP_MES_P,PP_MES_MIC
real :: Ptot,  CHLt, Pl, dCHL,rl
real :: pCHL(4) = 0D0

integer, parameter :: M     = 10 !discretize the normal distribution
real               :: dL    = 0d0,dI = 0., dT = 0. 
real               :: L_(M) = 0.,I_(M)=0., T_(M)=0.
real,     external :: normal
real,    parameter :: PMU0  = log(1d1)
real,    parameter :: PMU10 = log(pi/6*1d1**3)
real,    parameter :: Kp2   = .5 
real,    parameter :: gmax1 = 1.35 !Chai et al. (2002)
real,    parameter :: gmax2 =  .53 !Chai et al. (2002)
real,    parameter :: gb1   = -0.05 ! Feeding selectivity of microzoo.
real,    parameter :: gb2   =  0.02 ! Feeding selectivity of mesozoo.
real,    parameter :: unass1=  0.24 ! Fraction of unassimilated material of microzoo.
real,    parameter :: unass2=  0.31 ! Fraction of unassimilated material of mesozoo.
real,    parameter :: RDN   =  0.1
real,    parameter :: alphaK=  0.27
!------------------------------------------------------------------
VTR_L = params(iVTRL)
VTR_T = params(iVTRT)
VTR_I = params(iVTRI)
Kp1 = params(iKPHY)  ! Grazing half-saturation constant for MIC
mz2 = params(imz)
DO k = nlev, 1, -1   

   ! Retrieve current (local) state variable values.
   tC     = Temp(k)
   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif
   Varout(oPAR_,k)=par_
   NO3    = Vars(iNO3,    k)
   MIC    = Vars(iZOO,    k)   !Microzooplankton
   MES    = Vars(iZOO2,   k)   !Mesozoo.
   DET    = Vars(iDET,    k)
   DETFe  = Vars(iDETFe,  k)
   PHY    = Vars(iPHY(1), k)
   Varout(oPON,k)=MIC+MES+DET+PHY
   PMUPHY = Vars(iPMU,    k)
   VARPHY = Vars(iVAR,    k)
   MIoPHY = Vars(iMIo,    k)
   VIoPHY = Vars(iVIo,    k)
   MToPHY = Vars(iMTo,    k)
   VToPHY = Vars(iVTo,    k)

   PMU    = PMUPHY/PHY
   VAR    = VARPHY/PHY-PMU**2  !Correct VAR to the real one
   VAR    = max(VAR, eps)
   PMU    = PMU - PMU0         !Correct PMU to the real one
   MIo    = MIoPHY/PHY
   VIo    = VIoPHY/PHY-MIO**2
   VIo    = max(VIo, eps)
   MTo    = MToPHY/PHY
   VTo    = VToPHY/PHY-MTo**2
   VTo    = max(VTo, eps)
   Fe     = Vars(ifer,k)

   !x: ln optimal light (umol photons m-2 s-1)
   !Z: optimal temperature (ºC)
   !L: log cell volume (um^3)

   call citrate_PHY3(NO3,PAR_,tC, Fe, PMU, MIo, MTo, muNet,   &
     dmuNetdl,d2muNetdl2, dMudI, d2MudI2, dMudT, d2MudT2,     &
     theta, QN, &
     dthetadl, d2thetadl2, dQNdL, d2QNdL2,  &
     dthetadT, d2thetadT2,&
     dthetadI, d2thetadI2,SI, Lno3, .TRUE.)

   Varout(oTheta(1),k)= theta+0.5*(VAR*d2thetadl2+VIo*d2thetadI2&
                      + VTo*d2thetadT2)
   Varout(oSI(1)   ,k)= SI
   Varout(oLno3(1) ,k)= Lno3
   Varout(oQN(1)   ,k)= QN + 0.5*VAR*d2QNdL2
   Varout(odMudI   ,k)= dMudI
   Varout(odMudT   ,k)= dMudT
   Varout(od2mudI2 ,k)= d2MudI2
   Varout(od2mudT2 ,k)= d2MudT2
!=============================================================
! ZOOplankton section:
   tf_z = TEMPBOL(Ez,tC)
    
   IF (mod(it, nsave) .EQ. 1) THEN  !Calculate the below only when necessary
     ! Calculate the integration  
     L_min  = PMU-3D0*sqrt(VAR)  ! Minimal size
     L_max  = PMU+3D0*sqrt(VAR)  ! MadLimal size
     call GEN_SEQ(L_min,L_max,M,L_,dL)
     T_min  = MTo-3.*sqrt(VTo)
     T_max  = MTo+3.*sqrt(VTo)
     call GEN_SEQ(T_min,T_max,M,T_,dT)
     I_min  = MIo-3.*sqrt(VIo)
     I_max  = MIo+3.*sqrt(VIo)
     call GEN_SEQ(I_min,I_max,M,I_,dI)

     CHLt   =0D0
     pCHL(:)=0D0  ! Size fractionated Chl

     do i=1,M
       do j=1,M
         do k1=1,M
       !Calculate Chl in this size class
           call citrate_PHY3(NO3,PAR_,tC,Fe,L_(i),I_(j),T_(k1),mu2,&
                   dmudL2,d2mudL22, dMudx2, d2Mudx22, dMudZ2,&
                   d2MudZ2, theta2, QN2, dthetadl2, d2thetadl22,&
                   dQNdL2, d2QNdL22, dthetadZ2,&
                   d2thetadZ22, dthetadx2, d2thetadx22, &
                   SI, Lno3,.FALSE.)

           !Probability function
           Pl = PHY*normal(PMU,VAR,L_(i))*normal(MIo,VIo,I_(j)) &
                   *normal(MTo,VTo,T_(k1))

           ! Chl in this size class:
           dCHL = Pl*dL*dI*dT*theta2/QN2

           ! Calculate total Chl:
           CHLt = CHLt + dCHL

       !   Calculate size-fractionated Chl:
       !   Here the PMU_1 is already the normal PMU + log(10) (to avoid negative values)
           if (L_(i) .le. (PMU_1-PMU0)) then
              pCHL(4)=pCHL(4)+dCHL
           elseif (L_(i) .le. (PMU_3-PMU0)) then
              pCHL(3)=pCHL(3)+dCHL
           elseif (L_(i) .le. (PMU_10-PMU0)) then
              pCHL(2)=pCHL(2)+dCHL
           else
              pCHL(1)=pCHL(1)+dCHL
           endif
         enddo
       enddo
     enddo

     do j = 1, 4
        Varout(oCHLs(j),k)=pCHL(j)/max(CHLt,eps) !Get the percentage of size-fractionated CHL
     enddo

     ! Total NPP (use light as the in situ incubation):
     call citrate_PHY3(NO3,PAR(k),tC, Fe, PMU, MIo, MTo, muNet1,   &
       dmuNetdl1,d2muNetdl21, dMudI1, d2MudI21, dMudT1, d2MudT21,  &
       theta1, QN1, &
       dthetadl1, d2thetadl21, dQNdL1, d2QNdL21,  &
       dthetadT1, d2thetadT21,&
       dthetadI1, d2thetadI21,SI1,Lno31, .TRUE.)

     d2muNet_QNdl2 = d2Y_Xdl2(muNet1,QN1,dmuNetdl1, dQNdL1,&
                              d2muNetdl21, d2QNdl21) 
     d2muNet_QNdx2 = d2mudI21/QN1**2
     d2muNet_QNdZ2 = d2mudT21/QN1**2

     ! Primary production (ug C L-1 d-1)
     Varout(oPPt,k)= 12.*PHY*(muNet1/QN1 + 0.5*VAR*d2muNet_QNdl2 &
                                 + 0.5*VIo*d2muNet_QNdx2 &
                                 + 0.5*VTo*d2muNet_QNdZ2) 
     Varout(oPPt,k)= max(Varout(oPPt,k),0d0)
     Varout(oCHLt,  k)= CHLt
     Varout(oCHL(1),k)= CHLt
   Endif

! Microzoo.
! Feeding selectivity at the mean size:
   Cf   = gb1
   rl   = exp(Cf*PMU)

! Calculate total platable prey:
   Ptot = PHY*rl*(1d0+VAR/2d0*Cf**2)

! The grazing dependence on total prey
   gbar = grazing(grazing_formulation,Kp1,Ptot)

   !MicroZooplankton per capita total ingestion rate
   INGES1 = tf_z*gmax1*gbar

   !Zooplankton excretion rate (-> DIN)
   RES1 = INGES1*(1D0-GGE-unass1)

   !ZOOPLANKTON EGESTION (-> DETRITUS)
   EGES1 = INGES1 * unass1

   ! Microzoo Grazing rate on the mean size
   gbar = INGES1*MIC*rl/Ptot

! First derivative of microzoo. grazing rate evaluated at the mean size (negative)

   dgdlbar1  = Cf*gbar*2d0*Kp1**2/(Kp1**2+Ptot**2) 
   d2gdl2bar1= 4d0*Cf**2*Kp1**2*gbar*(Kp1**2-Ptot**2)/((Kp1**2+Ptot**2)**2)

! Mesozoo.
! Feeding selectivity at the mean size:
   Cf   = gb2
   rl   = exp(Cf*(PMU-PMU10))

! Calculate total platable prey:
   Ptot = PHY*rl*(1d0+VAR/2d0*Cf**2)

! The grazing dependence on total prey
   gbar = grazing(grazing_formulation,Kp2,(Ptot+MIC))

   !MesoZooplankton per capita total ingestion rate
   INGES2 = tf_z*gmax2*gbar

   ! Mesozoo Grazing rate on the mean size
   gbar = INGES2*MES*rl/(Ptot+MIC)

! First derivative of mesozoo. grazing rate evaluated at the mean size (negative)

   ! An intermediate variable to facilitate computation:
   cff = (Kp2**2-(Ptot+MIC)**2)/(Kp2**2+(Ptot+MIC)**2)
   cff1= cff*Ptot/(Ptot+MIC)+1d0

   dgdlbar2  = Cf*gbar*cff1
   d2gdl2bar2= Cf*gbar*( cff1**2                             &
       + Cf*Ptot*(-4d0*Kp2**2*Ptot/((Ptot+MIC)**2+Kp2**2)**2 &
       + cff*MIC/((Ptot+MIC)**2)) )

   !Zooplankton excretion rate (-> DIN)
   RES2 = INGES2*(1D0-GGE-unass2)

   !ZOOPLANKTON EGESTION (-> DETRITUS)
   EGES2 = INGES2*unass2

! End of zooplankton section
!=============================================================
!! Solve ODE functions:
   Zmort2 = MES*MES*mz2*tf_z  !Mortality term for MES

   VAR1  = VAR+dtdays   &
   * (VAR*(VAR*(d2muNetdl2-d2gdl2bar1-d2gdl2bar2))+2d0*VTR_L*muNet)

   VIo1  = VIo+dtdays*(VIo*VIo*d2mudI2+2d0*VTR_I*muNet)
   VTo1  = VTo+dtdays*(VTo*VTo*d2mudT2+2d0*VTR_T*muNet)
 
   !Eq. 18, Update PMU:
   PMU = PMU + PMU0  ! Restore to positive values
   PMU1= PMU + dtdays*(VAR*(dmuNetdl-dgdlbar1-dgdlbar2))
   MTo1= MTo + dtdays*VTo*dmudT
   MIo1= MIo + dtdays*VIo*dMudI

!All rates have been multiplied by dtdays to get the real rate correponding to the actual time step
   PP_ND= dtdays*RDN*DET*tf_z   
   PP_NZ= (MIC*RES1+MES*RES2)*dtdays   !Sum of MIC and MES to DIN   
   PP_DZ= (MIC*EGES1+MES*EGES2+Zmort2)*dtdays  !Sum of MIC and MES to detritus
   PP_MIC_P  =MIC*INGES1*dtdays      
   PP_MES_P  =MES*INGES2*dtdays*Ptot/(Ptot+MIC)
   PP_MES_MIC=MES*INGES2*dtdays* MIC/(Ptot+MIC)
   PP_PN= PHY*dtdays*(muNet+0.5*(VAR*d2muNetdl2 + VIo*d2MudI2 &
                                                + VTo*d2MudT2))
   PP_PN= max(PP_PN,0d0)

!Update tracers:
   DET1 = DET + PP_DZ - PP_ND 
   NO3  = NO3 + PP_ND + PP_NZ - PP_PN
   PHY1 = PHY + PP_PN - PP_MIC_P - PP_MES_P
   MIC  = MIC + PP_MIC_P*GGE - PP_MES_MIC
   MES  = MES + (MES*INGES2*GGE - Zmort2)*dtdays
   
   call IRONCYCLE(tC, DET, PP_NZ, PP_ND, PP_PN, PP_DZ, Fescav, DETFe, Fe)

   Varout(oFescav,k) = Fescav

   PMUPHY = PMUPHY + PHY*(PMU1-PMU) + PMU*(PHY1-PHY)
   VARPHY = VARPHY + PHY*(VAR1-VAR+2.*PMU*(PMU1-PMU)) &
          + (VAR+PMU**2)*(PHY1-PHY)
   MIoPHY = MIoPHY + PHY*(MIo1-MIo) + MIo*(PHY1-PHY)
   VIoPHY = VIoPHY + PHY*(VIo1-VIo+2.*MIo*(MIo1-MIo)) &
          + (VIo+MIo**2)*(PHY1-PHY)

   MToPHY = MToPHY + PHY*(MTo1-MTo) + MTo*(PHY1-PHY)
   VToPHY = VToPHY + PHY*(VTo1-VTo+2.*MTo*(MTo1-MTo)) &
          + (VTo+MTo**2)*(PHY1-PHY)

   Varout(oNO3,k)      = NO3
   Varout(oPHY(1),k)   = PHY1
   Varout(oPHYt,  k)   = PHY1
   Varout(oZOO,k)      = MIC
   Varout(oZOO2,k)     = MES
   Varout(oDET,k)      = DET1
   Varout(oDETFe,k)    = DETFe
   Varout(ofer,k)      = Fe
   Varout(oPMU,k)      = PMUPHY
   Varout(oMTo,k)      = MToPHY
   Varout(oMIo,k)      = MIoPHY
   Varout(oVAR,k)      = VARPHY
   Varout(oVTo,k)      = VToPHY
   Varout(oVIo,k)      = VIoPHY
   Varout(omuNet(1),k) = muNet               !Growth rate of mean size
   Varout(omuAvg,   k) = PP_PN/dtdays/PHY    !Avg. growth rate of the total community
   Varout(oGraz(1) ,k) = PP_MIC_P/PHY/dtdays        !Microzoo. grazing rate on PHY
   Varout(oMESg    ,k) = PP_MES_P/PHY/dtdays        !Mesozoo. grazing rate on PHY
   Varout(oMESgMIC ,k) = PP_MES_MIC/MIC/dtdays      !Mesozoo.  grazing on MIC
   Varout(odgdl1   ,k) = dgdlbar1
   Varout(odgdl2   ,k) = dgdlbar2
   Varout(od2gdl1  ,k) = d2gdl2bar1
   Varout(od2gdl2  ,k) = d2gdl2bar2
   Varout(oZ2N     ,k) = PP_NZ/dtdays
   Varout(oD2N     ,k) = PP_ND/dtdays
   Varout(odmudl   ,k) = dmuNetdl
   Varout(od2mu    ,k) = d2muNetdl2
ENDDO
END SUBROUTINE Citrate3_MOD

! The subroutine for 3-trait phytoplankton with derivatives calculated
subroutine citrate_PHY3(NO3,PAR_,Temp_, Fe, L, x, Z, mu,   &
     dmudL,d2mudL2, dMudx, d2Mudx2, dMudZ, d2MudZ2, theta, QN, &
     dthetadl, d2thetadl2, dQNdL, d2QNdL2, dthetadZ, d2thetadZ2,&
     dthetadx, d2thetadx2, SI, fN,cal_deriv)

!x: ln optimal light (umol photons m-2 s-1)
!Z: optimal temperature (ºC)
!L: log cell volume (um^3)

use BIO_MOD, only: thetamax,thetamin
use BIO_MOD, only: ScaleTrait,dY_Xdl,d2Y_Xdl2
use BIO_MOD, only: K0Fe,alphaFe,imu0,params

implicit none
real, intent(in)  :: L, x, Z, NO3, PAR_,Temp_, Fe 
logical,intent(in):: cal_deriv
real, intent(out) :: mu, dmudL, d2mudL2
real, intent(out) :: dMudx, d2Mudx2, dMudZ, d2MudZ2 
real, intent(out) :: theta, QN, dthetadl, d2thetadl2
real, intent(out) :: dQNdL, d2QNdL2,dthetadZ,d2thetadZ2
real, intent(out) :: dthetadx, d2thetadx2
real, intent(out) :: fN, SI
real :: umax, alpha,H, G, Iopt, hat, dhatdl, d2hatdl2,cff 
real :: Kn,dfNdl,d2fNdl2
real :: KFe,fFe
real :: dGdL, d2GdL2, Y, dYdx, d2Ydx2, dHdx, d2Hdx2
real :: P, ht, f, dPdZ, d2PdZ2, dfdZ, d2fdZ2
real :: cff1,Qmin,Qmax
real :: Z1, Z2
real ::  mu0b= .3
real, parameter   ::  b   = .21
real, parameter   ::  alphaK = .27, K0N = 0.2
real, parameter   ::  alphamu= .2,  betamu = -0.01
real, parameter   ::  Eb  = .0633
real, parameter   ::  a0  = .34, k = -.47 
real, parameter   ::  w   = 14., T0= 15d0
real, parameter   ::  Q0min = 0.15, Q0max= 0.3
real, parameter   ::  alphaQmin = -0.16, alphaQmax= -0.07

!Initialize outputs:
mu   =0.; dmudL=0.; d2mudL2=0.;
dMudx=0.; d2Mudx2=0.; dMudZ=0.; d2MudZ2=0.; 
theta=0.; QN=0.; dthetadl=0.; d2thetadl2=0.;
dQNdL=0.; d2QNdL2=0.;dthetadZ=0.;d2thetadZ2=0.;
dthetadx=0.; d2thetadx2=0.;

!Light term
if (PAR_ <= 0.) stop("PAR <= 0!")
mu0b = params(imu0)
umax = mu0b*exp(b*x)
alpha=   a0*exp(k*x)
Iopt =      exp(x)

H = 1./(PAR_/(alpha*Iopt*Iopt) + 1./umax - 2./(alpha*Iopt) + &
          1./(PAR_*alpha))
SI= H/umax
!Note here suboptimal is optimal (even if I < Iopt, because umax is higher, mu is also higher)

!Size (nutrient) term (L: lnvolume)
hat      = exp(alphamu*L + betamu*L**2)

if (cal_deriv) then
  dhatdl   = hat*(alphamu + 2.*betamu * L)
  d2hatdl2 = 2.*hat*betamu+hat*(alphamu+2.*betamu*L)**2
endif

! Nitrogen limitation:
Kn  = ScaleTrait(L, K0N, alphaK)
fN  = NO3/(NO3 + Kn)  !Nitrogen limitation index

! Iron limitation:
KFe = ScaleTrait(L, K0Fe, alphaFe)
fFe = Fe/(Fe + KFe)

if (fN .ge. fFe) then !Iron limitation
   call MM_2der(Fe, K0Fe,alphaFe, L, Kn,fN,dfNdl,d2fNdl2)
else                  !Nitrogen limitation
   call MM_2der(NO3,K0N,  alphaK, L, Kn,fN,dfNdl,d2fNdl2)
endif

G = hat * fN
if (cal_deriv) then
  !Calculate the derivatives of G
  dGdL   = hat*dfNdl + dhatdl*fN
  d2GdL2 = hat*d2fNdl2 + 2.*dfNdl*dhatdl + fN*d2hatdl2
  
  !Calculate the derivatives of H
  Y      = 1./H
  cff    = PAR_/(alpha*Iopt*Iopt)
  dYdx   = -(k+2.)*cff &
     - b/umax +2.*(k+1.)/(alpha*Iopt) - k/(PAR_*alpha)
  
  d2Ydx2 = (k+2.)**2*cff + b*b/umax  &
     - 2.*(k+1.)**2/(alpha*Iopt) + k**2/(PAR_*alpha)
  
  dHdx   = -H*H*dYdx
  d2Hdx2 = (2.*H*H*H)*dYdx*dYdx - d2Ydx2*H*H
endif
!Temperature term
ht      = exp(Eb*(Temp_ - T0))
P       = 1.-((Temp_ - Z)/w)**2
f       = ht * P

if (cal_deriv) then
   dPdZ    = 2.*(Temp_  - Z)/(w*w)
   d2PdZ2  = -2./(w*w)
   dfdZ    = ht*dPdZ 
   d2fdZ2  = ht*d2PdZ2 
endif

f       = max(f, 0d0)
Mu      = H*G*f
if (cal_deriv) then
  dMudx   = G*f*dHdx
  d2Mudx2 = G*f*d2Hdx2
  dMudL   = f*H*dGdL
  d2MudL2 = f*H*d2GdL2
  dMudZ   = H*G*dfdZ
  d2MudZ2 = H*G*d2fdZ2
endif
!N:C ratio at avg. size (assumed independent of light and temp.)

!Make Qmin and Qmax size dependent
Qmin=Q0min*exp(alphaQmin*L)
Qmax=Q0max*exp(alphaQmax*L)

cff1= 1.-Qmin/Qmax
QN  = Qmin/(1.-cff1*fN)

if (cal_deriv) then
  Z1  = alphaQmin - QN*(fN*(alphaQmin-alphaQmax)/Qmax +  &
                                (1./Qmax - 1./Qmin)*dfNdl)
  dQNdL = QN*Z1
  
  cff = (dfNdl-alphaQmax*fN)/Qmax*(alphaQmin-alphaQmax) + &
      (1./Qmax - 1./Qmin)*d2fNdl2 + &
      dfNdl*(alphaQmin/Qmin - alphaQmax/Qmax)
  
  Z2 = (alphaQmin-alphaQmax)*fN/Qmax + (1./Qmax-1./Qmin)*dfNdl
  Z2 = Z2*dQNdL
  d2QNdL2 = dQNdL*Z1 - QN*(Z2 + QN*cff)
endif

!Chl:C ratio
cff     = (thetamax - thetamin)/PAR_

!Chl:C ratio at avg. traits
theta   = thetamin+Mu/alpha*cff   !Unit: gChl/molC
  
if (cal_deriv) then
  !alpha only depends on Iopt
  !Derivatives at mean size
  dthetadl   = cff/alpha*dMudL
  d2thetadl2 = cff/alpha*d2MudL2
  
  !Derivatives at mean temp.
  dthetadZ   = cff/alpha*dMudZ
  d2thetadZ2 = cff/alpha*d2MudZ2
  
  !Derivatives at mean Iopt.
  dthetadx   = cff*dY_Xdl(mu, alpha, dMudx, k*alpha) 
  d2thetadx2 = cff*d2Y_Xdl2(mu,alpha,dMudx, k*alpha, d2Mudx2, &
                                               k**2*alpha)
endif
return
end subroutine 

!Mechaelis-Mention functions and derivatives
subroutine MM_2der(N,K0,alphaK,L,Kn,fN, dfNdl, d2fNdl2)
use BIO_MOD, only: ScaleTrait
implicit none
real,   intent(in)  :: N, K0, alphaK, L
real,   intent(out) :: fN, dfNdl, d2fNdl2
real,   intent(out) :: Kn
real :: cff
   ! Half saturation constant for growth at avg. size
   Kn = ScaleTrait(L, K0, alphaK)
   cff= N+Kn
   fN = N/(N + Kn) ! Nutrient limitation index at avg. size
   dfNdl = -alphaK*Kn*N/(N+Kn)**2
 d2fNdl2 = -alphaK**2*N*Kn*(1./(N+Kn)**2 - 2.*Kn/(cff*cff*cff))
 return
end subroutine MM_2der

subroutine GEN_SEQ(xmin, xmax, length, x, dx)
    !Generate a numeric vector with min, max and length
    implicit none
    real,   intent(in)  :: xmin, xmax
    integer,intent(in)  :: length
    real,   intent(out) :: x(length),dx
    integer :: i

    if (xmin >= xmax) then
      write(*,'(F12.3,A20,F12.3)') xmin," not smaller than ",xmax
      stop
    endif
    dx    = (xmax - xmin)/float(length-1)
    x(1)  = xmin
    do i  = 2, length
      x(i)= x(i-1) + dx  !Log size of each class
    enddo
end subroutine GEN_SEQ
