SUBROUTINE NPZD_CONT
use bio_MOD
use MOD_1D, only: it, nsave
implicit none

integer :: k,j,i
!INPUT PARAMETERS:
real :: tC,par_
!LOCAL VARIABLES of phytoplankton:
real, parameter :: thetamin = 0.02
real :: PMUPHY,VARPHY,PHY1
real :: NO3, Fe, PHY, MIC, MES, DET, DETFe, DET1
real :: PMU,VAR,PMU1,VAR1
real :: mu0,aI0
real :: QN,Qnmax,Qnmin  ! cell quota related variables
real :: muNet,dmuNetdl,d2muNetdl2, d3mudl3, d4mudl4
real :: dmuNetdl1,d2muNetdl21, d3mudl31, d4mudl41
real :: alphaI,SI,Lno3,KFe_,Fescav
real :: muNet1, SI1, Lno31, QN1
real :: rmax_T ! a scratch variable to temporally store phyto. growth rate
real :: theta,theta1,dthetadl,d2thetadl2,dQNdL,d2QNdL2
real :: dthetadl1,d2thetadl21,dQNdL1,d2QNdL21
real :: VTR,d2muNet_QNdl2

!Declarations of zooplankton:
real :: dgdlbar1,d2gdl2bar1,dgdlbar2,d2gdl2bar2
real :: INGES1,INGES2,RES1,RES2,EGES1,EGES2,Zmort2,mz2
real :: dx  ! size interval
real :: Kp1,Cf, cff, cff1  !Zooplankton variables
real :: pp_PN, PPpn, PP_MIC_P, PP_MES_P,PP_MES_MIC
real :: Ptot,  CHLt,KN, Pl, dCHL,NPPt,rl
real :: pCHL(4) = 0D0

integer, parameter :: M     = 60    !discretize the continous normal distribution
real               :: x(M)  = 0d0

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

!-----------------------------------------------------------------------
VTR    = params(iVTR)
alphaI = params(ialphaI)
Kp1    = params(iKPHY)              ! Grazing half-saturation constant for MIC
mz2    = params(imz)
thetamax=.47
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

   PMU    = PMUPHY/PHY
   VAR    = VARPHY/PHY-PMU**2  !Correct VAR to the real one
   PMU    = PMU - PMU0         !Correct PMU to the real one
   Fe     = Vars(ifer,k)

   call PHY_NPZDCONT(NO3,par_,tC,Fe,PMU,muNet,dmuNetdl,d2muNetdl2,d3mudl3,d4mudl4,SI,Lno3, theta, &
    QN,dQNdL,d2QNdL2,dthetadL,d2thetadl2)

   Varout(oTheta(1),k)= theta+0.5*VAR*d2thetadl2
   Varout(oQN(1)   ,k)= QN   +0.5*VAR*d2QNdl2

   Varout(oSI(1),k)  =SI
   Varout(oLno3(1),k)=Lno3

!=============================================================
! Calculate the community sinking rate of phytoplankton
! Phytoplankton sinking rate at the average size
!    w_p    = ScaleTrait(PMU,abs(w_p0),alphaW)  !Positive
!    dwdl   = alphaW*w_p
!    d2wdl2 = alphaW*dwdl
! 
!  ! Sinking rate (m) of phytoplankton
!    w_pAvg = w_p+0.5*VAR*d2wdl2
!  ! sinking rate must be negative!
!    w_pAvg = -min(w_pAvg, Wmax)
!=============================================================
!! ZOOplankton section:
   tf_z = TEMPBOL(Ez,tC)
    
   IF (mod(it, nsave) .EQ. 1) THEN  !Calculate the below only when necessary
     ! Calculate the integration  
     x(1)  = PMU-6D0*sqrt(VAR)  ! Minimal size
     x(M)  = PMU+6D0*sqrt(VAR)  ! Maximal size
     dx    = (x(M) - x(1))/float(M-1)
     do i  = 2,M-1
       x(i) = x(i-1) + dx  !Log size of each class
     enddo

     ! Calculate total phytoplankton palatable prey (Ptot,Eq. 13):
     ! Ptot roughly equals to Phy
     CHLt  = 0D0
     pCHL(:)=0D0  ! Size fractionated Chl

     do i=1,M
       !Calculate Chl in this size class

       !Pertain only meaningful size class
       if (x(i) .gt. -3. .and. x(i) .le. 15.) then 
         mu0=params(imu0)  *exp(alphamu*x(i)+betamu*x(i)**2)
         aI0=params(iaI0_C)*exp(alphaI*x(i))
         KN =params(iKN)   *exp(alphaK*x(i))
        KFe_=params(iKFe)  *exp(alphaFe*x(i))

        ! Minimal N:C ratio
        Qnmin=0.06
        ! Maximal N:C ratio
        Qnmax=3.*Qnmin

        tf_p=TEMPBOL(Ep,tC)
        ! The maximal growth rate (rmax_T) under temperature tC 
        rmax_T = mu0*tf_p

        !The light limitation index (SI)
        ! Unit of aI0_C: (W m-2)-1, different from other models (aI0_Chl)
        ! Include photoinhibition (Platt et al. J Phycol 1980)
        SI = 1. - exp(-aI0*par_/params(imu0)/tf_p)

        ! Growth rate (d-1) is the product function of temperature, light, and nutrient   
        Lno3  = NO3/(NO3 + KN)

        if (DO_IRON) then
           Lno3  = min(Lno3, Fe/(Fe + KFe_))
        endif
        muNet1 = rmax_T*Lno3*SI

        QN = Qnmin/(1d0-(1d0-Qnmin/Qnmax)*Lno3)

        theta  = thetamin+muNet1/PAR_/aI0*(thetamax-thetamin)   !Unit: gChl/molC
        !Probability function
        Pl   = PHY*normal(PMU,VAR,x(i))

        ! Chl in this size class:
        dCHL = Pl*dx*theta/QN

         ! Calculate total Chl:
        CHLt = CHLt + dCHL

     !   Calculate size-fractionated Chl:
     !   Here the PMU_1 is already the normal PMU + log(10) (to avoid negative values)
         if (x(i) .le. (PMU_1-PMU0)) then
            pCHL(4)=pCHL(4)+dCHL
         elseif (x(i) .le. (PMU_3-PMU0)) then
            pCHL(3)=pCHL(3)+dCHL
         elseif (x(i) .le. (PMU_10-PMU0)) then
            pCHL(2)=pCHL(2)+dCHL
         else
            pCHL(1)=pCHL(1)+dCHL
         endif
       endif
     enddo

     do j = 1, 4
        Varout(oCHLs(j),k)=pCHL(j)/max(CHLt,eps) !Get the percentage of size-fractionated CHL
     enddo

     ! Total NPP (use light as the in situ incubation):
     call PHY_NPZDCONT(NO3,PAR(k),tC,Fe,PMU,muNet1,dmuNetdl1,d2muNetdl21,&
         d3mudl31,d4mudl41,SI1,Lno31, theta1, &
      QN1,dQNdL1,d2QNdL21,dthetadL1,d2thetadl21)

     d2muNet_QNdl2 = d2Y_Xdl2(muNet1,QN1,dmuNetdl1, dQNdL1, d2muNetdl21, d2QNdl21) 
     NPPt          = PHY*(muNet1/QN1+0.5*VAR*d2muNet_QNdl2) 

     Varout(oCHLt,k) = CHLt
     Varout(oPPt, k) = NPPt*12d0 ! Correct the unit to ug C/L
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
   * (VAR*(VAR*(d2muNetdl2-d2gdl2bar1-d2gdl2bar2 + VTR*d4mudl4) &
   - 5d0*VTR*d2muNetdl2)+2d0*VTR*muNet)
  
   !Eq. 18, Update PMU:
   PMU = PMU + PMU0  ! Restore to positive values
   PMU1= PMU + dtdays*(VAR*(dmuNetdl-dgdlbar1-dgdlbar2+VTR*d3mudl3) - 3d0*VTR*dmuNetdl)

!All rates have been multiplied by dtdays to get the real rate correponding to the actual time step
   PP_ND= dtdays*RDN*DET*tf_z   
   PP_NZ= (MIC*RES1+MES*RES2)*dtdays     !Sum of MIC and MES to DIN   
   PP_DZ= (MIC*EGES1+MES*EGES2+Zmort2)*dtdays  !Sum of MIC and MES to detritus
   PP_MIC_P = MIC*INGES1*dtdays      
   PP_MES_P = MES*INGES2*dtdays*Ptot/(Ptot+MIC)
   PP_MES_MIC=MES*INGES2*dtdays* MIC/(Ptot+MIC)
   PP_PN= PHY*dtdays*(muNet+0.5*VAR*(d2muNetdl2+VTR*d4mudl4)-1.5*VTR*d2muNetdl2)

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

   Varout(odVAR,k)     = (VAR1-VAR)/dtdays
   Varout(oNO3,k)      = NO3
   Varout(oPHY(1),k)   = PHY1
   Varout(oPHYt,  k)   = PHY1
   Varout(oZOO,k)      = MIC
   Varout(oZOO2,k)     = MES
   Varout(oDET,k)      = DET1
   Varout(oDETFe,k)    = DETFe
   Varout(ofer,k)      = Fe
   Varout(oPMU,k)      = PMUPHY
   Varout(oVAR,k)      = VARPHY
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
   Varout(od3mu    ,k) = d3mudl3 
   Varout(od4mu    ,k) = d4mudl4
ENDDO
END SUBROUTINE NPZD_CONT 

! The subroutine only for phytoplankton
subroutine PHY_NPZDCONT(NO3,PAR_,Temp_,Fe, PMU, muNet,dmudl,d2mudl2,d3mudl3,d4mudl4,SI,fN, theta, &
    QN,dQNdL,d2QNdL2,dthetadL,d2thetadl2)
use bio_MOD, only : ScaleTrait, TEMPBOL, params, dY_Xdl, d2Y_Xdl2
use bio_MOD, only : ialphaI,  imu0, iaI0_C, iKN
use bio_MOD, only : Ep, K0Fe, alphaFe, KFe, do_IRON,betamu, alphamu, iKFe
implicit none
real, intent(in)  :: PMU, NO3, PAR_,Temp_, Fe 
real, intent(out) :: muNet, dmudl, d2mudl2, theta, QN, SI, fN
real, intent(out) :: d3mudl3, d4mudl4, dQNdL,d2QNdL2,dthetadL,d2thetadl2
real :: alphaK = 0.27, alphaI, mu0hat
real :: dmu0hatdl, d2mu0hatdl2, aI, cff,cff1,daI_mu0hatdl
real :: daI_mu0hat2dl
real :: d3mu0hatdl3, d4mu0hatdl4
real :: dSIdl,dcffdl
real :: mu0hatSI, mu0hat1 
real :: dmu0hatSIdl, dmu0hat_aIdl,d2mu0hat_aIdl2
real :: d2aI_mu0hatdl2,d3aI_mu0hatdl3,d2SIdl2,d2mu0hatSIdl2
real :: d3SIdl3,daI_mu0hat3dl,d2aI_mu0hat2dl2,d4SIdl4
real :: daI_mu0hat4dl,d2aI_mu0hat3dl2,d3aI_mu0hat2dl3, d4aI_mu0hatdl4
real :: d3muIhatdl3, d4muIhatdl4
real :: Kn,K0N,dfNdl,d2fNdl2,d3fNdl3, d4fNdl4, Qmin,Qmax,tf
real :: fFe
real, parameter :: thetamin = 0.02, thetamax = 0.47

alphaI     =params(ialphaI)
tf         =TEMPBOL(Ep,Temp_)
cff1       =alphamu + 2.* betamu * PMU
mu0hat     =tf*params(imu0)*exp(alphamu*PMU + betamu*PMU**2)
dmu0hatdl  =mu0hat*cff1
d2mu0hatdl2=mu0hat*2.*betamu+mu0hat*cff1**2
d3mu0hatdl3=(2.*betamu+cff1**2)*dmu0hatdl+4.*betamu*mu0hat*cff1
d4mu0hatdl4=dmu0hatdl*8.*betamu*cff1+(2.*betamu+cff1**2)*d2mu0hatdl2+8.*betamu**2*mu0hat

! Initial slope of P-I curve
aI=ScaleTrait(PMU, params(iaI0_C), alphaI)

!The light limitation index (SI)
mu0hat1 =tf*params(imu0)  !To simplify the size dependence of aI
!cff=exp(-aI*PAR_/mu0hat)

!mu0hat1 = mu0hat when alphamu and betamu = 0
cff=exp(-aI*PAR_/mu0hat1)
SI =1.-cff

!daI_mu0hatdl=aI*(alphaI-alphamu-2D0*betamu*PMU)/mu0hat
daI_mu0hatdl=aI*alphaI/mu0hat1

dSIdl=cff*PAR_*daI_mu0hatdl

mu0hatSI    = mu0hat*SI
dmu0hatSIdl = mu0hatSI*cff1 + mu0hat*dSIdl
!dmu0hat_aIdl= (dmu0hatdl - alphaI*mu0hat)/aI
!d2mu0hat_aIdl2 = mu0hat/aI*2D0*betamu+(alphamu-alphaI+2D0*betamu*PMU)*dmu0hat_aIdl

!daI_mu0hat2dl  = daI_mu0hatdl/mu0hat - aI/mu0hat**3*dmu0hatdl !Correct

!d2aI_mu0hatdl2 = alphaI*daI_mu0hatdl - aI/mu0hat**2*d2mu0hatdl2 - daI_mu0hat2dl*dmu0hatdl

!d3aI_mu0hatdl3 = d2aI_mu0hatdl2*(alphaI-alphamu-2d0*betamu*PMU)-4d0*betamu*daI_mu0hatdl

!d2SIdl2 = -PAR_*dSIdl*daI_mu0hatdl+PAR_*cff*d2aI_mu0hatdl2
!d2SIdl2=PAR_*alphaI**2*aI*cff/mu0hat1*(1d0-PAR_*aI/mu0hat1)
d2SIdl2=PAR_*alphaI/mu0hat1*aI*(cff*alphaI-dSIdl)

d2mu0hatSIdl2 = cff1*dmu0hatSIdl + 2.*betamu*mu0hatSI + mu0hat*cff1*dSIdl + mu0hat*d2SIdl2  !Correct

!d3SIdl3 = PAR_*(-2.*dSIdl*d2aI_mu0hatdl2 - d2SIdl2*daI_mu0hatdl + (1.-SI)*d3aI_mu0hatdl3)  !Correct
d3SIdl3=PAR_/mu0hat1*alphaI**3*aI*cff*((1.-PAR_/mu0hat1*aI)**2-aI*PAR_/mu0hat1)

!daI_mu0hat3dl  = daI_mu0hat2dl/mu0hat - aI*dmu0hatdl/mu0hat**4 !Correct
!
!d2aI_mu0hat2dl2= d2aI_mu0hatdl2/mu0hat - daI_mu0hatdl/mu0hat**2*dmu0hatdl - daI_mu0hat3dl * dmu0hatdl - aI/mu0hat**3*d2mu0hatdl2 !Correct 
!
!d3aI_mu0hatdl3 = alphaI*d2aI_mu0hatdl2 - aI/mu0hat**2*d3mu0hatdl3 -2d0*daI_mu0hat2dl*d2mu0hatdl2 -d2aI_mu0hat2dl2*dmu0hatdl !Correct
!
!daI_mu0hat4dl  = daI_mu0hat3dl/mu0hat - aI/mu0hat**5*dmu0hatdl  !Correct
!
!d2aI_mu0hat3dl2= d2aI_mu0hat2dl2/mu0hat - daI_mu0hat2dl/mu0hat**2*dmu0hatdl - aI/mu0hat**4*d2mu0hatdl2 - daI_mu0hat4dl*dmu0hatdl   !Correct
!
!d3aI_mu0hat2dl3= (d3aI_mu0hatdl3*mu0hat - daI_mu0hatdl*d2mu0hatdl2)/mu0hat**2 -2d0/mu0hat**3*(d2aI_mu0hatdl2*mu0hat - daI_mu0hatdl*dmu0hatdl)*dmu0hatdl - (aI/mu0hat**3*d3mu0hatdl3+2d0*daI_mu0hat3dl*d2mu0hatdl2+d2aI_mu0hat3dl2*dmu0hatdl)  !Correct
!
!d4aI_mu0hatdl4 = alphaI*d3aI_mu0hatdl3 - (aI/mu0hat**2*d4mu0hatdl4 + 3.* daI_mu0hat2dl*d3mu0hatdl3+3.*d2aI_mu0hat2dl2 * d2mu0hatdl2)- d3aI_mu0hat2dl3*dmu0hatdl  !Correct

!d4SIdl4 = PAR_*((1.-SI)*d4aI_mu0hatdl4 - 3.*dSIdl*d3aI_mu0hatdl3 - 3.*d2aI_mu0hatdl2*d2SIdl2-daI_mu0hatdl*d3SIdl3)   !Correct
d4SIdl4 = PAR_/mu0hat1*alphaI**4*aI*cff*((1.-PAR_/mu0hat1*aI) &
    * ((1.-PAR_/mu0hat1*aI)**2 - aI*PAR_/mu0hat1)  &
    + aI*PAR_/mu0hat1*(2.*PAR_/mu0hat1*aI-3.))

d3muIhatdl3 = d2mu0hatSIdl2*cff1 + dmu0hatSIdl*4.*betamu + d2mu0hatdl2*dSIdl + 2.*dmu0hatdl*d2SIdl2 + mu0hat*d3SIdl3  !Correct

d4muIhatdl4 = mu0hat*d4SIdl4 + 4.*dmu0hatdl*d3SIdl3 + 6.*d2mu0hatdl2*d2SIdl2 + 4.*dSIdl*d3mu0hatdl3+SI*d4mu0hatdl4  !Correct

K0N = params(iKN)
Kn  = ScaleTrait(PMU, K0N, alphaK)
fN  = NO3/(NO3 + Kn)  !Nitrogen limitation index

!Add iron limitation:
if (DO_IRON) then
   K0Fe    = params(iKFe)
   alphaFe = alphaK
   KFe     = ScaleTrait(PMU, K0Fe, alphaFe)
   fFe     = Fe/(Fe + KFe)
endif

!Evaluate whether N or Fe is limiting (Liebig's law):
!Only needs to call MM once:
if (DO_IRON .and. (fFe < fN)) then  !Fe is limiting
   call MM(Fe,  K0Fe, alphaFe, PMU, fN, dfNdl, d2fNdl2, d3fNdl3, d4fNdl4)
else
   call MM(NO3, K0N,  alphaK,  PMU, fN, dfNdl, d2fNdl2, d3fNdl3, d4fNdl4)
endif 

! Phytoplankton growth rate at the mean size:
 muNet = mu0hat*SI*fN
 dmudl = dmu0hatSIdl*fN + mu0hatSI*dfNdl
d2mudl2=2.*dmu0hatSIdl*dfNdl+d2mu0hatSIdl2*fN+mu0hatSI*d2fNdl2 
d3mudl3=3.*(d2mu0hatSIdl2*dfNdl+dmu0hatSIdl*d2fNdl2) +fN*d3muIhatdl3 + mu0hatSI*d3fNdl3 !Correct
d4mudl4=4.*d3muIhatdl3*dfNdl+6.*d2mu0hatSIdl2*d2fNdl2+4.*dmu0hatSIdl*d3fNdl3 + fN*d4muIhatdl4 + mu0hatSI*d4fNdl4  !Correct

Qmin=0.06
Qmax=3.*Qmin

!N:C ratio at avg. size

cff1    = 1d0-Qmin/Qmax
cff     = 1d0-cff1*fN
QN      = Qmin/cff
dcffdl  = -cff1*dfNdl
dQNdL   = cff1/cff**2 * dfNdl * Qmin
d2QNdL2 = Qmin*cff1*(d2fNdl2/cff**2 - 2./cff**3*dfNdl*dcffdl)  !Correct


!Chl:C ratio at avg. size
cff     = (thetamax - thetamin)/PAR_
theta   = thetamin+muNet/aI*cff   !Unit: gChl/molC
dthetadl = cff*dY_Xdl(muNet, aI, dmudl, aI*alphaI)  !Correct

d2thetadl2 = cff*d2Y_Xdl2(muNet,aI,dmudl,aI*alphaI, d2mudl2, &
                                               aI*alphaI**2) !Correct

!! Calculate the community sinking rate of phytoplankton
!      ! Phytoplankton sinking rate at the average size
!      w_p = ScaleTrait(PMU,self%pars(iw_p),self%pars(ialphaW))
!      d2wpdl2 = self%pars(ialphaW)**2*w_p
!      self%w_pAvg = w_p+0.5*VAR*d2wpdl2 ! Sinking rate of phytoplankton
return
end subroutine

!Density function of normal distribution
pure real function normal(mean,var,l)
  implicit none
  real, intent(in) :: mean,var,l
  real, parameter  :: pi=3.1415926535897932384633
  normal = 1D0/sqrt(2D0*var*pi)*exp(-(l-mean)**2/var/2D0) 
  return
end function normal

!Mechaelis-Mention functions and derivatives
subroutine MM(N, K0, alphaK, L, fN, dfNdl, d2fNdl2, d3fNdl3, d4fNdl4)
use bio_MOD, only : ScaleTrait
implicit none
real,   intent(in)  :: N, K0, alphaK, L
real,   intent(out) :: fN, dfNdl, d2fNdl2, d3fNdl3, d4fNdl4
real                :: Kn

   ! Half saturation constant for growth at avg. size
   Kn = ScaleTrait(L, K0, alphaK)
   fN = N/(N + Kn) ! Nutrient limitation index at avg. size
   dfNdl = -alphaK*Kn*N/(N+Kn)**2
 d2fNdl2 = -alphaK**2*N*Kn*(1./(N+Kn)**2 - 2.*Kn/(N+Kn)**3)
 d3fNdl3 = alphaK**3*N*Kn*(2.*N*Kn-(Kn-N)**2)/(Kn+N)**4
 d4fNdl4 = alphaK**4*N*Kn*(11.*Kn*N*(N-Kn)+Kn**3-N**3)/(N+Kn)**5  !Correct
end subroutine
!
subroutine IRONCYCLE(Temp, DET, PP_NZ,PP_ND, PP_PN,PP_DZ, Fe_scav, DETFe, DFe)
use bio_MOD, only : TEMPBOL,Ez,dtdays,Fe_N
implicit none
real, intent(in)    :: Temp,DET, PP_NZ, PP_PN, PP_DZ, PP_ND
real, intent(inout) :: DFe, DETFe
real, intent(out)   :: Fe_scav     !Iron scavenging for diagonosis
real                :: keq, cff
real, parameter     :: Kscm= 5D-3  !Minimal scavenging rate (d-1)
real, parameter     :: Ksc = 3D-2  !Particle dependent scavenging rate (umolN-1 d-1)
!Iron:Nitrogen ratio, Aumont et al. (2003) set Fe/C = 4D-6 (mol:mol). 
!assume redfield ratio of C/N. Times 1000 to convert umol N to nmol Fe
real, parameter     :: lFe  = 0.6    ! Iron ligand concentration (nM). (TOM10 P. 19)
!dFedt = -phytouptake + Zoo excretion - scavenging + remineralization + dust deposition
! (From TOM10 and PISCES)

!lFe: total ligand conc. (Nikelsen et al. Geosci. Model. Dev. 2015)
!     When iron concentration is above 0.6 nM, it is scavenged by DET
!To close the iron cycle, need to model the Fe in Detritus
!Following TOM10Appendix, Eq. 46
!The equilibrium constant between free iron and ligands and organic complexes. 
keq  = 10**(17.27-1565.7/(273.15 + Temp))
!lFe: 0.6 nM
!Following TOM10Appendix, Eq. 45
!Iron scavenging rate = (Basal scavenging rate + particle asorbtion)*FEprime
cff     = 1D0+(lFe-DFe)*keq 

Fe_scav = (Kscm + Ksc*DET*TEMPBOL(Ez,Temp))             &
               * (-cff + sqrt(cff**2 + 4D0*DFe*keq))/2D0/keq  

cff = Fe_N*PP_ND !The regeneration flux from DETFe ==> DFe

DFe = DFe + cff + ((PP_NZ-PP_PN)*Fe_N - Fe_scav*dtdays)

! dDETFe/dt = Zooplankton defecation and mortality + scavenging - regeneration
DETFe = DETFe + PP_DZ*Fe_N + dtdays*Fe_scav - cff

! In this way, although we do not explicitly model the iron contents in PHY and ZOO, the total mass of iron is conserved (excluding dust deposition which should balance with detritus sinking)
end subroutine
