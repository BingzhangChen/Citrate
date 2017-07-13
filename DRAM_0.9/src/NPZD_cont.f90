SUBROUTINE NPZD_CONT
use bio_MOD
implicit none

integer :: k,j,i
!INPUT PARAMETERS:
real :: tC,par_
!LOCAL VARIABLES of phytoplankton:
real :: PMUPHY,VARPHY,NO31,PHY1,ZOO1
real :: PMUPHY1,VARPHY1  !State variables
real :: PMU,VAR,PMU1,VAR1
real :: mu0,aI0
real :: QN,QP  ! cell quota related variables
real :: muNet,dmuNetdl,d2muNetdl2, d3mudl3, d4mudl4
real :: alphaK,alphaI, alphaG,SI,Lno3,KFe_
real :: muNet1, SI1, Lno31, QN1
real :: mu ! a scratch variable to temporally store phyto. growth rate
real :: theta,theta1,dgdlbar,d2gdl2bar
real :: NO3, Fe, PHY, ZOO, DET, DET1

!Declarations of zooplankton:
real :: dx  ! size interval
real :: gmax,RDN,mz,Kp,Cf  !Zooplankton variables
real :: pp_PN, PPpn
real :: Ptot, PCtot, Ptot1,CHLt,KN, Pl, dCHL,NPPt,rl
real :: pCHL(4) = 0D0
real :: VTR,Fe_scav
real :: x(M)
real,     external :: pnorm, normal
real,    parameter :: PMU0  = log(1d1)
real,    parameter :: bI0   = 0D0
real,    parameter :: eps   = 1d-20
real,    parameter :: Qpmin = 0.05/16d0
real,    parameter :: KPO4  = 0.5/16d0

integer, parameter :: M=50    !discretize the continous normal distribution

!-----------------------------------------------------------------------
VTR    = params(iVTR)
alphamu= params(ialphamu)
betamu = params(ibetamu)
alphaI = params(ialphaI)
alphaK = params(ialphaKN)
gmax   = params(igmax)
RDN    = 0.1
mz     = params(imz)
Kp     = 0.5
alphaG = 1D0+10**params(ialphaG)
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
   ZOO    = Vars(iZOO,    k)
   DET    = Vars(iDET,    k)
   PHY    = Vars(iPHY(1), k)
   Varout(oPON,k)=ZOO+DET+PHY
   PMUPHY = Vars(iPMU,    k)
   VARPHY = Vars(iVAR,    k)

   PMU    = PMUPHY/PHY
   PMU    = PMU - PMU0     !Correct PMU to the real one
   VAR    = VARPHY/PHY**2  !Correct VAR to the real one
   Fe     = max(Vars(ifer,k),Femin)

   call PHY_NPZDCONT(NO3,par_,tC,Fe,PMU,muNet,dmuNetdl,d2muNetdl2,d3mudl3,d4mudl4,SI,Lno3, theta, QN)

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
    
   ! Calculate the integration  
   x(1)  = PMU-3D0*sqrt(VAR)  ! Minimal size
   x(M)  = PMU+3D0*sqrt(VAR)  ! Maximal size
   dx    = (x(M) - x(1))/float(M-1)
   do i  = 2,M-1
     x(i) = x(i-1) + dx  !Log size of each class
   enddo

   ! Calculate total phytoplankton palatable prey (Ptot,Eq. 13):
   ! Ptot roughly equals to Phy
   Cf    = params(igb)
   Ptot  = 0D0  ! Total palatable prey
   NPPt  = 0D0
   PCtot = 0D0  ! Phytoplankton carbon
   Ptot1 = 0D0
   CHLt  = 0D0
   pCHL(:)=0D0  ! Size fractionated Chl

   do i=1,M
     !Calculate Chl in this size class

     !Pertain only meaningful size class
     if (x(i) .gt. -3 .and. x(i) .le. 15.) then 
       mu0=params(imu0)  *exp(alphamu*x(i)+betamu*x(i)**2)
       aI0=params(iaI0_C)*exp(alphaI*x(i))
       KN =params(iKN)   *exp(alphaK*x(i))
      KFe_=params(iKFe)  *exp(alphaFe*x(i))

   ! MONOD(Temp, PAR,NO3,PO4,mu0,Qnmin,Qpmin,aI0,bI0,KN,KP,DFe,KFe,muNet,QN,QP,theta,SI,Lno3)
      call MONOD(tC,par_,NO3,1d0,mu0,params(iQ0N),Qpmin,aI0,bI0,KN,KPO4,&
                 Fe,KFe_,mu, QN, QP, theta, SI, Lno3)
   
      !Probability function
      Pl=PHY*normal(PMU,VAR,x(i))

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

       ! Total Phytoplankton carbon
       PCtot= PCtot + Pl*dx/QN

       ! Total NPP (use light as the in situ incubation):
       call MONOD(tC,PAR(k),NO3,1d0,mu0,params(iQ0N),Qpmin,aI0,bI0,KN,KPO4,&
                  Fe,KFe_,muNet1, QN1,QP, theta1, SI1, Lno31)

       NPPt = NPPt + Pl*dx*muNet1/QN

       !Feeding preference
       rl   = exp(Cf*x(i))
       Ptot = Ptot +Pl*dx*rl   ! Total palatable prey

    !   Considering killing-the-winner strategy
       Ptot1= Ptot1+Pl**alphaG*dx*rl
     endif
   enddo
   do j = 1, 4
      Varout(oCHLs(j),k)=pCHL(j)/max(CHLt,eps) !Get the percentage of size-fractionated CHL
   enddo

   PCtot           = max(PCtot, eps)
   Varout(oCHLt,k) = CHLt
   Varout(oPPt, k) = NPPt*12d0 ! Correct the unit to ug C/L

! The grazing dependence on total prey
   gbar = grazing(grazing_formulation,Kp,PHY)

   !Zooplankton per capita total ingestion rate
   INGES = tf_z*gmax*gbar

   !Zooplankton excretion rate (-> DOM)
   RES = INGES*(1D0-GGE-unass)

   !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass

! Grazing rate on the mean size of PHY (specific to N-based Phy biomass, unit: d-1) (Eq. 12)
   rl   = exp(Cf*PMU)

   ! The probability density at avg. size
   Pl   = normal(PMU,VAR,PMU)*PHY

   ! Grazing rate on the mean size
   gbar = INGES*ZOO*rl*Pl**(alphaG-1D0)/max(Ptot1,eps) 

! First derivative of grazing rate evaluated at the mean size (negative)

   dgdlbar  = Cf*gbar 
   d2gdl2bar= ((1D0-alphaG)/VAR+Cf*Cf)*gbar

!!End of zooplankton section
!=============================================================
!! Solve ODE functions:
   Zmort = ZOO*ZOO*mz*tf_z  !Mortality term for ZOO

   VAR1  = VAR+dtdays   &
   * (VAR*(VAR*(d2muNetdl2-d2gdl2bar+VTR*d4mudl4)-5d0*VTR*d2muNetdl2)+2d0*VTR*muNet)
  
   !Eq. 18, Update PMU:
   PMU = PMU + PMU0  ! Restore to positive values
   PMU1= PMU + dtdays*(VAR*(dmuNetdl-dgdlbar+VTR*d3mudl3) - 3d0*VTR*dmuNetdl)

!   Constrain the PMU and VAR: 
!   PMU1 = min(max(PMU1,0.01),PMUmax)
!   VAR1 = min(max(VAR1,0.01),VARmax)
          
!All rates have been multiplied by dtdays to get the real rate correponding to the actual time step
   PP_ND= dtdays*RDN*DET*tf_z   
   PP_NZ= ZOO*RES*dtdays        
   PP_DZ= (ZOO*EGES+Zmort)*dtdays 
   PP_ZP= ZOO*INGES*dtdays      
   PPpn = PHY*dtdays*(muNet+0.5*VAR*(d2muNetdl2+VTR*d4mudl4)-1.5*VTR*d2muNetdl2)
   PP_PN= max(PPpn,0D0)
!
   DET1 = DET + PP_DZ         - PP_ND 
   NO31 = NO3 + PP_ND + PP_NZ - PP_PN
   PHY1 = PHY + PP_PN         - PP_ZP
   ZOO1 = ZOO + PP_ZP         - PP_NZ - PP_DZ
   
   call IRONCYCLE(tC, DET, Fe ,PP_ND,PP_NZ,PP_PN,Fe_scav)

   PMUPHY1 = PMUPHY + PHY*(PMU1-PMU) + PMU*(PHY1-PHY)

   ! Use VAR*PHY**2 as the tracer
   VARPHY1 = VARPHY + PHY**2*(VAR1-VAR) + 2.*PHY*VAR*(PHY1-PHY)

   Varout(oNO3,k)      = NO31
   Varout(oPHY(1),k)   = PHY1
   Varout(oPHYt,  k)   = PHY1
   Varout(oZOO,k)      = ZOO1
   Varout(oDET,k)      = DET1
   Varout(ofer,k)      = Fe
   Varout(oPMU,k)      = PMUPHY1
   Varout(oVAR,k)      = VARPHY1
   Varout(oTheta(1),k) = CHLt/PCtot 
   Varout(oQN(1)   ,k) = PHY/PCtot
   Varout(omuNet(1),k) = muNet               !Growth rate of mean size
   Varout(omuAvg,   k) = PP_PN/dtdays/PHY    !Avg. growth rate of the total community
   Varout(oGraz(1) ,k) = PP_ZP/dtdays/PHY
   Varout(oZ2N     ,k) = PP_NZ/dtdays
   Varout(oD2N     ,k) = PP_ND/dtdays
   Varout(odmudl   ,k) = dmuNetdl
   Varout(od2mu    ,k) = d2muNetdl2
   Varout(od3mu    ,k) = d3mudl3 
   Varout(od4mu    ,k) = d4mudl4
ENDDO
END SUBROUTINE NPZD_CONT 

! The subroutine only for phytoplankton
subroutine PHY_NPZDCONT(NO3,PAR_,Temp_,Fe, PMU, muNet,dmudl,d2mudl2,d3mudl3,d4mudl4,SI,fN, theta, QN)
use bio_MOD, only : ScaleTrait, TEMPBOL, params
use bio_MOD, only : ialphaI, ialphaKN,imu0, ialphamu,ibetamu, iaI0_C, iKN, iQ0N  
use bio_MOD, only : Ep, K0Fe, alphaFe, KFe, do_IRON,betamu, alphamu, iKFe, ialphaFe
implicit none
real, intent(in)  :: PMU, NO3, PAR_,Temp_, Fe 
real, intent(out) :: muNet, dmudl, d2mudl2, theta, QN, SI, fN
real, intent(out) :: d3mudl3, d4mudl4
real :: muNet_, dmudl_, d2mudl2_,d3mudl3_
real :: alphaK, alphaI, mu0hat
real :: dmu0hatdl, d2mu0hatdl2, aI, cff,daI_mu0hatdl
real :: daI_mu0hat2dl
real :: d3mu0hatdl3, d4mu0hatdl4
real :: dSIdl
real :: mu0hatSI 
real :: dmu0hatSIdl, dmu0hat_aIdl,d2mu0hat_aIdl2
real :: d2aI_mu0hatdl2,d3aI_mu0hatdl3,d2SIdl2,d2mu0hatSIdl2
real :: d3SIdl3,daI_mu0hat3dl,d2aI_mu0hat2dl2,d4SIdl4
real :: daI_mu0hat4dl,d2aI_mu0hat3dl2,d3aI_mu0hat2dl3, d4aI_mu0hatdl4
real :: d3muIhatdl3, d4muIhatdl4
real :: Kn,K0N,dfNdl,d2fNdl2,d3fNdl3, d4fNdl4, Qmin,Qmax,tf
real :: fFe, dfFedl, d2fFedl2, d3fFedl3, d4fFedl4
real, parameter :: thetamin = 0.02, thetamax = 0.5

alphaK     =params(ialphaKN) 
alphamu    =params(ialphamu)
betamu     =params(ibetamu)
alphaI     =params(ialphaI)
tf         =TEMPBOL(Ep,Temp_)
mu0hat     =tf*params(imu0)*exp(alphamu*PMU + betamu*PMU**2)
dmu0hatdl  =mu0hat*(alphamu + 2D0 * betamu * PMU)
d2mu0hatdl2=mu0hat*2D0*betamu+mu0hat*(alphamu+2D0*betamu*PMU)**2
d3mu0hatdl3=(2.*betamu+(alphamu+2.*betamu*PMU)**2)*dmu0hatdl+4.*betamu*mu0hat*(alphamu+2.*betamu*PMU)

d4mu0hatdl4=dmu0hatdl*8.*betamu*(alphamu+2.*betamu*PMU)+(2.*betamu+(alphamu+2.*betamu*PMU)**2)*d2mu0hatdl2+8.*betamu**2*mu0hat

! Initial slope of P-I curve
aI=ScaleTrait(PMU, params(iaI0_C), alphaI)

!The light limitation index (SI)
cff=exp(-aI*PAR_/mu0hat)
SI =1D0-cff

daI_mu0hatdl=aI*(alphaI-alphamu-2D0*betamu*PMU)/mu0hat

dSIdl=cff*PAR_*daI_mu0hatdl

mu0hatSI=mu0hat*SI
dmu0hatSIdl = mu0hatSI*(alphamu+2D0*betamu*PMU) + mu0hat*dSIdl
dmu0hat_aIdl= (dmu0hatdl - alphaI*mu0hat)/aI
d2mu0hat_aIdl2 = mu0hat/aI*2D0*betamu+(alphamu-alphaI+2D0*betamu*PMU)*dmu0hat_aIdl

daI_mu0hat2dl  = daI_mu0hatdl/mu0hat - aI/mu0hat**3*dmu0hatdl !Correct

d2aI_mu0hatdl2 = alphaI*daI_mu0hatdl - aI/mu0hat**2*d2mu0hatdl2 - daI_mu0hat2dl*dmu0hatdl

d3aI_mu0hatdl3 = d2aI_mu0hatdl2*(alphaI-alphamu-2d0*betamu*PMU)-4d0*betamu*daI_mu0hatdl

d2SIdl2 = -PAR_*dSIdl*daI_mu0hatdl+PAR_*cff*d2aI_mu0hatdl2

d2mu0hatSIdl2 = (alphamu+2d0*betamu*PMU)*dmu0hatSIdl + 2d0*betamu*mu0hatSI + mu0hat*(alphamu+2d0*betamu*PMU)*dSIdl + mu0hat*d2SIdl2  !Correct

d3SIdl3 = PAR_*(-2.*dSIdl*d2aI_mu0hatdl2 - d2SIdl2*daI_mu0hatdl + (1.-SI)*d3aI_mu0hatdl3)  !Correct

daI_mu0hat3dl  = daI_mu0hat2dl/mu0hat - aI*dmu0hatdl/mu0hat**4 !Correct

d2aI_mu0hat2dl2= d2aI_mu0hatdl2/mu0hat - daI_mu0hatdl/mu0hat**2*dmu0hatdl - daI_mu0hat3dl * dmu0hatdl - aI/mu0hat**3*d2mu0hatdl2 !Correct 

d3aI_mu0hatdl3 = alphaI*d2aI_mu0hatdl2 - aI/mu0hat**2*d3mu0hatdl3 -2d0*daI_mu0hat2dl*d2mu0hatdl2 -d2aI_mu0hat2dl2*dmu0hatdl !Correct

daI_mu0hat4dl  = daI_mu0hat3dl/mu0hat - aI/mu0hat**5*dmu0hatdl  !Correct

d2aI_mu0hat3dl2= d2aI_mu0hat2dl2/mu0hat - daI_mu0hat2dl/mu0hat**2*dmu0hatdl - aI/mu0hat**4*d2mu0hatdl2 - daI_mu0hat4dl*dmu0hatdl   !Correct

d3aI_mu0hat2dl3= (d3aI_mu0hatdl3*mu0hat - daI_mu0hatdl*d2mu0hatdl2)/mu0hat**2 -2d0/mu0hat**3*(d2aI_mu0hatdl2*mu0hat - daI_mu0hatdl*dmu0hatdl)*dmu0hatdl - (aI/mu0hat**3*d3mu0hatdl3+2d0*daI_mu0hat3dl*d2mu0hatdl2+d2aI_mu0hat3dl2*dmu0hatdl)  !Correct

d4aI_mu0hatdl4 = alphaI*d3aI_mu0hatdl3 - (aI/mu0hat**2*d4mu0hatdl4 + 3.* daI_mu0hat2dl*d3mu0hatdl3+3.*d2aI_mu0hat2dl2 * d2mu0hatdl2)- d3aI_mu0hat2dl3*dmu0hatdl  !Correct

d4SIdl4 = PAR_*((1.-SI)*d4aI_mu0hatdl4 - 3.*dSIdl*d3aI_mu0hatdl3 - 3.*d2aI_mu0hatdl2*d2SIdl2-daI_mu0hatdl*d3SIdl3)   !Correct

d3muIhatdl3 = d2mu0hatSIdl2*(alphamu + 2.*betamu*PMU) + dmu0hatSIdl*4.*betamu + d2mu0hatdl2*dSIdl + 2.*dmu0hatdl*d2SIdl2 + mu0hat*d3SIdl3  !Correct

d4muIhatdl4 = mu0hat*d4SIdl4 + 4.*dmu0hatdl*d3SIdl3 + 6.*d2mu0hatdl2*d2SIdl2 + 4.*dSIdl*d3mu0hatdl3+SI*d4mu0hatdl4  !Correct

K0N = params(iKN)
call MM(NO3, K0N, alphaK, PMU, fN, dfNdl, d2fNdl2, d3fNdl3, d4fNdl4)
Kn = ScaleTrait(PMU, K0N, alphaK)

! Phytoplankton growth rate at the mean size:
 muNet = mu0hat*SI*fN
 dmudl = dmu0hatSIdl*fN - mu0hatSI*alphaK*Kn*NO3/(NO3 + Kn)**2
d2mudl2=2*dmu0hatSIdl*dfNdl+d2mu0hatSIdl2*fN+mu0hatSI*d2fNdl2 
d3mudl3=3.*(d2mu0hatSIdl2*dfNdl+dmu0hatSIdl*d2fNdl2) +fN*d3muIhatdl3 + mu0hatSI*d3fNdl3 !Correct
d4mudl4=4.*d3muIhatdl3*dfNdl+6.*d2mu0hatSIdl2*d2fNdl2+4.*dmu0hatSIdl*d3fNdl3 + fN*d4muIhatdl4 + mu0hatSI*d4fNdl4  !Correct

!Add iron limitation:
if (DO_IRON) then
   K0Fe    = params(iKFe)
   alphaFe = params(ialphaFe)
   call MM(Fe, K0Fe, alphaFe, PMU, fFe, dfFedl, d2fFedl2, d3fFedl3, d4fFedl4)
   KFe     = ScaleTrait(PMU, K0Fe, alphaFe)
   muNet_  = muNet*fFe
   dmudl_  = muNet*dfFedl   + dmudl*fFe
   d2mudl2_= muNet*d2fFedl2 + 2.*dmudl*dfFedl + fFe*d2mudl2
   d3mudl3_= muNet*d3fFedl3 + 3.*(dmudl*d2fFedl2+d2mudl2*dfFedl)+fFe*d3mudl3
   d4mudl4 = muNet*d4fFedl4 + 4.*dmudl*d3fFedl3 + 6.*d2mudl2*d2fFedl2 + 4.*d3mudl3*dfFedl + fFe*d4mudl4  
   muNet   = muNet_
   dmudl   = dmudl_
   d2mudl2 = d2mudl2_
   d3mudl3 = d3mudl3_
endif  !All correct
Qmin=params(iQ0N)
Qmax=3D0*params(iQ0N)

!N:C ratio at avg. size
 QN  =Qmin/(1D0-(1D0-Qmin/Qmax)*fN)
theta=thetamin+muNet/PAR_/aI*(thetamax-thetamin)   !Unit: gChl/molC

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
 d2fNdl2 = -alphaK**2*N*Kn*(1/(N+Kn)**2 - 2*Kn/(N+Kn)**3)
 d3fNdl3 = alphaK**3*N*Kn*(2.*N*Kn-(Kn-N)**2)/(Kn+N)**4
 d4fNdl4 = alphaK**4*N*Kn*(11.*Kn*N*(N-Kn)+Kn**3-N**3)/(N+Kn)**5  !Correct
end subroutine
!
subroutine IRONCYCLE(Temp, DET, DFe, PP_ND, PP_NZ,PP_PN,Fe_scav)
use bio_MOD, only : Femin,TEMPBOL,Ez,dtdays
implicit none
real, intent(in)    :: Temp,DET, PP_ND, PP_NZ, PP_PN
real, intent(inout) :: DFe
real, intent(out)   :: Fe_scav
real                :: keq, cff
real, parameter     :: Kscm= 3D-5  !Minimal scavenging rate
real, parameter     :: Ksc = 3D-2  !Particle dependent scavenging rate (umolN-1 d-1)
!Iron:Nitrogen ratio, Aumont et al. (2003) set Fe/C = 4D-6 (mol:mol). 
!assume redfield ratio of C/N. Times 1000 to convert umol N to nmol Fe
real, parameter     :: Fe_N = 0.0265 ! Fe:Nitrogen molar ratio
real, parameter     :: lFe  = 0.6    ! Iron ligand concentration (nM). (TOM10 P. 19)

DFe  = max(DFe,Femin)

!dFedt = -phytouptake + Zoo excretion - scavenging + remineralization + dust deposition
! (From TOM10 and PISCES)

!lFe: total ligand conc. (Nikelsen et al. Geosci. Model. Dev. 2015)
!     When iron concentration is above 0.6 nM, it is scavenged by DET
!Following TOM10Appendix, Eq. 46
!The equilibrium constant between free iron and ligands and organic complexes. 
keq  = 10**(17.27-1565.7/(273.15 + Temp))
!lFe: 0.6 nM
!Following TOM10Appendix, Eq. 45
!Iron scavenging rate = (Basal scavenging rate + particle asorbtion)*FEprime
cff     = 1D0+(lFe-DFe)*keq 
Fe_scav = KScm + Ksc*DET*TEMPBOL(Ez,Temp)             &
               * (-cff + sqrt(cff**2 + 4.*DFe*keq))/2D0/keq  

DFe = DFe + ((PP_ND+PP_NZ-PP_PN)*Fe_N- Fe_scav*dtdays)

end subroutine


