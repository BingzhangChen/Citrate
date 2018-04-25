Module BIO_MOD
! This module contains some common functions and variables for all biological models
USE    param_MOD
implicit none

! MPI variables:
integer            :: numtasks, ierr, MPIRUN

! Parameters for phytoplankton size fractional Chl
real, parameter :: pi      = 3.1415926535897932384633D0
real, parameter :: eps     = 1d-20
real, parameter :: PMU_min = log(1d1*pi/6d0*0.6**3)
real, parameter :: PMU_max = log(1d1*pi/6d0*4d1**3)
real, parameter :: PMU_1   = log(1d1*pi/6d0)
real, parameter :: PMU_3   = log(1d1*pi/6d0*3d0**3)
real, parameter :: PMU_10  = log(1d1*pi/6d0*1d1**3)
integer  :: AllocateStatus
integer  :: N_MLD     ! vertical grid index at the bottom of MLD
real     :: Temp(nlev), PAR(nlev), dtdays, Ntot, PARavg, wstr0(1) 
real     :: DFe(nlev)                           ! Dissolved iron concentration
real     :: Z_r(1:nlev), Z_w(0:nlev), Hz(nlev)  ! Grid variables
real     :: I_zero
!integer  :: iZOO, iZOO2, iDET,iDET2,iDETp,iDETFe, iPMU, iVAR, iPO4,iDIA
!integer  :: iMTo, iVTo, iMIo, iVIo, iVPHY,iVNO3,iCOVNP
!integer  :: oVNO3,oVPHY,oCOVNP
!integer  :: oZOO, oZOO2,oDET, oPON, oFER, oZ2N, oD2N, oPHYt,oCHLt,oPPt,omuAvg
!integer  :: oPO4, oPOP, oDIA, oDIAu,oDETp, oDET2, oDETFe
!integer  :: oPMU, oVAR, oMTo, oVTo, oMIo, oVIo
!integer  :: odmudl,odgdl,od2mu,od2gdl,odmudT,od2mudT2,odmudI,od2mudI2  
!integer  :: oD_NO3,oD_ZOO,oD_ZOO2,oD_DET,oD_DET2,oD_fer
!integer  :: oD_PMU,oD_VAR,oD_MTo,oD_MIo,oD_VTo,oD_VIo
!integer  :: oPAR_,oD_DETp,oD_DETFe,oD_PO4,oD_DIA,oD_VPHY,oD_VNO3,oD_COVNP
!integer  :: oMESg,oMESgMIC,odgdl1,odgdl2,od2gdl1,od2gdl2,odVAR
!integer  :: oCHLs(4)   ! Four size fractions of CHL
!
!! Indices for parameters used in DRAM
!integer  :: imu0,iaI0,igmax,iKN,iKP,iKPHY,iKPnif,iLnifp,iKFe,iRDN_N,iRDN_P
!integer  :: ialphamu,ibetamu,ialphaKN,irhom,iDp, iIopt, ibeta
!integer  :: imu0B, iaI0B, iA0N2, iRL2,iKN2,ibI0B
!integer  :: iVTR,iVTRL,iVTRT,iVTRI,ifer,od3mu,od4mu
!integer  :: izetaN,izetaChl, iaI0_C 
!integer  :: ialphaI,iA0N,ialphaA,ialphaG,ialphaK, ialphaFe
!integer  :: iQ0N,ialphaQ,iPenfac,iLref,iwDET,irdN,imz
!integer  :: ithetamin,iQNmin,itau,iwDET2,igb,oFescav,odstdep
!integer, allocatable :: iPHY(:), oPHY(:),oTheta(:),oQN(:),oQp(:)
!integer, allocatable :: iCHL(:), oCHL(:),omuNet(:),oD_PHY(:),oD_CHL(:)
!integer, allocatable :: iPHYC(:), oPHYC(:),oD_PHYC(:)
!integer, allocatable :: oGraz(:),oSI(:), oLno3(:)
!integer, allocatable :: Windex(:)

! Some common Model parameters:
real, parameter :: PMUmax =1.5D1, VARmax=50D0
real            :: Femin  =0.02,K0Fe  =0.8, alphaFe=0.14
real, parameter :: GGE    =0.3, unass =0.24
real, parameter :: Fe_N   =0.0265 ! Fe:Nitrogen molar ratio (nmol : umol)
real, parameter :: RMchl0 =0.1
real, parameter :: Ep     =0.5,  Ez   =0.65 
real            :: alphamu=0.2,  betamu=-0.01
real, parameter :: zetaChl=0.6,  zetaN =0.8
real, parameter :: thetamax = 0.63, thetamin=0.02
integer, parameter :: nutrient_uptake=1, grazing_formulation=3   
real :: KFe    =0.08     !Unit: nM. Gregg (2003) used ~0.1 nM for phyto

! Size and weights for discrete models
real, allocatable          :: PMU_(:)
real, allocatable          :: wtCHL(:,:)  ! weight for each size class
logical, parameter :: kill_the_winner=.TRUE.
logical, public    :: N2fix          =.FALSE.
logical            :: singlerun      =.FALSE.
integer,   parameter :: namlst=8

! Some common local variables:
real :: INGES,gbar,EGES,Zmort,RES
real :: pp_ZP, pp_NZ, pp_ND, pp_DZ, tf_p, tf_z

integer, parameter :: Yes = 1, No = 0

CONTAINS
!========================================================
subroutine assign_PMU
implicit none
integer :: i
real    :: dx
integer :: AllocateStatus
! Calculate the average size 
allocate(PMU_(NPHY), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
PMU_(:)    = 0d0

allocate(wtCHL(4,NPHY), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

wtCHL(:,:) = 0d0

PMU_(1)= PMU_min
dx     = (PMU_max-PMU_min)/dble(NPHY-1)

do i=2,NPHY
   PMU_(i) = PMU_(i-1)+dx
enddo

! Calculate the weight for each size class 
! (largest size class being the first to be consistent with observational data)
do i=1,NPHY
   if (PMU_(i) .le. PMU_1) then
      wtCHL(4,i) = 1d0  ! < 1 um
   elseif (PMU_(i) .le. PMU_3) then
      wtCHL(3,i) = 1d0  ! 1-3 um
   elseif (PMU_(i) .le. PMU_10) then
      wtCHL(2,i) = 1d0  ! 3-10 um
   else
      wtCHL(1,i) = 1d0  ! >10 um
   endif
enddo
end subroutine assign_PMU
!========================================================
real function TEMPBOL(Ea,tC)
implicit none
!DESCRIPTION:
!The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
! tC: in situ temperature
! Tr: reference temperature
!
!INPUT PARAMETERS:
real, intent (in) :: Ea, tC
! boltzman constant constant [ eV /K ]
real, parameter   :: kb = 8.62d-5, Tr = 15D0

TEMPBOL = exp(-(Ea/kb)*(1D0/(273.15 + tC)-1D0/(273.15 + Tr)))
return 
end function TEMPBOL
!====================================================
real function ScaleTrait( logsize, star, alpha ) 
implicit none
real, intent(IN) :: logsize, star, alpha

! Calculate the size-scaled value of a trait
! for the given log (natural, base e) of cell volume as pi/6*ESD**3 (micrometers). 

ScaleTrait = star * exp( alpha * logsize )

return
end function ScaleTrait
!====================================================
real function PenAff( logsize, alpha, Pfac, lmin ) 
implicit none
real, intent(IN) :: logsize, alpha, Pfac, lmin 

!A 'penalty' function to reduce the value of affinity for nutrient at very small cell sizes
!in order to avoid modeling unrealistically small cell sizes.  This is needed because affnity
!increases with decreasing cell size, which means that under low-nutrient conditions, without
!such a penalty, unrealistically small cell sizes could be predicted.
!This penalty function becomes zero at logsize = lmin.   
   
  PenAff = 1.0 - exp(Pfac*alpha*(logsize - lmin))
end function PenAff
!====================================================
real function grazing(Hollingtype, Ksat, Prey)
implicit none
real,    intent(in) :: Ksat, Prey
integer, intent(in) :: Hollingtype
! kp relates to the capture coefficient
SELECT CASE(Hollingtype)
  ! Holling Type I
  case (1)
    grazing = min(Prey/2.0/Ksat,1.0)
  ! Holling Type II
  case (2)
    grazing = Prey/(Ksat + Prey)  
  ! Holling Type III
  case (3) 
    grazing = min(Prey*Prey/(Ksat*Ksat + Prey*Prey), 1D0)
 ! Ivlev
  case (4)
 !To be consistent with other functions  
    grazing = 1d0-exp(-log(2d0)*Prey/Ksat) 

END SELECT
return
end function grazing
!===================================================
!Functions deriving the first and second derivatives of X/Y ~ L 
pure real function dY_Xdl(Y, X, dYdl, dXdl) 
implicit none
real, intent(in)    :: Y, X, dYdl, dXdl
   dY_Xdl = dYdl/X - Y/X**2*dXdl
end function

pure real function d2Y_Xdl2(Y, X, dYdl, dXdl, d2Ydl2, d2Xdl2) 
implicit none
real, intent(in)    :: Y, X, dYdl, dXdl, d2Ydl2, d2Xdl2
 d2Y_Xdl2 = d2Ydl2/X - 2.*dYdl*dXdl/(X**2) - Y/(X**2)*d2Xdl2 +  &
                                          2.*Y/(X**3)*(dXdl**2)
end function

! Light ~ growth function from Edwards et al. L&O 2015
pure real function mu_Edwards2015(I_, lnIopt, mumax, alpha) 
implicit none
! lnIopt is log(Iopt)
real, intent(in)    :: I_, lnIopt, mumax, alpha
mu_Edwards2015 =  1./( I_/(alpha*exp(2.*lnIopt)) + 1./mumax -2./(alpha*exp(lnIopt))  &
            + 1./(I_*alpha)  ) 
end function mu_Edwards2015

! The subroutine for phytoplankton closure model
SUBROUTINE PHY_NPCLOSURE(NO3,PAR_,Temp_,PHY,VPHY,VNO3, COVNP,mu0hatSI, muNet,SI,theta,QN, Snp, Chl, NPP, SVPHY, SVNO3, SCOVNP)
implicit none
real, intent(in)  :: NO3, PAR_,Temp_,PHY, VPHY, VNO3, COVNP 
real, intent(out) :: muNet, theta, QN, SI, Snp, Chl, NPP, SVPHY, SVNO3, &
SCOVNP,mu0hatSI

! muNet: mean growth rate at <N>
real :: mu0hat
real :: alphaI, cff,cff1, eta, dYdN,d2YdN2, dEta_dN, d2ChldN2, Q
real :: KN, tf, fN, dmuQ_dN, d2NPPdN2
real, parameter :: Qmin = 0.06, Qmax=0.18

alphaI   = exp(params(iaI0_C))
tf       = TEMPBOL(Ep,Temp_)   !Temperature effect

!The temperature and light dependent component
mu0hat   = tf*exp(params(imu0))
mu0hatSI = mu_Edwards2015(PAR_, params(iIopt),mu0hat, alphaI) 
SI       = mu0hatSI/mu0hat
KN       = exp(params(iKN))
fN       = NO3/(NO3 + Kn)  !Nitrogen limitation index

! Phytoplankton growth rate at the mean nitrogen:
muNet = mu0hatSI*fN

! Snp: ensemble mean production (nitrogen based)
Snp   = PHY*muNet + mu0hatSI*(Kn*COVNP/(Kn+NO3)**2 - Kn*PHY*VNO3/(Kn+NO3)**3)

!N:C ratio at <N>
cff1  = 1.-Qmin/Qmax
cff   = 1.-cff1*fN
QN    = Qmin/cff
Q     = 1./QN

!Chl:C ratio at <N>
cff   = (thetamax - thetamin)/PAR_
theta = thetamin+muNet/alphaI*cff   !Unit: gChl/molC

!Chl:N ratio at <N>
eta      = theta*Q
dYdN     = Kn/(NO3 + Kn)**2
d2YdN2   = -dYdN*2./(NO3 + Kn) 
cff1     = 1./Qmax - 1./Qmin
dEta_dN  = dYdN*(theta*cff1 + Q*mu0hatSI/alphaI*cff)
d2ChldN2 = 2.*PHY*(dYdN**2*(cff1*cff*mu0hatSI/alphaI)-dEta_dN/(NO3+Kn))

!Ensemble mean Chl
Chl      = PHY*eta + .5*(2.*COVNP*dEta_dN + VNO3*d2ChldN2)
Chl      = max(Chl, eps)

dmuQ_dN  = dYdN*(muNet*cff1 + Q*mu0hatSI)
d2NPPdN2 = PHY*(d2YdN2*dmuQ_dN/dYdN + 2.*dYdN**2*mu0hatSI*cff1)

! NPP: ensemble mean carbon based primary production
NPP      = PHY*muNet*Q + .5*(2.*COVNP*dmuQ_dN + VNO3*d2NPPdN2)
NPP      = max(NPP, 0.)

! Calculate sources and sinks of variances of N, P, and covariance of NP
SVPHY    = 2.*mu0hatSI*(fN*VPHY         + Kn*PHY*       COVNP/(Kn+NO3)**2)
SVNO3    =-2.*mu0hatSI*(fN*COVNP        + Kn*PHY*        VNO3/(Kn+NO3)**2) 
SCOVNP   =    mu0hatSI*(fN*(COVNP-VPHY) + Kn*PHY*(VNO3-COVNP)/(Kn+NO3)**2) 
return
END SUBROUTINE PHY_NPCLOSURE
END Module BIO_MOD
