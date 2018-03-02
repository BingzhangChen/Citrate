Module BIO_MOD
USE    param_MOD
implicit none

! MPI variables:
integer            :: numtasks, taskid, ierr, MPIRUN

! Number of vertical layers
integer, parameter :: nlev        = 30  

! Options of biological models
integer, parameter :: EFTdisc     = 1
integer, parameter :: EFTcont     = 2
integer, parameter :: Geiderdisc  = 3
integer, parameter :: NPZDcont    = 4
integer, parameter :: EFTsimple   = 5
integer, parameter :: Geidersimple= 6
integer, parameter :: NPZDFix     = 7
integer, parameter :: NPZDdisc    = 8
integer, parameter :: EFTsimIRON  = 9
integer, parameter :: NPZDFixIRON = 10
integer, parameter :: GeidsimIRON = 11
integer, parameter :: NPZDdiscFe  = 12
integer, parameter :: EFTdiscFe   = 13
integer, parameter :: EFTDarwin   = 14
integer, parameter :: EFT2sp      = 15
integer, parameter :: NPZD2sp     = 16
integer, parameter :: NPPZDD      = 17
integer, parameter :: EFTPPDD     = 18
integer, parameter :: NPZDN2      = 19
integer, parameter :: CITRATE3    = 20
integer, parameter :: GeiderDroop = 21

! Parameters for phytoplankton size fractional Chl
real, parameter :: pi=3.1415926535897932384633D0
real, parameter :: eps     = 1d-30
real, parameter :: PMU_min = log(1d1*pi/6d0*0.6**3)
real, parameter :: PMU_max = log(1d1*pi/6d0*4d1**3)
real, parameter :: PMU_1   = log(1d1*pi/6d0)
real, parameter :: PMU_3   = log(1d1*pi/6d0*3d0**3)
real, parameter :: PMU_10  = log(1d1*pi/6d0*1d1**3)
integer  :: AllocateStatus
integer  :: N_MLD     ! vertical grid index at the bottom of MLD
integer  :: bot_bound = 1 ! Option of bottom boundary condition
integer, parameter :: Dirichlet      = 0
integer, parameter :: Neumann        = 1

real     :: Temp(nlev), PAR(nlev), dtdays, Ntot, PARavg, wstr0(1) 
real     :: DFe(nlev)                           ! Dissolved iron concentration
real     :: Z_r(1:nlev), Z_w(0:nlev), Hz(nlev)  ! Grid variables
real     :: I_zero
integer  :: NVAR, Nout, iZOO, iZOO2, iDET,iDET2,iDETp,iDETFe, iPMU, iVAR, iPO4,iDIA
integer  :: iMTo, iVTo, iMIo, iVIo
integer  :: NVsinkterms,NPHY, NPar
integer  :: oZOO, oZOO2,oDET, oPON, oFER, oZ2N, oD2N, oPHYt,oCHLt,oPPt,omuAvg
integer  :: oPO4, oPOP, oDIA, oDIAu,oDETp, oDET2, oDETFe
integer  :: oPMU, oVAR, oMTo, oVTo, oMIo, oVIo
integer  :: odmudl,odgdl,od2mu,od2gdl,odmudT,od2mudT2,odmudI,od2mudI2  
integer  :: oD_NO3,oD_ZOO,oD_ZOO2,oD_DET,oD_DET2,oD_fer
integer  :: oD_PMU,oD_VAR,oD_MTo,oD_MIo,oD_VTo,oD_VIo
integer  :: oPAR_,oD_DETp,oD_DETFe,oD_PO4,oD_DIA
integer  :: oMESg,oMESgMIC,odgdl1,odgdl2,od2gdl1,od2gdl2,odVAR
integer  :: oCHLs(4)   ! Four size fractions of CHL

! Indices for parameters used in DRAM
integer  :: imu0,iaI0,igmax,iKN,iKP,iKPHY,iKPnif,iLnifp,iKFe,iRDN_N,iRDN_P
integer  :: ialphamu,ibetamu,ialphaKN,irhom
integer  :: imu0B, iaI0B, iA0N2, iRL2,iKN2,ibI0B
integer  :: iVTR,iVTRL,iVTRT,iVTRI,ifer,od3mu,od4mu
integer  :: izetaN,izetaChl, iaI0_C 
integer  :: ialphaI,iA0N,ialphaA,ialphaG,ialphaK, ialphaFe
integer  :: iQ0N,ialphaQ,iPenfac,iLref,iwDET,irdN,imz
integer  :: ithetamin,iQNmin,itau,iwDET2,igb,oFescav,odstdep
integer, allocatable :: iPHY(:), oPHY(:),oTheta(:),oQN(:),oQp(:)
integer, allocatable :: iCHL(:), oCHL(:),omuNet(:),oD_PHY(:),oD_CHL(:)
integer, allocatable :: iPHYC(:), oPHYC(:),oD_PHYC(:)
integer, allocatable :: oGraz(:),oSI(:), oLno3(:)
integer, allocatable :: Windex(:)
real,    allocatable :: Vars(:,:),Varout(:,:), params(:)
character(LEN=6), allocatable :: ParamLabel(:)

! Fixed Model parameters:
real, parameter :: PMUmax =1.5D1, VARmax=50D0
real            :: Femin  =0.02,K0Fe  =0.8, alphaFe=0.14
real, parameter :: GGE    =0.3, unass =0.24
real, parameter :: Fe_N   =0.0265 ! Fe:Nitrogen molar ratio (nmol : umol)
real, parameter :: thetm  =0.65
real, parameter :: RMchl0 =0.1
real, parameter :: Ep     =0.5,  Ez   =0.65 
real            :: alphamu=0.2,  betamu=-0.01
real, parameter :: zetaChl=0.6,  zetaN =0.8
real, parameter :: thetamax = 0.63, thetamin=0.02

!Temperature senstivity tuned by the algorithm
real :: KFe    =0.08     !Unit: nM. Gregg (2003) used ~0.1 nM for phyto

!These two parametes also to be tuned by the algorithm

! Size and weights for discrete models
real, allocatable          :: PMU_(:)
real, allocatable          :: wtCHL(:,:)  ! weight for each size class
! Indices for state variables
integer, parameter :: iNO3=1,oTemp=1,oPAR=2,oAks=3,oDust=4,oNO3=1,ow=5
integer            :: nutrient_uptake=1, grazing_formulation=3   
logical, parameter :: kill_the_winner=.TRUE.
logical, public    :: N2fix          =.FALSE.
logical            :: do_IRON        =.FALSE.
logical            :: singlerun      =.FALSE.

! Output indices:
character(LEN=10), allocatable  :: Labelout(:)  

integer,   parameter :: namlst=8
character(LEN=3)     :: Stn(Nstn)

! Model options:
integer              :: Model_ID = NPZDFix

! Local variables:
real :: INGES,gbar,EGES,Zmort,RES
real :: pp_ZP, pp_NZ, pp_ND, pp_DZ, tf_p, tf_z

integer, parameter :: Yes = 1, No = 0

CONTAINS
!========================================================
! Calculate total nitrogen in the system
subroutine Cal_TN
implicit none
integer :: i,k

Ntot = 0d0

do k = 1,nlev
   do i = 1,iDET
      Ntot = Ntot + Vars(i,k) * Hz(k)
   enddo
enddo

end subroutine Cal_TN
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


END Module BIO_MOD
