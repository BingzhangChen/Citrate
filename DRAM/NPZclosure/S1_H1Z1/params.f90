MODULE PARAM_MOD
USE MPI
USE BIOCOM_MOD
implicit none

! Number of stations to run
integer, parameter      :: Nstn      = 1

! Station name:
character(3), parameter :: Stn(Nstn) = 'S1'

integer, parameter :: bot_bound  = Neumann ! Option of bottom boundary condition

! Current biological MODEL!
integer, parameter :: Model_ID   = NPZclosure

! Number of phytoplankton groups
integer, parameter :: NPHY       = 1
logical, parameter :: DO_IRON    = .FALSE.  ! Whether involve iron or not
integer, parameter :: N_fer      = 33    ! Number of vertical layers for Iron
! Total of observation times in forcing data
character(LEN=5), parameter :: LabelForc(TNFo) &
     = (/'temp ','NO3  ','Aks  ','wROMS','par  ', 'solfe','fer  ', 'PO4  ', 'wstr '/)
!! End of definding external forcing variables

! Indices for state variables
integer, parameter :: iNO3   = 1
integer, parameter :: iPHY(1)= 2
integer, parameter :: iZOO   = 3
integer, parameter :: iVNO3  = 4
integer, parameter :: iVPHY  = 5
integer, parameter :: iVZOO  = 6
integer, parameter :: iCOVNP = 7
integer, parameter :: iCOVPZ = 8
integer, parameter :: iCOVNZ = 9
integer, parameter :: NVAR   = iCOVNZ   ! Total number of biological tracers
integer            :: iPHYC(NPHY)=0     ! Needed for compilation
! Define tracer matrix:
real               :: Vars(NVAR,nlev) = 0d0

! Define the number of sinking tracers:
integer, parameter :: NVsinkterms =  2  ! Include PHY and VPHY

! Define the index of sinking tracers in the Vars matrix: 
integer, parameter :: Windex(NVsinkterms) = [iPHY(1), iVPHY]

! Indices for output variables
integer, parameter :: oCHLt  = iCOVNZ  + 1
integer, parameter :: oNPP   = oCHLt   + 1
integer, parameter :: oPON   = oNPP    + 1
integer, parameter :: oPAR_  = oPON    + 1
integer, parameter :: omuNet = oPAR_   + 1
integer, parameter :: oSI    = omuNet  + 1
integer, parameter :: otheta = oSI     + 1
integer, parameter :: oQN    = otheta  + 1

! The diffusion output order must be consistent with the tracer order!
integer, private   :: i
integer, parameter :: oD_VARS(NVAR) = [(oQN + i, i = 1, NVAR)]
integer, parameter :: Nout   = oD_VARS(NVAR)

! Initialize Varout matrix:
real               :: Varout(Nout, nlev) = 0d0

! Initialize labels for Varout:
character(LEN=10)  :: Labelout(Nout+ow)  = 'Unknown'
! Define parameter indices:
integer, parameter :: imu0    = 1
integer, parameter :: iIopt   = imu0  + 1
integer, parameter :: iaI0_C  = iIopt + 1
integer, parameter :: iKN     = iaI0_C+ 1
integer, parameter :: iDp     = iKN   + 1
integer, parameter :: iwDET   = iDp   + 1  ! Index for phytoplankton sinking rate
integer, parameter :: ibeta   = iwDET + 1  ! Beta: ratio of total variance to squares of mean concentration
integer, parameter :: iKPHY   = ibeta + 1  ! Maximal grazing rate
integer, parameter :: igmax   = iKPHY + 1  ! Zooplankton maximal grazing rate
integer, parameter :: imz     = igmax + 1  ! Zooplankton mortality rate
integer, parameter :: iVPHY0  = imz   + 1  ! Initial fraction of VPHY over B
integer, parameter :: iVNO30  = iVPHY0+ 1  ! Initial fraction of VNO3 over B
integer, parameter :: NPar    = iVNO30     ! Total number of parameters
real               :: params(NPar)     = 0d0  ! Define parameters
character(LEN=8)   :: ParamLabel(NPar) = 'Unknown' !Define parameter labels

integer, parameter :: LINEAR              = 1       ! Linear zooplankton mortality and linear grazing function
integer, parameter :: QUADRATIC           = 2       ! Quadratic zooplankton mortality
integer, parameter :: H2T                 = 2       ! Holling Type II with threshold
integer, parameter :: H3                  = 3       ! Holling Type III

integer, parameter :: grazing_formulation = LINEAR     ! Type of grazing function
integer, parameter :: ZOO_MORT            = LINEAR  ! Type of zooplankton mortality
!  Maxima and minima  of parameter values 
real               :: MaxValue(NPar)   = 0d0, MinValue(NPar) = 0d0
CONTAINS

SUBROUTINE choose_model
implicit none
real cff !a scratch variable

if(taskid==0) write(6,*) 'Nutrient-Phytoplankton-Zooplankton (NPZ) closure model selected!'

! Assign label names for Varout:
Labelout(oTemp       ) = 'Temp'
Labelout(oPAR        ) = 'PAR0'
Labelout(oAks        ) = 'Aks '
Labelout(oDust       ) = 'Dust'
Labelout(ow          ) = 'w   '
Labelout(iNO3    + ow) = 'NO3 '
Labelout(iPHY    + ow) = 'PHY'
Labelout(iZOO    + ow) = 'ZOO'
Labelout(iVPHY   + ow) = 'VPHY'
Labelout(iVNO3   + ow) = 'VNO3'
Labelout(iVZOO   + ow) = 'VZOO'
Labelout(iCOVNP  + ow) = 'COVNP'
Labelout(iCOVPZ  + ow) = 'COVPZ'
Labelout(iCOVNZ  + ow) = 'COVNZ'
Labelout(oSI     + ow) = 'SI'
Labelout(oQN     + ow) = 'QN'
Labelout(otheta  + ow) = 'Theta'
Labelout(oPAR_   + ow) = 'PAR'
Labelout(omuNet  + ow) = 'mu'
Labelout(oCHLt   + ow) = 'CHL'
Labelout(oNPP    + ow) = 'NPP'
Labelout(oPON    + ow) = 'PON'

DO i = 1, NVAR
   Labelout(oD_VARS(i) + ow) = 'D_'//trim(Labelout(i+ow))
ENDDO

! Write out Varout labels to check:
if (taskid == 0) then
   do i = 1, Nout+ow
      write(6,*) 'Labelout(',i,') = ',trim(Labelout(i))
   enddo
endif

! Initialize parameters
if (taskid==0) write(6,'(I2,1x,A30)') NPar,'parameters to be estimated.'

ParamLabel(imu0)   = 'mu0hat'
ParamLabel(imz)    = 'mz'
ParamLabel(iVPHY0) = 'VPHY0'
ParamLabel(iVNO30) = 'VNO30'
ParamLabel(iKPHY)  = 'KPHY'
ParamLabel(igmax)  = 'gmax'
ParamLabel(iaI0_C) = 'aI0_C'
ParamLabel(iKN)    = 'KN'
ParamLabel(iwDET)  = 'wPHY'
ParamLabel(ibeta)  = 'beta'
ParamLabel(iIopt)  = 'Iopt'
ParamLabel(iDp)    = 'DPHY'

!Scaling factor of real gmax to mumax
MaxValue(igmax)    =  10.
MinValue(igmax)    =  0.2
  params(igmax)    =  2.0

!Scaling factor of real mz to gmax (assume linear mortality)
MaxValue(imz)      =  0.8
MinValue(imz)      =  0.01
  params(imz)      =  0.1

MaxValue(iKPHY)    =  2d0
MinValue(iKPHY)    =  0.02
  params(iKPHY)    =  0.5

MaxValue(iVPHY0)   =  0.35
MinValue(iVPHY0)   =  0.01
  params(iVPHY0)   =  0.3

MaxValue(iVNO30)   =  0.35
MinValue(iVNO30)   =  0.01
  params(iVNO30)   =  0.3

! KN:
! Fennel et al. (2006): 0.007~1.5
! Chai et al. (2002): 0.05~1
! Franks (2009): 0.005~3
MaxValue(iKN)    =  3.0
MinValue(iKN)    =  0.05
  params(iKN)    =  0.2

! The ranges of Iopt and aI follow Edwards et al. L&O (2015)
MaxValue(iIopt)  = 2500.
MinValue(iIopt)  = 50.
  params(iIopt)  = 1d3

MaxValue(ibeta)  = 1d1
MinValue(ibeta)  = 0.001
  params(ibeta)  = .1

! The ratio of phytoplankton death rate to mumax
MaxValue(iDp)    = 0.9
MinValue(iDp)    = 0.01
  params(iDp)    = 0.1

MaxValue(iaI0_C) = 0.1
MinValue(iaI0_C) = 0.01
  params(iaI0_C) = 0.05

! Ranges of mu0 follow Chen & Laws L&O (2017)
MaxValue(imu0)   = 2.5
MinValue(imu0)   = 0.1
  params(imu0)   = 0.8

! Detritus sinking rate
! Fennel et al. (2006) gave range of 0.009-25 m/d
! Kishi et al. (2007) gave sinking rate of POC of 40 m/d
! But this is phytoplankton sinking rate, so based on Smayda (1970)
MaxValue(iwDET)  =  40.
MinValue(iwDET)  =  0.01
  params(iwDET)  =  0.1

!Log transform:
do i = 1, NPar
   if (MaxValue(i) .le. 0d0) &
   stop "Parameter MaxValues uninitialized!"
   MaxValue(i) = log(MaxValue(i))
   MinValue(i) = log(MinValue(i))
   params(i)   = log(params(i))
enddo
end subroutine choose_model
END MODULE
