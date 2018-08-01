MODULE PARAM_MOD
USE BIOCOM_MOD
implicit none

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
! Define tracer matrix:
real               :: Vars(NVAR,nlev) = 0d0

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

integer, parameter :: LINEAR              = 1       ! Linear zooplankton mortality and linear grazing function
integer, parameter :: QUADRATIC           = 2       ! Quadratic zooplankton mortality
integer, parameter :: H2T                 = 2       ! Holling Type II with threshold
integer, parameter :: H3                  = 3       ! Holling Type III

integer, parameter :: grazing_formulation = LINEAR  ! Type of grazing function
integer, parameter :: ZOO_MORT            = LINEAR  ! Type of zooplankton mortality
!  Maxima and minima  of parameter values 
real               :: MaxValue(NPar)   = 0d0, MinValue(NPar) = 0d0
END MODULE
