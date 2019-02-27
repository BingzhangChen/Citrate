MODULE PARAM_MOD
USE mpi
USE BIOCOM_MOD
implicit none

! Number of stations to run
integer, parameter      :: Nstn      = 1

! Station name:
character(5), parameter :: Stn(Nstn) = 'SEATS'

integer, parameter :: bot_bound  = Neumann ! Option of bottom boundary condition

! Current biological MODEL!
integer, parameter :: Model_ID   = NPZDFix

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
integer, parameter :: iPHY(1)= iNO3   + 1
integer, parameter :: iZOO   = iPHY(1)+ 1
integer, parameter :: iDET   = iZOO   + 1
integer, parameter :: NVAR   = iDET     ! Total number of biological tracers
integer            :: iPHYC(NPHY)=0     ! Needed for compilation
! Define tracer matrix:
real               :: Vars(NVAR,nlev) = 0d0

! Define the number of sinking tracers:
integer, parameter :: NVsinkterms =  1

! Define the index of sinking tracers in the Vars matrix: 
integer, parameter :: Windex(NVsinkterms) = [ iDET]

! Indices for output variables
integer, parameter :: oNO3   = 1
integer, parameter :: oPHY(1)= oNO3    + 1
integer, parameter :: oZOO   = oPHY(1) + 1
integer, parameter :: oDET   = oZOO    + 1
integer, parameter :: oCHLt  = oDET    + 1
integer, parameter :: oPPt   = oCHLt   + 1
integer, parameter :: oPON   = oPPt    + 1
integer, parameter :: oPAR_  = oPON    + 1
integer, parameter :: oSI    = oPAR_   + 1
integer, parameter :: oLno3  = oSI     + 1
integer, parameter :: omuNet = oLno3   + 1
integer, parameter :: oGraz  = omuNet  + 1
integer, parameter :: oZ2N   = oGraz   + 1
integer, parameter :: oD2N   = oZ2N    + 1
integer, parameter :: otheta = oD2N    + 1
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
integer, parameter :: iaI0_C  = imu0  + 1
integer, parameter :: iKN     = iaI0_C+ 1
integer, parameter :: iwDET   = iKN   + 1  ! Index for detrital sinking rate
integer, parameter :: igmax   = iwDET + 1  ! Zooplankton maximal grazing rate
integer, parameter :: imz     = igmax + 1  ! Zooplankton mortality rate
integer, parameter :: NPar    = imz        ! Total number of parameters
real               :: params(NPar)     = 0d0        ! Define parameters
character(LEN=8)   :: ParamLabel(NPar) = 'Unknown'  !Define parameter labels

!  Maxima and minima  of parameter values 
real               :: MaxValue(NPar)      = 0d0, MinValue(NPar) = 0d0
integer            :: iVPHY, iVNO3, ibeta, iCOVNP, iVZOO, iCOVPZ, oPO4, iCOVNZ
integer            :: iKFE
real               :: MaxN, MaxVar
integer, parameter :: grazing_formulation = 3
CONTAINS

SUBROUTINE choose_model
IMPLICIT NONE

if(taskid==0) write(6,*) 'Nutrient-Phytoplankton-Zooplankton-Detritus (NPZD) model selected!'

! Assign label names for Varout:
Labelout(oTemp      ) = 'Temp'
Labelout(oPAR       ) = 'PAR0'
Labelout(oAks       ) = 'Aks '
Labelout(oDust      ) = 'Dust'
Labelout(ow         ) = 'w   '
Labelout(iNO3   + ow) = 'NO3 '
Labelout(iPHY(1)+ ow) = 'PHY'
Labelout(iZOO   + ow) = 'ZOO'
Labelout(iDET   + ow) = 'DET'
Labelout(oSI    + ow) = 'SI'
Labelout(oLNO3  + ow) = 'LNO3'
Labelout(oQN    + ow) = 'QN'
Labelout(otheta + ow) = 'Theta'
Labelout(oPAR_  + ow) = 'PAR'
Labelout(omuNet + ow) = 'mu'
Labelout(oGraz  + ow) = 'Graz'
Labelout(oCHLt  + ow) = 'CHL'
Labelout(oPPt   + ow) = 'NPP'
Labelout(oPON   + ow) = 'PON'
Labelout(oQN    + ow) = 'QN'
Labelout(oZ2N   + ow) = 'ZTON'
Labelout(oD2N   + ow) = 'DTON'

DO i = 1, NVAR
   Labelout(oD_VARS(i) + ow) = 'D_'//trim(Labelout(i+ow))
ENDDO

! Write out Varout labels to check:
IF (TASKID == 0) THEN
   do i = 1, Nout+ow
      write(6,*) 'Labelout(',i,') = ',trim(Labelout(i))
   enddo
ENDIF

! Initialize parameters
if (taskid==0) write(6,'(I2,1x,A30)') NPar,'parameters to be estimated.'

ParamLabel(imu0)   = 'mu0hat'
ParamLabel(imz)    = 'mz'
ParamLabel(igmax)  = 'gmax'
ParamLabel(iaI0_C) = 'aI0_C'
ParamLabel(iKN)    = 'KN'
ParamLabel(iwDET)  = 'wPHY'

! Ranges of mu0 follow Chen & Laws L&O (2017)
MaxValue(imu0)   = 3.5
MinValue(imu0)   = 0.1
  params(imu0)   = 0.95


!Scaling factor of real gmax to mumax
MaxValue(igmax)    =  4.
MinValue(igmax)    =  0.1
  params(igmax)    =  2.

!Scaling factor of real mz to mu0 (assume linear mortality)
MaxValue(imz)      =  0.8
MinValue(imz)      =  0.01
  params(imz)      =  0.2

! KN:
! Fennel et al. (2006): 0.007~1.5
! Chai et al. (2002): 0.05~1
! Franks (2009): 0.005~3
MaxValue(iKN)    =  5.0
MinValue(iKN)    =  0.1
  params(iKN)    =  .4

MaxValue(iaI0_C) = 0.1
MinValue(iaI0_C) = 0.01
  params(iaI0_C) = 0.05

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
