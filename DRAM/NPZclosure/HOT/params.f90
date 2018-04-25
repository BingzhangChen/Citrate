MODULE param_MOD
USE MPI
implicit none

! Number of stations to run
integer, parameter :: Nstn       = 1

! Station name:
character(len=3), parameter :: Stn(Nstn) = 'HOT'

! Number of vertical layers
integer, parameter :: nlev        = 40  

! Bottom boundary condition:
integer, parameter :: Dirichlet   = 0
integer, parameter :: Neumann     = 1
integer, parameter :: bot_bound   = 1 ! Option of bottom boundary condition

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
integer, parameter :: NPclosure   = 22
integer, parameter :: NPZclosure  = 23

! Current biological MODEL!
integer, parameter :: Model_ID    = NPZclosure

! Number of phytoplankton groups
integer, parameter :: NPHY       = 1

logical, parameter :: DO_IRON = .FALSE.  ! Whether involve iron or not

!  Indices for external forcing variables
integer, parameter :: etemp      = 1
integer, parameter :: eNO3       = 2
integer, parameter :: eAks       = 3
integer, parameter :: ew         = 4
integer, parameter :: ePAR       = 5
integer, parameter :: eDust      = 6
integer, parameter :: eFer       = 7
integer, parameter :: ePO4       = 8
integer, parameter :: ewstr      = 9
integer, parameter :: TNFo       = ewstr ! Total number of forcings
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

! Define tracer matrix:
real               :: Vars(NVAR,nlev) = 0d0

! Define the number of sinking tracers:
integer, parameter :: NVsinkterms =  2  ! Include PHY and VPHY

! Define the index of sinking tracers in the Vars matrix: 
integer, parameter :: Windex(NVsinkterms) = [iPHY(1), iVPHY]

! Indices for output variables
integer, parameter :: oTemp  = 1,oPAR=2,oAks=3,oDust=4,ow=5
integer, parameter :: oCHL   = iCOVNZ + 1
integer, parameter :: oNPP   = oCHL    + 1
integer, parameter :: oPAR_  = oNPP    + 1
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
integer, parameter :: iIopt   = imu0  +1
integer, parameter :: iaI0_C  = iIopt +1
integer, parameter :: iKN     = iaI0_C+1
integer, parameter :: iDp     = iKN   +1
integer, parameter :: iwDET   = iDp   +1  ! Index for phytoplankton sinking rate
integer, parameter :: ibeta   = iwDET +1  ! Beta: ratio of total variance to squares of mean concentration
integer, parameter :: igmax   = ibeta +1  ! Maximal grazing rate
integer, parameter :: imz     = igmax +1  ! Zooplankton mortality rate
integer, parameter :: NPar    = imz       ! Total number of parameters
real               :: params(NPar) = 0d0  ! Define parameters
character(LEN=8)   :: ParamLabel(NPar) = 'Unknown' !Define parameter labels

!  Maxima and minima  of parameter values 
real               :: MaxValue(NPar)   = 0d0, MinValue(NPar) = 0d0
integer            :: taskid

! Indices in other models (not used, but needed when compiling)
integer            :: iPHYC(NPHY), iDET, iZOO2, iPMU, iVAR, iMTo, iVTo, iMIo
integer            :: iVIo, ifer, iDETFe,iPO4,  iDIA, iDETp 

contains

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
Labelout(oCHL    + ow) = 'CHL'
Labelout(oNPP    + ow) = 'NPP'

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
if (taskid==0) write(6,'(I2,1x,A20)') NPar,'parameters to be estimated.'

ParamLabel(imu0)   = 'mu0hat'
ParamLabel(imz)    = 'mz'
ParamLabel(igmax)  = 'gmax'
ParamLabel(iaI0_C) = 'aI0_C'
ParamLabel(iKN)    = 'KN'
ParamLabel(iwDET)  = 'wPHY'
ParamLabel(ibeta)  = 'beta'
ParamLabel(iIopt)  = 'Iopt'
ParamLabel(iDp)    = 'DPHY'

MaxValue(imz)   =  0.2
MinValue(imz)   =  0.05
MaxValue(igmax) =  2d0
MinValue(igmax) =  0.5

! KN:
! Fennel et al. (2006): 0.007~1.5
! Chai et al. (2002): 0.05~1
! Franks (2009): 0.005~3
MaxValue(iKN) =  3.0
MinValue(iKN) =  0.05

! The ranges of Iopt and aI follow Edwards et al. L&O (2015)
MaxValue(iIopt)  = 2500.
MinValue(iIopt)  = 50.

MaxValue(ibeta)  = 7.0
MinValue(ibeta)  = 0.001

MaxValue(iDp)    = 0.5
MinValue(iDp)    = 0.01

MaxValue(iaI0_C) = 0.1
MinValue(iaI0_C) = 0.01

! Ranges of mu0 follow Chen & Laws L&O (2017)
MaxValue(imu0)   = 2.5
MinValue(imu0)   = 0.1

! Detritus sinking rate
! Fennel et al. (2006) gave range of 0.009-25 m/d
! Kishi et al. (2007) gave sinking rate of POC of 40 m/d
! But this is phytoplankton sinking rate, so based on Smayda (1970)
MaxValue(iwDET)  =  50.
MinValue(iwDET)  =  0.01

!Log transform:
do i = 1, NPar
   if (MaxValue(i) .le. 0d0) &
   stop "Parameter MaxValues uninitialized!"
   MaxValue(i) = log(MaxValue(i))
   MinValue(i) = log(MinValue(i))

   !Initial values of params randomly selected from the range between min and max
10 call random_number(cff)
   params(i)   = MinValue(i) + cff*(MaxValue(i) - MinValue(i))

   !To avoid zero in initial params
   if (params(i) == 0d0) goto 10

enddo
end subroutine choose_model
END MODULE
