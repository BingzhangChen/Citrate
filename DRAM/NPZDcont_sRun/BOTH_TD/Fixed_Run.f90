Program Single_Run
USE Interface_MOD
implicit none
real :: start, finish

! The model output to match with observational data: 
real, allocatable    :: Ymod(:)
integer              :: i
real                 :: mu0hat, KN, wDET,aI0, KPHY
real                 :: VTR, mz, alphaI

NAMELIST /parameters/ mu0hat,KN,  wDET, aI0,  &
    VTR, KPHY, mz, alphaI, KFe, bot_bound

!  open the namelist file and read station name.
OPEN(namlst,file='param.nml',status='old',action='read')
READ(namlst,nml=parameters)
close(namlst)

MPIRUN   = 0
taskid   = 0
numtasks = 1

call cpu_time(start) 
singlerun = .TRUE.

!Initialize the Arrays of model parameters with the biological model chosen
call SetUpArrays

! Initialize the 1D model:
call Model_setup

NDAYS    = 1800
savefile = .TRUE.
! Read model parameters:
params(imu0)      = mu0hat
params(iVTR)      = VTR
params(iKPHY)     = KPHY
params(iwDET)     = wDET
params(iaI0_C)    = aI0
params(iKN)       = KN
params(ialphaI)   = alphaI
params(imz)       = mz
params(iKFe)      = KFe

do i = 1, NPar
   write(6, 101) trim(ParamLabel(i)), params(i)
enddo
101 format(A6, 2x, F10.3)

allocate(Ymod(ANobs))
Ymod(:) = 0d0
call CalcYSSE(params, Ymod, SSqE)

do i = 1, size(SSqE)
   write(6, 102) 'SSqE(',i,')=', SSqE(i)
enddo
102 format(A5,I2,A2,1x,1pe10.3)
call cpu_time(finish) 
print '("Time = ",f8.3," mins.")', (finish-start)/60.0 
End program
