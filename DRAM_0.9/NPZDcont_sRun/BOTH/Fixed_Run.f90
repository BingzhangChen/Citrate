Program Single_Run
use sub_mod
implicit none
real :: start, finish

! The model output to match with observational data: 
real, allocatable    :: Ymod(:)
integer              :: i
real                 :: mu0hat, KN, wDET, Q0N,aI0, KPHY,alphaG
real                 :: alphaKN, VTR, gmax, mz, alphaI,dustso,gb

NAMELIST /parameters/ mu0hat,KN,  wDET, Q0N,aI0,alphaKN, alphamu, &
    VTR, gmax,KPHY, mz, alphaI,alphaG, betamu, KFe, alphaFe, dustso,gb, bot_bound

!  open the namelist file and read station name.
OPEN(namlst,file='param.nml',status='old',action='read')
READ(namlst,nml=parameters)
close(namlst)

call cpu_time(start) 
singlerun = .TRUE.

!Initialize the Arrays of model parameters with the biological model chosen
call SetUpArrays

! Initialize the 1D model:
call Model_setup

NDAYS    = 1440
savefile = .TRUE.
! Read model parameters:
params(imu0)      = mu0hat
params(iVTR)      = VTR
params(iKPHY)     = KPHY
params(iwDET)     = wDET
params(iaI0_C)    = aI0
params(iQ0N)      = Q0N
params(iKN)       = KN
params(ialphaKN)  = alphaKN
params(ialphamu)  = alphamu
params(ialphaI)   = alphaI
params(igmax)     = gmax
params(imz)       = mz
params(ibetamu)   = betamu
params(iKFe)      = KFe
params(ialphaFe)  = alphaFe
params(idustsol)  = dustso
params(ialphaG)   = alphaG
params(igb)       = gb

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
