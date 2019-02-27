Program Single_Run
use sub_mod
implicit none
real :: start, finish

! The model output to match with observational data: 
real, allocatable    :: Ymod(:)
integer              :: i
real                 :: mu0hat, KN, A0N, wDET, aI0_C,Q0N,aI0
real                 :: alphaKN, alphaA

namelist /parameters/    mu0hat, KN, wDET, aI0_C,Q0N,aI0

!  open the namelist file and read station name.
open(namlst,file='param.nml',status='old',action='read')
read(namlst,nml=parameters)
close(namlst)

call cpu_time(start) 
singlerun = .TRUE.

!Initialize the Arrays of model parameters with the biological model chosen
call SetUpArrays

! Initialize the 1D model:
call Model_setup

NDays    = 1080
savefile = .TRUE.
! Read model parameters:
params(imu0)=mu0hat
params(iwDET)=wDET
params(iaI0_C)=aI0_C
params(iQ0N)=Q0N
params(iaI0)=aI0
params(iKN) =KN

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
102 format(A5,I1,A2,1x,1pe10.3)
call cpu_time(finish) 
print '("Time = ",f8.3," mins.")', (finish-start)/60.0 

End program
