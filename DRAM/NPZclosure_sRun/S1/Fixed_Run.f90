Program Single_Run
USE Interface_MOD
implicit none
real :: start, finish

! The model output to match with observational data: 
real, allocatable    :: Ymod(:), obs(:)
integer              :: i

MPIRUN   = 0
taskid   = 0
numtasks = 1

call cpu_time(start) 
singlerun = .TRUE.

!Initialize the Arrays of model parameters with the biological model chosen
call SetUpArrays

! Initialize the 1D model:
call Model_setup

NDAYS    = 1080
savefile = .TRUE.

! Assign parameter values
params(imu0)   = log(2.5)
params(iIopt)  = log(1000.)
params(iaI0_C) = log(0.05)
params(iKN)    = log(0.2 )
params(iDp)    = log(0.1 )
params(iwDET)  = log(0.01)
params(ibeta)  = log(0.0001)
params(iKPHY)  = log(0.5)
params(igmax)  = log(2.0)
params(imz)    = log(0.05)
params(iVNO30) = log(0.1)
params(iVPHY0) = log(0.1)

do i = 1, NPar
   write(6, 101) trim(ParamLabel(i)), exp(params(i))
enddo
101 format(A6, 2x, F10.3)

allocate(Ymod(ANobs))
allocate(obs(ANobs))
Ymod(:) = 0d0; obs(:)=0d0
call CalcYSSE(params, Ymod, obs, SSqE)

do i = 1, size(SSqE)
   write(6, 102) 'SSqE(',i,')=', SSqE(i)
enddo
102 format(A5,I2,A2,1x,1pe10.3)
call cpu_time(finish) 
print '("Time = ",f8.3," mins.")', (finish-start)/60.0 
End program
