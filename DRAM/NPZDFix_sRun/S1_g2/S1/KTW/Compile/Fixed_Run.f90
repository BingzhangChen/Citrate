Program Single_Run
use sub_mod
implicit none
real :: start, finish

! The model output to match with observational data: 
real, allocatable    :: Ymod(:)
integer              :: i

call cpu_time(start) 
singlerun = .TRUE.

!Initialize the Arrays of model parameters with the biological model chosen
call SetUpArrays

! Initialize the 1D model:
call Model_setup

NDays    = 1440
savefile = .TRUE.
do i = 1, NPar
   write(6, 101) trim(ParamLabel(i)), params(i)
enddo
101 format(A6, 2x, F10.3)

allocate(Ymod(ANobs))
Ymod(:) = 0d0
call CalcYSSE(params, Ymod, SSqE)

write(6, *) 'SSqE: ', SSqE

call cpu_time(finish) 
print '("Time = ",f8.3," mins.")', (finish-start)/60.0 

End program
