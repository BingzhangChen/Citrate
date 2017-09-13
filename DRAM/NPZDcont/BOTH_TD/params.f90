Module param_MOD
use MPI
implicit none

! Number of stations to run
integer, parameter :: Nstn       = 2

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
integer, parameter :: TNFo       = ewstr ! TOtal number of forcings

! Total of observation times in forcing data
character(LEN=5), parameter :: LabelForc(TNFo) &
     = (/'temp ','NO3  ','Aks  ','wROMS','par  ', 'solfe','fer  ', 'PO4  ', 'wstr '/)

integer,          parameter :: NFobs(TNFo)     & 
     = (/    12,   12, 360,     12,  12,      12,   12,    12, 12 /)

integer, parameter :: N_fer=33
end module
