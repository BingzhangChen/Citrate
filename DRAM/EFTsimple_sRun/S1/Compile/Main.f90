PROGRAM AMAssim
USE sub_mod
implicit none
real(4) :: start,finish,t2,t1
logical :: there
! All parameters varied in an Identical Twin Test assimilation.
integer :: i,j,k,row,col
integer :: Readfile   = NO  ! Read parameter set from "enssig" and "enspar", and start from Subpcurr
real, allocatable  :: enspar1(:,:)    ! scratch matrix to store previous runs of parameters
integer            :: NR_enspar=100
integer            :: N_cvm
namelist /MCMCrun/ nruns, EnsLen, NDays, Readfile, NR_enspar, MPIRUN

 !  open the namelist file and read station name.
open(namlst,file='Model.nml',status='old',action='read')
read(namlst,nml=MCMCrun)
close(namlst)
     
if (MPIRUN .eq. 1) then
   ! ***** Initialize MPI *****
   call MPI_INIT(ierr)
   
   ! Returns the size of the group associated with a communicator.
   call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
   
   ! Determines the rank of the calling process in the communicator
   call MPI_COMM_RANK(MPI_COMM_WORLD, taskid, ierr)
   
   ! Notes: each mpi process proposes a new point based on the current position
   ! and the updated Pcvm. In each subprocess, the new Loglikelihood is calculated
   ! and compared with previous one. Acceptance is calculated and the new position
   ! is estimated. After that, the new position is sent to the root process, which
   ! combines all child processes, writes to the enspar and enssig file, calculates
   ! the new Pcvm and sends to all child processes
else
   taskid  = 0
   numtasks= 1
endif

call cpu_time(start) 

!Initialize the Arrays of model parameters with the biological model chosen
call SetUpArrays

! Initialize the 1D model:
call Model_setup

! Initialize the random number generator with a seed
! call sgrnd(17001)
call random_seed()

! An initial guess that is some factor times the inital parameter estimates 
subpguess = Npv

!  Set current to guess
subpcurr  = subpguess
subppro   = subpguess
subppro2  = subpguess
BestLogLike = -1d12  
startrun  = 0
      
! Estimate the priors based on initial parameter values
call EstimatePriors(PriorCvm, InvPriorCvm, error)

! Set the labels for the standard deviations for each type of observation 
call SetSigmaLabels
      
! Take the intial Covariance matrix to be the Prior Covariance Matrix
! (these are in compacted form)
! PriorCvm is the prior covariance matrix generated by the subroutine EstimatePriors 
! Cvm is the Covariance Matrix for parameters
Cvm = PriorCvm
Rwtold = 1d0  ! Set the initial weight equal to One
      
! Set the Proposal covariance matrix, by scaling the Parameter Covariance matrix
Pcvm = Cvm*Spcvm/Np2Vary  ! Np2Vary = d in p. 11 of Laine 2008

! Add a small term to the diagonal entries, so that the matrix will not be singular. 
do k = 1, Np2Vary
   Pcvm(k*(k+1)/2) = Pcvm(k*(k+1)/2) + CvEpsilon*Spcvm/Np2Vary
enddo

! Calculate the Cholesky factor for Pcvm, which is called Rchol here. 
! The Pcvm and Rchol will not vary until the burn-in finishes.
call cholesky(Pcvm,Np2Vary,Np2Vary*(Np2Vary+1)/2,Rchol,nullty,error)

if (taskid .eq. 0) then
   write(6, *) 'Initial Covariance matrix = '
   do row = 1, Np2Vary
      write(6,3000) (Cvm(row*(row-1)/2+col), col = 1, row )
   end do

   write(6, *) 'Initial Proposal Covariance matrix = '
   
   DO row = 1, Np2Vary
      write(6,3000) (Pcvm(row*(row-1)/2+col), col = 1, row )
   End do
   
   write(6, *) 'Cholesky factor for Proposal Covariance matrix = '
   DO row = 1, Np2Vary
      write(6,3000) ( Rchol(row*(row-1)/2+col), col = 1, row )
   Enddo
endif


! Reset subpmean to Zero for the new run
subpmean     = 0d0
subpcurrmean = 0d0
sigmamean    = 0d0

IF(nruns .gt. 1) THEN
  ! Read the ensemble files
  IF (Readfile) then
     inquire(FILE=epfn,exist=there)
     if (.not. there) then
        write(6,*) 'Cannot find the file ',epfn,'! Quit!'
        stop
     else
     ! Read the enspar file
       i = NPar+1
       allocate(enspar1(NR_enspar,i))
       enspar1(:,:) = 0d0
       call Readcsv(epfn,NR_enspar,i, enspar1)

     ! Calculate CVM (Take the final N_cvm chains, otherwise Pcvm too large)
       N_cvm = NR_enspar-BurnInt
       if (NR_enspar > N_cvm) then
          k = NR_enspar-N_cvm+1
          Rwtold=real(N_cvm)
       else
          k = 1
          Rwtold=real(NR_enspar)
       endif

       do i=k,NR_enspar
          do j=1,NPar
             Apvcurr(j)=enspar1(i,j+1)
          enddo

          subpcurr = Npv_(Apvcurr)

          call UpdateCVM(Cvm,subpcurrmean,dble(i-k+1),subpcurr,subpcurrmean,Cvm)
       enddo

    
       Pcvm = Cvm*Spcvm/NPar
     ! add a small term to the diagonal entries,
     ! so that the matrix will not be singular. 
       do k = 1, NPar
         Pcvm(k*(k+1)/2)=Pcvm(k*(k+1)/2) + CvEpsilon*Spcvm/NPar
       enddo

       if (taskid .eq. 0) then
          write(6, *) 'Updated Covariance matrix = '
          do row = 1, NPar
             write(6,3000) (Cvm(row*(row-1)/2+col), col = 1, row )
          end do
          write(6, *) 'Updated Proposal Covariance matrix = '
          DO row = 1, Np2Vary
             write(6,3000) (Pcvm(row*(row-1)/2+col), col = 1, row )
          End do
       endif

       startrun=NR_enspar

     ! Obtain current position:
       do i = 1, NPar
          Apvcurr(i) = enspar1(NR_enspar,i+2)
       enddo

       subpcurr    = Npv_(Apvcurr)
       subppro     = subpcurr
       subpbest    = subpcurr
       CurrLogLike = enspar1(NR_enspar,2)
       deallocate(enspar1)

       if (taskid .eq. 0) then
        ! Write out current parameters:
        do i = 1, Np2Vary
           write(6,101) ParamLabel(i), subpcurr(i), Apvcurr(i)
        enddo
       endif

     ! Read the enssig file (containing two stations)
       i = 1+NDTYPE*Nstn*2
       allocate(enspar1(NR_enspar,i))
       enspar1(:,:) = 0d0
       call Readcsv(esfn,NR_enspar,i, enspar1)

       do i=1,NDTYPE * Nstn
          sigma(i)=enspar1(NR_enspar,1            +i)
           SSqE(i)=enspar1(NR_enspar,1+NDTYPE*Nstn+i)
       enddo

       sigmabest= sigma
       CurrSSqE = SSqE
       deallocate(enspar1)
     endif
  ELSE 
     if (taskid .eq. 0) then
     ! Create new ensemble files
     ! Parameter Ensemble file (only one file needed)
       open(epfint,file=epfn,status='replace',action='write')
       write(epfint,1800) (ParamLabel(i), i = 1, Np2Vary)
       close(epfint)

     ! Sigma (standard error) Ensemble file (One file needed)
       open(esfint,file=esfn,status='replace',action='write')
       write(esfint,1900)  (SigmaLabel(i), i = 1,NDTYPE*Nstn), &
                           ( SSqELabel(i), i = 1,NDTYPE*Nstn)
       close(esfint)
    
     ! Output Ensemble files to store the ensemble of simulated values
     ! Be consistent with subroutine modelensout
       open(eofint, file=eofn,status='replace',action='write')
       write(eofint,'(5(A8))')  'RunNo   ',   &
                                'DOY     ',   &
                                'Depth   ',   &
                                'Name    ',   &
                                'Value   '
       close(eofint)
     endif
  ENDIF
ENDIF

! First, one simulation, to initialize all routines
! This avoids the problem that
! Parameter Values will be read from the data file on the initial run;
! if they are different than the values in this program for Assimilated Params,
! this could return a much different cost. 
! If the cost with the Values in the Parm file is much lower,
! the assimilation may almost never accept any new parameter sets)

MeanLogLike = 0d0

if (taskid .eq. 0) then
  write(6, *) ' Testing running time for each model run:'
  ! Write output 
  ! run with the Best Parameters, writing to the best output file
  call cpu_time(t1)
  open(bofint, file=bofn,action='write',status='replace')
endif

savefile = .FALSE.
call model(bofint, subppro, SSqE )
CurrLogLike = CalcLogLike(SSqE,sigma,subppro)

if (taskid .eq. 0) then
   write(6, 1001) CurrLogLike
   close(bofint)
   call cpu_time(t2)
   print '("One model run takes ",f8.3," seconds.")', t2-t1 
endif

CurrSSqE = SSqE
!------------------------------------------
!	HERE STARTS THE MAIN LOOP
!------------------------------------------
! The counter for the jth run
jrun        = startrun + 1

! Total number of runs
nruns       = nruns + startrun
BestLogLike = CurrLogLike 
call MCMC_adapt

if (taskid .eq. 0) then
   write(6, *) ' Finished the main loop for assimilation! '
   write(6, *)
   write(6, *) ' Final Proposal Covariance matrix = '
   DO row = 1, Np2Vary
      write(6,3000) ( Pcvm( row*(row-1)/2 + col ), col = 1, row )
   ENDDO
   
   write(6, *) ' Writing the last entry in',          &
               ' the ensembles of simulated values.'
   !write output of Ensemble file for simulated values
   call modelensout(eofint, jrun, subpcurr, dumE)

   write(6,*) ' Write out the simulated results with best parameters:'
   call modelensout(eofint, jrun, subpbest, dumE)

   close(eofint)

101  format(A8, 2(1x, 1pe12.2))
1001 format(/,'LogL = ',1pe13.3,/)

1800 format('LogL     ', 100(a15) )
1900 format('LogL     ', 100(a12,2x) )
3000 format(5x,<NPar>(1pe8.1,1x))
   call cpu_time(finish)
   print '("Time = ",f8.3," hours.")', (finish-start)/3600.0 
endif

!End mpi
if (MPIRUN .eq. 1) call MPI_finalize(ierr)

END PROGRAM AMAssim
