MODULE MOD_1D
use BIO_MOD
implicit none
public
! Grid parameters
real, private, parameter :: hmax   = 2.5d2   ! Total water depth
real, private, parameter :: thetaS = 2d0    ! surface stretching parameter
real, private, parameter :: dtsec  = 6D2    ! time step in seconds
real, private, parameter ::d_per_s = 864d2   ! how many seconds in one day
real, private, parameter :: Yup    = 0d0
real, private, parameter :: Ydw    = 0d0 
                      !how many seconds of one year
integer, private, parameter :: y_per_s = INT(d_per_s*360), &

!  Number of vertical points of NO3,Temp, and Aks profile
                      N_NO3   = 37,               &  
                      N_Temp  = 57,               &
                      N_w     = 40,               &
                      N_par   = 1 ,               &
                      N_Dust  = 1,                &
                      N_fer   = 33,               &
                      nsave   = INT(d_per_s)/INT(dtsec) ! Timesteps to save

integer, private            :: N_Aks(Nstn)              ! Station dependent

! Use Taketo's Aks data for K2 and my own ROMS data for S1 

! Forcing data time indices 
real, private, target :: obs_time_temp(NFobs(etemp), Nstn)
real, private, target :: obs_time_NO3( NFobs(eNO3) , Nstn)
real, private, target :: obs_time_Aks( NFobs(eAks) , Nstn)
real, private, target :: obs_time_w(   NFobs(ew  ) , Nstn)
real, private, target :: obs_time_par( NFobs(ePAR) , Nstn)
real, private, target :: obs_time_Dust(NFobs(eDust), Nstn)
real, private, target :: obs_time_Fer( NFobs(eFer) , Nstn)

! Forcing data
real, private :: obs_NO3( N_NO3 ,  1+NFobs(eNO3 ))
real, private :: obs_fer( N_fer ,  1+NFobs(efer ))
real, private :: obs_Temp(N_Temp,  1+NFobs(etemp))
real, private :: obs_PAR( N_par ,  1+NFobs(ePAR ))
real, private :: obs_w(   N_w   ,  1+NFobs(ew   ))
real, private :: obs_Dust(N_Dust,  1+NFobs(eDust))
real, private, allocatable :: obs_Aks(:,:) ! (N_Aks ,  1+NFobs(eAks ))

real, pointer :: pb(:), pc(:,:)
integer       :: nn  ! a scratch integer for the dimensions of matrix

! Vertically interpolated temperature and Aks at each obs timing
real, private, target :: VTemp(nlev, NFobs(etemp), Nstn)
real, private, target :: VAks(0:nlev,NFobs(eAks ), Nstn)
real, private, target :: Vw(  0:nlev,NFobs(ew )  , Nstn)
real, private, target :: VPAR(1,     NFobs(ePAR) , Nstn)
real, private, target :: VDust(1,    NFobs(eDust), Nstn)

logical, public  :: savefile
logical, public  :: INCLUDESIZE = .FALSE.

! New calculated state variables
real, private, allocatable    :: ww(:,:)  
!!$$-----------------------------------------------------------------
!$ The declaration of the data part:
! length of string for labels 
integer, parameter   :: LabelLen = 7 

! NDTYPE is the number of types of data (or observations), 
integer              :: NDTYPE
character(LabelLen), allocatable :: DataLabel(:)

! Number of observations (data points) of each data type
integer, allocatable :: NDPTS(:,:)
integer              :: TNobs(Nstn), ANobs

! Indeces for data type:
integer, parameter   :: itNO3 = 1, itCHL = 2, itNPP = 3
integer, parameter   :: itPON = 4
integer, parameter   :: itP10 = 5, itP03 = 6, itP01 = 7, itP_1 = 8

! The number of days for model run, needs to change at different stages during MCMC
integer              :: NDays = 1080

! The data arrays are read from csv files. 
! Dimensions of the input data
integer, allocatable :: ncol(:,:)
integer, allocatable :: nrow(:,:)
integer, allocatable :: Tnrow(:)

! Declare the dimensions of the data:
real,    allocatable ::  TINData(:,:)
real,    allocatable ::  CHLData(:,:)
real,    allocatable ::  NPPData(:,:)
real,    allocatable ::  PONData(:,:)
real,    allocatable ::  C2NData(:,:)
real,    allocatable :: SizeData(:,:)
real,    allocatable ::  OBSData(:,:)

! DOY and Depth for Observational data assembled as a single matrix 
real,    allocatable ::  OBS_DOY(:), OBS_Depth(:)
! Data label for OBS data
character(LabelLen), allocatable :: OBS_Label(:)

! Initial profile of NO3:
real, private               :: NO3(nlev,1,Nstn)

! Initial profile of fer:
real, private               :: fer(nlev,1,Nstn)

CONTAINS
!-----------------------------------------------------------------------
! Call other subroutines to read and set-up data, as needed 
! This subroutine is called before the real model run
subroutine Setup_OBSdata
implicit none
character(LEN=20)    :: TIN_OBS_file(Nstn)
character(LEN=20)    :: CHL_OBS_file(Nstn)
character(LEN=20)    :: NPP_OBS_file(Nstn)
character(LEN=20)    :: PON_OBS_file(Nstn)
character(LEN=20)    :: C2N_OBS_file(Nstn)
character(LEN=20)    :: SIZE_OBS_file(Nstn)
real,    allocatable :: DOY(:), Depth(:)
integer              :: k, i,oi,j
integer              :: NN  ! Number of data types in observational data
integer              :: NL  ! for counting of OBS_DOY et al.
integer, allocatable :: kk(:)

Select case(Model_ID)
  case(NPZDFix,NPPZDD,EFTPPDD,Geidersimple, EFTsimple, NPZDFixIRON, GeidsimIRON,EFT2sp,NPZD2sp, EFTsimIRON)
     ! Data types must be the same for different stations.
     ! But can be different for different models
     NDTYPE = 4  !TIN, CHL, PP, PON
     INCLUDESIZE = .FALSE.
  case(NPZDdisc, EFTdisc,EFTcont, Geiderdisc,NPZDCONT)
     NDTYPE = 8  !TIN, CHL, PP, PON, CHL>10, CHL3-10,CHL1-3,CHL<1
     INCLUDESIZE = .TRUE.
     if (trim(Stn(1))  .eq. 'HOT' .and. Nstn .eq. 1) then
          NDTYPE = 4
     INCLUDESIZE = .FALSE.
     endif
  case default
     write(6,*) 'Model option incorrect! Quit!'
     stop
End select

allocate(DataLabel(NDTYPE),  STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating DataLabel ***"

DataLabel(itNO3) = 'TIN'
DataLabel(itCHL) = 'CHL'
DataLabel(itNPP) = 'NPP'
DataLabel(itPON) = 'PON'

If (INCLUDESIZE) then
  DataLabel(itP10) = 'P10'
  DataLabel(itP03) = 'P03'
  DataLabel(itP01) = 'P01'
  DataLabel(itP_1) = 'P_1'
  allocate(nrow(5,Nstn), STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Error in allocating nrow ***"
  allocate(ncol(5,Nstn), STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Error in allocating ncol ***"
Else
  allocate(nrow(NDTYPE,Nstn), STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Error in allocating nrow ***"
  allocate(ncol(NDTYPE,Nstn), STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Error in allocating ncol ***"
Endif

allocate(NDPTS(NDTYPE,Nstn),  STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Error in allocating NDPTS ***"

Do j = 1, Nstn
  write(6,*) 'Stn name: ',Stn(j)
! Assign the data dimension, must be consistent with external file:
  if (trim(Stn(j)) .eq. 'S1') then
      N_Aks(j)       = 40
      NDPTS(itNO3,j) = 902
      NDPTS(itCHL,j) = 426
      NDPTS(itNPP,j) = 128
      NDPTS(itPON,j) = 32

      if (INCLUDESIZE) then
         NDPTS(itP10,j) = 166
         NDPTS(itP03,j) = NDPTS(itP10,j)
         NDPTS(itP01,j) = NDPTS(itP10,j)
         NDPTS(itP_1,j) = NDPTS(itP10,j)
      endif
  else if (trim(Stn(j)) .eq. 'K2') then
      N_Aks(j)       = 25
      NDPTS(itNO3,j) = 974
      NDPTS(itCHL,j) = 470
      NDPTS(itNPP,j) = 112
      NDPTS(itPON,j) = 29
      if (INCLUDESIZE) then
         NDPTS(itP10,j) = 143
         NDPTS(itP03,j) = NDPTS(itP10,j)
         NDPTS(itP01,j) = NDPTS(itP10,j)
         NDPTS(itP_1,j) = NDPTS(itP10,j)
      endif
  else if (trim(Stn(j)) .eq. 'HOT') then
      N_Aks(j)       = 40
      NDPTS(itNO3,j) = 3910
      NDPTS(itCHL,j) = 8180
      NDPTS(itNPP,j) = 1659
      NDPTS(itPON,j) = 2251
  else
      write(6,*) 'Station number incorrect! Stop!'
      stop
  endif
  TNobs(j)  = sum(NDPTS(:,j))
  ANobs     = sum(TNobs(:))   !All data points of all stations
  ! Assign obs. data matrix: 
  IF (INCLUDESIZE) THEN
     do i = 1,5
        nrow(i,j)=NDPTS(i,j)
        if (i < 5) then
           ncol(i,j)=3
        else
           ncol(i,j)=6
        endif
     enddo
  ELSE
     do i = 1,NDTYPE
        nrow(i,j)=NDPTS(i,j)
        ncol(i,j)=3
     enddo
  ENDIF
  write(6,'(A31,1x,A3,1x,A1,1x,I5)') 'Total number of observations of',Stn(j),'=',TNobs(j)
Enddo  ! ==> End of assigning data numbers at different stations
  
! Total number of observations of each data type for all stations
NN = size(nrow,1)
allocate(Tnrow(NN))
Do i = 1, NN
   Tnrow(i)=sum(nrow(i,:))
Enddo
!Start to assign TIN,Chl, and NPP data:
allocate( TINData(Tnrow(1),ncol(1,1)), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating TINData ***"

allocate( CHLData(Tnrow(2),ncol(2,1)), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating CHLData ***"

allocate( NPPData(Tnrow(3),ncol(3,1)), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating NPPData ***"

allocate( PONData(Tnrow(4),ncol(4,1)), STAT = AllocateStatus)

allocate( OBSData(ANobs  ,3 ), STAT = AllocateStatus)
IF (AllocateStatus /= 0) STOP "*** Problem in allocating OBSData ***"
! Initialize observational data:
TINData(:,:) = 0d0
CHLData(:,:) = 0d0
NPPData(:,:) = 0d0
PONData(:,:) = 0d0
OBSData(:,:) = 0d0

! Initialize size data:
if (INCLUDESIZE) then
   allocate(SizeData(Tnrow(5),ncol(5,1)), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Problem in allocating SizeData ***"
   SizeData(:,:) = 0d0
endif

! Initialize observed DOY and Depth
! This is for writing out output
allocate(OBS_DOY(  ANobs))
allocate(OBS_Depth(ANobs))
allocate(OBS_Label(ANobs))

OBS_DOY(:)   = 0d0
OBS_Depth(:) = 0d0
OBS_Label(:) = ' NA    '

allocate(kk(NN))
kk(:) = 0
NL    = 0
Do j = 1, Nstn
  ! Read observational data:
  TIN_OBS_file(j) = trim(Stn(j))//'_TIN.dat' 
  CHL_OBS_file(j) = trim(Stn(j))//'_CHL.dat'
  NPP_OBS_file(j) = trim(Stn(j))//'_NPP.dat'
  PON_OBS_file(j) = trim(Stn(j))//'_PON.dat'
  
  call Readcsv(TIN_OBS_file(j),nrow(1,j),ncol(1,j), TINData((kk(1)+1):(kk(1)+nrow(1,j)),:))
  call Readcsv(CHL_OBS_file(j),nrow(2,j),ncol(2,j), CHLData((kk(2)+1):(kk(2)+nrow(2,j)),:))
  call Readcsv(NPP_OBS_file(j),nrow(3,j),ncol(3,j), NPPData((kk(3)+1):(kk(3)+nrow(3,j)),:))
  call Readcsv(PON_OBS_file(j),nrow(4,j),ncol(4,j), PONData((kk(4)+1):(kk(4)+nrow(4,j)),:))

  if (INCLUDESIZE) then
     ! Use percentage data
     SIZE_OBS_file(j) = trim(Stn(j))//'_size_Perc.dat'
     call Readcsv(SIZE_OBS_file(j),nrow(5,j),ncol(5,j),SizeData((kk(5)+1):(kk(5)+nrow(5,j)),:))
  endif

  ! Give DOY and Depth for all the obs. data
  k = 0
  do i = 1, NDTYPE
    allocate(  DOY(NDPTS(i,j)))
    allocate(Depth(NDPTS(i,j)))

    selectcase(i)
    case(itNO3)   ! Nitrate
      DOY   = TINData((kk(1)+1):(kk(1)+nrow(1,j)),1)
      Depth = TINData((kk(1)+1):(kk(1)+nrow(1,j)),2)
    case(itCHL)   ! CHL
      DOY   = CHLData((kk(2)+1):(kk(2)+nrow(2,j)),1)
      Depth = CHLData((kk(2)+1):(kk(2)+nrow(2,j)),2)
    case(itNPP)   ! NPP
      DOY   = NPPData((kk(3)+1):(kk(3)+nrow(3,j)),1)
      Depth = NPPData((kk(3)+1):(kk(3)+nrow(3,j)),2)
    case(itPON)   ! PON
      DOY   = PONData((kk(4)+1):(kk(4)+nrow(4,j)),1)
      Depth = PONData((kk(4)+1):(kk(4)+nrow(4,j)),2)
    case(itP10,itP03,itP01,itP_1)   ! Size-fractionated Chl
      DOY   = SizeData((kk(5)+1):(kk(5)+nrow(5,j)),1)
      Depth = SizeData((kk(5)+1):(kk(5)+nrow(5,j)),2)
    case default
      print *, 'Errors in selecting data types! Quit!'
      stop
    endselect

    do oi = 1, NDPTS(i,j)
         OBS_DOY(NL+k+oi) =   DOY(oi)
       OBS_Depth(NL+k+oi) = Depth(oi)
       OBS_Label(NL+k+oi) = DataLabel(i)
    enddo
    k = k + NDPTS(i,j)
    deallocate(DOY)
    deallocate(Depth)
  enddo

  !Update kk (the index in TINData etc.)
  do oi = 1, NN
     kk(oi) = kk(oi) + nrow(oi,j)
  enddo
  !Update NL (for OBS_DOY et al.)
  NL = NL + k
  write(6,'(A27,1x,A3)') 'Observation data set up for', Stn(j)
Enddo
end subroutine Setup_OBSdata
!========================================================
subroutine Model_setup
! Called at the beginning of the program, only once
implicit none
integer            :: i,j,k
! Input files
character(LEN=20)  :: forcfile(TNFo), forcfiletime(TNFo)

!the fraction of a time step in one day
dtdays = dtsec/d_per_s
write(6,'(A13,1x,F12.3,A12)') 'Timestepping: ',dtdays,'of one day.'

! Setup grid
call setup_grid

do j = 1, Nstn
   ! Forcing files:
   do i = 1, TNFo
      forcfile(i)     = trim(Stn(j))//'_'//trim(LabelForc(i))//'.dat'
      forcfiletime(i) = trim(Stn(j))//'_'//trim(LabelForc(i))//'_time.dat'
   enddo
   
   ! Use Taketo's Aks data only for K2:
   if (trim(Stn(j)) .eq. 'K2') then
      forcfile(eAks)     = trim(Stn(j))//'_'//trim(LabelForc(eAks))//'T.dat'
      forcfiletime(eAks) = trim(Stn(j))//'_'//trim(LabelForc(eAks))//'_timeT.dat'
   endif
   
   ! Read NO3 data:
   call Readcsv(forcfile(eNO3), N_NO3, size(obs_NO3,2), obs_NO3) 
   
   ! Interpolate initial NO3:
   call gridinterpol(N_NO3,1,obs_NO3(:,1),obs_NO3(:,2),                 &
                     nlev, Z_r, NO3(:,1,j)) 
   if (do_IRON) then
     ! Read fer data:
     call Readcsv(forcfile(efer), N_fer, size(obs_fer,2), obs_fer) 
     
     ! Interpolate initial fer:
     call gridinterpol(N_fer,1,obs_fer(:,1),obs_fer(:,2),                 &
                       nlev, Z_r, fer(:,1,j)) 
   endif
  
   ! Calculate obs. time indices in seconds:
   ! Read obs. time file:
   call Readcsv(forcfiletime(ePAR), 1,NFobs(ePAR), obs_time_par(:,j)) 
   call Readcsv(forcfiletime(eNO3), 1,NFobs(eNO3), obs_time_NO3(:,j)) 
   call Readcsv(forcfiletime(etemp),1,NFobs(etemp),obs_time_temp(:,j)) 
   call Readcsv(forcfiletime(eAks), 1,NFobs(eAks), obs_time_Aks(:,j)) 
   call Readcsv(forcfiletime(ew  ), 1,NFobs(ew)  , obs_time_w(:,j)) 
   
   if (do_IRON) then  ! Read dust deposition time data:
      call Readcsv(forcfiletime(eDust), 1,NFobs(eDust), obs_time_Dust(:,j)) 
      obs_time_Dust(:,j)  = obs_time_Dust(:,j) *3d1*dble(d_per_s)  ! Dust deposition
   endif

   ! Convert obs_time to the unit of seconds:
   obs_time_par(:,j)  = obs_time_par(:,j) *3d1*dble(d_per_s)
   obs_time_w(:,j)    = obs_time_w(:,j)   *3d1*dble(d_per_s)
   obs_time_Aks(:,j)  = obs_time_Aks(:,j) *3d1*dble(d_per_s)
   obs_time_temp(:,j) = obs_time_temp(:,j)*3d1*dble(d_per_s)
   obs_time_NO3(:,j)  = obs_time_NO3(:,j) *3d1*dble(d_per_s)
   ! Read external w data:
   !call Readcsv(forcfile(ew), size(obs_w,1), size(obs_w,2), obs_w) 
   
   ! Interpolate external w data:
   
   !subroutine gridinterpol(N,cols,obs_z,obs_prof,nlev_,model_z,model_prof)
   !call gridinterpol(N_w,NFobs(ew),obs_w(:,1),       &
   !      obs_w(:,2:(NFobs(ew)+1)),                   &
   !      nlev+1, Z_w, Vw(:,:)) 
   
   ! Read external PAR data:
   call Readcsv(forcfile(ePAR), size(obs_PAR,1),size(obs_PAR,2),obs_PAR) 

   ! No need to vertically interpolate PAR data:
   do k = 1, (size(obs_PAR,2)-1)
      VPAR(1,k,j) = obs_PAR(1,k+1)
   enddo

   If (do_IRON) then
      ! Read external Dust data:
      call Readcsv(forcfile(eDust), size(obs_Dust,1),size(obs_Dust,2),obs_Dust) 

      ! No need to vertically interpolate Dust data:
      do k = 1, (size(obs_Dust,2)-1)
         VDust(1,k,j) = obs_Dust(1,k+1)
      enddo
   Endif

   ! Read external Temp data:
   call Readcsv(forcfile(etemp),size(obs_Temp,1),size(obs_Temp,2),obs_Temp) 
   
   ! Interpolate external Temp data:
   call gridinterpol(N_Temp,NFobs(etemp),obs_Temp(:,1),      &
         obs_Temp(:,2:(NFobs(etemp)+1)),                     &
         nlev, Z_r, VTemp(:,:,j)) 

   ! Read external Aks data:
   allocate(obs_Aks(N_Aks(j), 1+NFobs(eAks)),  STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Problem in allocating obs_Aks ***"
   call Readcsv(forcfile(eAks), N_Aks(j), size(obs_Aks,2), obs_Aks) 
   
   ! Interpolate external Aks data:
   call gridinterpol(N_Aks(j),NFobs(eAks),obs_Aks(:,1), &
         obs_Aks(:,2:(NFobs(eAks)+1)),               &
         nlev+1, Z_w, VAks(:,:,j)) 
   
   savefile = .FALSE.
   write(6,*) 'Model setup finished for ', Stn(j)
   deallocate(obs_Aks)
Enddo
End subroutine Model_setup
!========================================================
! This subroutine does the main work in the 1D model
! and give the output to match with the observational data
! Must be called after initialization
SUBROUTINE Timestep(TINout, CHLout, NPPout, PONout, Sizeout)
implicit none
real,    parameter  :: cnpar = 0.6
! Local scratch variables
integer  :: it,i,k,nm, j,jj, Nstep, current_day, current_DOY,DOY
integer  :: i_,j_

! Counting of the index of each data type
real     :: cff,current_sec
real     :: I_0(1), dust0(1),Aks(0:nlev),w(0:nlev)
real     :: depth(1)
real     :: CHL_(nlev),Vars1(nlev),Vars2(nlev),ww_(0:nlev)
real     :: Dust_solubility = 0.02  !Aumont et al. (2003)
real, allocatable  :: a(:,:)
integer, parameter :: mode0 = 0
integer, parameter :: mode1 = 1
real,    parameter :: Aks_th= 1D-3 ! Aks threshold for calculating average PAR in MLD 
real,    parameter :: Dust_ironfrac   = 0.035 !TOM10Appendix p. 24, unit: g/g
real,    parameter :: mon2sec         = 2592D3

character(LEN=20)  :: outfile
! The model output (final year) to match with observational data: 
!real,    intent(out)           :: dCChl
real,    intent(out)           ::  TINout(size( TINData,1),1)
real,    intent(out)           ::  CHLout(size( CHLData,1),1)
real,    intent(out)           ::  NPPout(size( NPPData,1),1)
real,    intent(out)           ::  PONout(size( PONData,1),1)
real,    intent(out), optional :: Sizeout(size(SizeData,1),4)
logical                        :: MLD_found

! For showing whether some unconstrained outputs are realistic or not
!dCChl       = 0d0    ! Difference between C:Chl and realistic boundaries
TINout(:,:) = 0d0
CHLout(:,:) = 0d0
NPPout(:,:) = 0d0
PONout(:,:) = 0d0

if (INCLUDESIZE) Sizeout(:,:) = 0d0

DO jj = 1, Nstn
  ! Initialize initial NO3:
  Vars(iNO3,:) = NO3(:,1,jj)


  ! Initialize other variables:
  do k = 1,nlev
     do i = 1,NPHY
        Vars(iPHY(i), k) = 0.1/float(NPHY)
     enddo
     Vars(iZOO,k) = 0.1
     Vars(iDET,k) = 0.1
  
     if (NVAR > iDET) then
        do i = (iDET+1), NVAR
           Vars(i,k) = 1D-2
        enddo
     endif
  enddo

  if (do_IRON) then
    ! Initialize initial fer:
     Vars(ifer,:) = fer(:,1,jj)
  endif
 
  ! Initialize output data:
  Varout(:,:) = 1D-3
  do i = 1,NVAR
     do k = 1,nlev
        Varout(i,k)= Vars(i,k)
     enddo
  enddo

  ! Sinking rate
  allocate(ww(0:nlev,NVsinkterms),  STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Problem in allocating ww ***"
  ! Initialize sinking rate (UNIT: m/s !):
  ww(:,:) = 0d0
  nn      = size(ww,2)
  do k = 1,nlev-1
   !Phytoplankton no sinking
   !Detritus sinking rate (convert to UNIT: m/s)
    if (Model_ID==NPPZDD .OR. Model_ID==EFTPPDD) then
       ww(k,nn-1)=-Params(iwDET )/dble(d_per_s) 
       ww(k,nn)  =-Params(iwDET2)/dble(d_per_s) 
    else
       ww(k,nn) = -Params(iwDET)/dble(d_per_s) 
    endif
  enddo

  ! Create output files:
  if (savefile) then

    ! Merge all data into one output file (Oct. 15 2016)
     outfile = trim(Stn(jj))//'.out'

    ! Create data out file
    open (unit=10, file = outfile, status = 'replace')
    write(10, 1) 'Type','Timestep','Days',Z_r
  endif
  
  1 format(2(A10,1x), A7, <nlev  >(2x,F12.5))
  3 format(A10,1x,A7, <nlev+1>(2x,F12.5))
  
  ! Number of time steps
  Nstep = NDays*INT(d_per_s)/INT(dtsec) 
  ! 'START TIME STEPPING'
  DO it = 1, Nstep+1
  
  ! Calculate current timing (zero is starting time):
     current_sec = float(it-1)*dtsec

!%%%%%%%%%%%%%%%%%%%%%%%%For each time step, read in external environmental data%%%%%%%%%%%%%%%
    ! Interpolate Temp data throughout the water column:
     pb => obs_time_temp(:,jj)
     nn = size(pb,1)
     pc => VTemp(:,:,jj)
     call time_interp(int(current_sec),nn, nlev,pb,pc,Temp)

    ! Interpolate temporal PAR data: 
     pb => obs_time_par(:,jj)
     nn = size(pb,1)
     pc => VPAR(:,:,jj)
     call time_interp(int(current_sec),nn,1,pb,pc,I_0)
    
    ! Convert the unit of par to W m-2:
     cff = I_0(1)/4d-1
  
    ! Calculate light field:
     CHL_(:) = Varout(oCHLt,:)
     call Calculate_PAR(cff, nlev, Hz, CHL_, PAR)
      
     if (do_IRON) then
        ! Interpolate temporal Dust data: 
        pb => obs_time_Dust(:,jj)
        nn = size(pb,1)
        pc => VDust(:,:,jj)
        call time_interp(int(current_sec),nn,1,pb,pc, dust0)
  
     ! Iron atmospheric deposition:
     ! Dust deposition unit: kg/m2/month. Assume 0.035g Fe per g of dust
     ! so need to convert into nM at each time step
     ! Deposition = Dust*10^12*solubility*Fe%/56*dt/mon2sec*surface_area/surface_grid_volume/1000
         
       !cff= Dust_solubility*dust0(1)*1D9*Dust_ironfrac/55.85/Hz(nlev)/mon2sec*dtsec                     
       Dust_solubility=params(idustsol)
       cff= Dust_solubility*dust0(1)*1D9/55.85/Hz(nlev)/mon2sec*dtsec                     
       
       !write(6,*) 'Annual added dissolved Fe (umol/m2) = ', cff * Hz(nlev)*mon2sec/dtsec*12.
       
       Vars(ifer,nlev) = Vars(ifer,nlev) + cff
     endif

    ! Interpolate Aks data throughout the water column:
     pb => obs_time_Aks(:,jj)
     nn = size(pb,1)
     pc => VAks(:,:,jj)
     call time_interp(int(current_sec), nn, nlev+1,pb, pc, Aks)
     nullify(pb)
     nullify(pc) 

    ! Calculate the vertical grid index (N_MLD) at the bottom of MLD:
     MLD_found = .FALSE.
     do k = nlev, 1, -1
        if (Aks(k) .gt. Aks_th .and. Aks(k-1) .le. Aks_th) then
           N_MLD     = k
           MLD_found = .TRUE.
           exit
        endif
     enddo
     if (.not. MLD_found) then
        if (Aks(nlev) .lt. Aks_th) then
           ! Even surface layer is not mixed enough
           N_MLD=nlev + 1
        else
           N_MLD=1  ! Mixed throughout the whole water column
        endif
     endif
    ! Calculate average PAR within the surface mixed layer (from nlev to N_MLD):
     if (N_MLD .le. nlev) then
        PARavg=0d0
        do k=nlev,N_MLD,-1
           PARavg=PARavg+PAR(k)*Hz(k)
        enddo 
        PARavg = PARavg/abs(Z_w(N_MLD-1))
     endif
  
    ! Interpolate w data throughout the water column:
    ! call time_interp(int(current_sec), size(obs_time_w,1), nlev+1,&
    !      obs_time_w, Vw, w)
    ! w(0   ) = 0d0
    ! w(nlev) = 0d0
    IF (mod(it, nsave) .EQ. 1) THEN
      ! Check whether the values are valid:
      do j = 1,NVAR
         do k=1,nlev
            if( (Vars(j,k) .ne. Vars(j,k)) .OR. (Vars(j,k) .le. 0.) ) then
                write(6,*) 'j = ',j
                write(6,*) 'k = ',k
                write(6,*) 'At time step ',it
                write(6,*) 'The variable ',trim(Labelout(j+ow)), &
                          ' is invalid at depth ',Z_r(k)
                write(6,*) 'Vars(j,k) =',  Vars(j,k)
                stop
            endif
         enddo 
      enddo
  
     ! Calculate model time in days:
       current_day = int(current_sec/d_per_s)
  
    ! Calculate DATE OF the YEAR (DOY)
       current_DOY = mod(current_day, 360)
  
      ! Save data to output files:
      if (savefile) then

      ! Save Temp data into the output file:
         write(10, 2) Labelout(oTemp), it, current_day, Temp

      ! Save PAR data into the output file:
         write(10, 2) Labelout(oPAR),  it, current_day, PAR
  
      ! Save Aks data into the output file:
         write(10, 2) Labelout(oAks),  it, current_day, Aks(1:nlev)

      ! Save w data into the output file:
      !   write(9+ow  , 4) it, current_day, w
  
      ! Save state variables and diagnostics:
         do i=(ow+1),(Nout+ow)
            Vars1(:) = Varout(i-ow,:)
            write(10, 2) Labelout(i),  it, current_day, Vars1
         enddo
      endif  !==> End of saving results
    
  2 format(A10,1x, I10,1x,I7, <nlev  >(2x,1pe12.3))
  !-------------------------------------------------------------------------
    ! Calculate model outputs (final year) to match with obs. data
   If ((NDays-current_day) .le. 360) Then
 
         ! Loop through the DOY of the observed data
     Do i=1, sum(nrow(:,jj))
        ! Calculate the index for the data
        if (jj .eq. 1) then
           nm = i
        else
           nm = 0
           do j = 1, (jj-1)
              nm = nm + sum(nrow(:,j))
           enddo
           nm = nm + i
        endif
        DOY      = INT(min(OBS_DOY(nm),360.0))
        depth(1) = -abs(OBS_Depth(nm))  !Consistent with ROMS convention
  
        IF (DOY .eq. current_DOY) then

         ! Check whether C:Chl (may add other parameters) realistic or not:
         ! allocate(a(nlev,NPHY))
         ! a(:,:) = 0d0
         ! do i_ = 1, nlev
         !   do j_ = 1, NPHY
         !    a(i_,j_) = 12d0/Varout(oTheta(j_),i_)   !unit: gC:gChl
         !     dCChl=dCChl+max((a(i_,j_)-8D2),0d0)+abs(min((a(i_,j_)-5D0),0d0))
         !   enddo
         ! enddo
         ! deallocate(a)

          allocate(a(nlev,1))
          a(:,:)=0d0

          if (i .le. nrow(1,jj)) then
             if (jj .eq. 1) then
                nm = i
             else
                nm = 0
                do j = 1, (jj-1)
                   nm = nm + nrow(1, j)
                enddo
                nm = nm + i
             endif

           ! Calculate TIN output:
           a(:,1) = Vars(iNO3,:)
           call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                             TINout(nm,1))       
          elseif (i .le. (nrow(1,jj)+nrow(2,jj)) ) then
             if (jj .eq. 1) then
                nm = i- nrow(1,jj)
             else
                nm = 0
                do j = 1, (jj-1)
                   nm = nm + nrow(2, j)
                enddo
                nm = nm + i - nrow(1,jj)
             endif
          ! Calculate CHL output:
           a(:,1) = Varout(oCHLt,:)
           call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                CHLout(nm,1))   
  
          elseif (i .le. (nrow(1,jj)+nrow(2,jj)+nrow(3,jj)) ) then
             if (jj .eq. 1) then
                nm = i- nrow(1,jj) - nrow(2,jj)
             else
                nm = 0
                do j = 1, (jj-1)
                   nm = nm + nrow(3, j)
                enddo
                nm = nm + i - nrow(1,jj) - nrow(2,jj)
             endif

          ! Calculate NPP output:
           a(:,1) = Varout(oPPt,:)   ! Carbon based PP
           call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                NPPout(nm,1))  

          elseif (i .le. (nrow(1,jj)+nrow(2,jj)+nrow(3,jj)+nrow(4,jj))) then
             if (jj .eq. 1) then
                nm = i- nrow(1,jj) - nrow(2,jj) - nrow(3,jj)
             else
                nm = 0
                do j = 1, (jj-1)
                   nm = nm + nrow(4, j)
                enddo
                nm = nm + i - nrow(1,jj) - nrow(2,jj) - nrow(3,jj)
             endif

          ! Calculate PON output:
           a(:,1) = Varout(oPON,:)   ! Total PON
           call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1),&
                PONout(nm,1))  

          else
           if(INCLUDESIZE) then
             if (jj .eq. 1) then
                nm = i- nrow(1,jj) - nrow(2,jj) - nrow(3,jj) - nrow(4,jj)
             else
                nm = 0
                do j = 1, (jj-1)
                   nm = nm + nrow(5, j)
                enddo
                nm=nm+i-nrow(1,jj)-nrow(2,jj)-nrow(3,jj)-nrow(4,jj)
             endif
          !  Calculate Size-fractionated Chl output:
             do j = 1, 4
                a(:,1) = Varout(oCHLs(j),:) 
                call gridinterpol(nlev,1,Z_r(:),a(:,1),1,depth(1), &
                     Sizeout(nm,j))
             enddo
            endif
          endif  ! ==> End of selecting data types
          deallocate(a)
        Endif !==> End of matching DOY
       Enddo  !==> End of loop of selecting obs. data
   Endif !==> End of checking if final year
  !---------------------------------------------------------------------
 ENDIF  !==> END of daily work
4 format(I10,1x,I7, <nlev+1>(2x,1pe12.3))

  ! Biological rhs:
  selectcase(Model_ID)
    case(NPZD2sp)
      call NPZD_2sp
    case(NPPZDD)
      call NPPZDD_MOD
    case(EFTPPDD)
      call EFT_PPZDD_MOD
    case(NPZDdisc)
  !    call NPZD_DISC
    case(Geiderdisc)
      call Geider_DISC
    case(EFTdisc)
      call FLEXEFT_DISC
    case(EFTcont)
      call FLEXEFT_CONT
    case(NPZDFix) 
      call NPZD_Fix
    case(NPZDFixIRON)
      call NPZD_Fix
    case(NPZDcont) 
      call NPZD_CONT
    case(EFTsimple)
      call FlexEFT_simple
    case(EFT2sp)
      call FLEXEFT_2SP
    case(EFTsimIRON)
      call FlexEFT_simple
    case(Geidersimple)
      call Geider_simple
    case(GeidsimIRON)
      call Geider_simple
    case default
      write(6,*) 'Error in choosing biological models! Quit.'
      stop
  endselect
  ! Pass the new state variables to Vars
  do j = 1,NVAR
     Vars(j,:)=Varout(j,:)
  enddo
  
  ! Diffusion:
  do j = 1,NVAR
     Vars1(:)     = Vars(j,:)
     call diff_center(nlev,dtsec,cnpar,1,Hz,Aks,Vars1,Vars2)
  
     ! Save diffusion fluxes (normalized to per day)
     Varout(oD_NO3+j-1,:) = (Vars2(:) - Vars1(:))/dtdays
  
     ! Update the state variables:
     do k = 1, nlev
        Vars(j,k)=max(Vars2(k), 1D-10)
     enddo
  enddo
  
  ! Sinking:

  do j = 1,NVsinkterms
     ww_(:)   = ww(:,j)
     Vars2(:) = Vars(Windex(j),:)
     call adv_center(nlev,dtsec,Hz,Hz,ww_(:),1,1,Yup,Ydw,6,mode1,Vars2(:))
     Vars(Windex(j),:) = Vars2(:)
  enddo
  
  ! Vertical advection (due to w, unit: m/s):
  ! Now do not consider vertical w
  !do j = 1,NVAR
  !   call adv_center(nlev,dtsec,Hz,Hz,w(:),  &
  !                   1,1,Yup,Ydw,6,mode0,Vars(j,:))
  !enddo
  
  !write(6,*) 'After advection, Vars(iCHL,1) = ',Vars(iCHL,1)
  
!  write(6,*) 'Everything OK for Sinking!' 
  ENDDO    ! ==> End of time stepping
  
  ! Close files:
  if (savefile) then
     do i = 1, (Nout+ow)
        close (unit=9+i)
     enddo
  endif
  deallocate(ww)
ENDDO  ! ==> End of Stn
!dCChl = dCChl * 1D-6
END subroutine Timestep
!=====================================================
subroutine End_model
  implicit none

!  ! Release memory:
  deallocate(Vars)
  deallocate(Varout)
  deallocate(params)
  deallocate(iPHY)
  deallocate(Labelout)
  deallocate(Windex)
  deallocate(oPHY)
  deallocate(oCHL)
  deallocate(oTheta)
  deallocate(oQN)
  deallocate(omuNET)
  deallocate(oGraz)
!  deallocate(ow_p)
  deallocate(oSI)
  deallocate(oLno3)
  deallocate(oD_PHY)
  deallocate(oD_CHL)

END subroutine End_model
!========================================================
subroutine setup_grid
implicit none

! Local scratch variable:
integer                :: i
real                   :: sc_r, C_sig

Z_w(0) = -hmax

!Following Song and Haidvogel (1994). sinh is the hyperbolic sin function
do i   = 1,nlev
    sc_r  = (float(i-nlev) - 0.5)/float(nlev)
   C_sig  = sinh(thetaS*sc_r)/sinh(thetaS)      ! -1 < C_sig < 0
   Z_r(i) = C_sig*hmax
    sc_r  = (float(i-nlev))/float(nlev)
   C_sig  = sinh(thetaS*sc_r)/sinh(thetaS)      ! -1 < C_sig < 0
   Z_w(i) = C_sig*hmax
    Hz(i) = Z_w(i) - Z_w(i-1)
enddo
end subroutine setup_grid
!===========
subroutine Calculate_PAR(I_0, nlev_, Hz, Chl, PAR)
  !top level is nlev_, bottom layer is 1, following ROMS convention
  implicit none
  real   , intent(in) :: I_0
  integer, intent(in) :: nlev_    ! Total number of vertical layers
  real   , intent(in) :: Hz(nlev_), Chl(nlev_)
  real   , intent(out):: PAR(nlev_)
  integer             :: i
  real                :: par0, attn    ! Scratch variable
  real   , parameter  :: kw = 0.04, kc = 0.025

  par0 = I_0   !Light at the grid surface
  do i = nlev_,1,-1
     attn   = exp(-0.5*Hz(i)*(kw + kc*Chl(i)))
     PAR(i) = par0*attn
     par0   = PAR(i)*attn
  enddo
end subroutine Calculate_PAR

!-----------------------------------------------------------------------
End MODULE MOD_1D
