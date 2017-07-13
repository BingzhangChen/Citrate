SUBROUTINE choose_model
USE bio_MOD
implicit none
integer             :: i
namelist /Model/    Stn, Model_ID, nutrient_uptake, grazing_formulation, bot_bound
character(len=100)  :: format_string

!  open the namelist file and read station name.
open(namlst,file='Model.nml',status='old',action='read')
read(namlst,nml=Model)
close(namlst)

if (Model_ID     == NPZDdisc) then
  write(6,*) 'Inflexible NPZD discrete model selected!'
  NPHY    = 20
else if (Model_ID== NPZDdiscFe) then
  write(6,*) 'Inflexible NPZD discrete model with Iron selected!'
  NPHY    = 20
  do_IRON = .TRUE.
else if (Model_ID== Geiderdisc) then
  write(6,*) 'Geider discrete model selected!'
  NPHY    = 20
else if (Model_ID== EFTdisc) then
  write(6,*) 'Flexible discrete model selected!'
  NPHY    = 20
else if (Model_ID== EFTdiscFe) then
  write(6,*) 'Flexible discrete model selected!'
  NPHY    = 20
  do_IRON = .TRUE.
else if (Model_ID==NPZDFix) then
  write(6,*) 'Inflexible NPZD model selected!'
  NPHY = 1
else if (Model_ID==NPZDN2) then
  write(6,*) 'NPZD model with N2 fixation selected!'
  NPHY = 1
  N2fix= .true.
else if (Model_ID==NPZDFixIRON) then
  write(6,*) 'Inflexible NPZD model with Iron selected!'
  do_IRON = .TRUE.
  NPHY    = 1
else if (Model_ID==Geidersimple) then
  write(6,*) 'Geider simple model selected!'
  NPHY = 1
else if (Model_ID==GeidsimIRON) then
  write(6,*) 'Geider simple model with Iron selected!'
  do_IRON = .TRUE.
  NPHY = 1
else if (Model_ID==EFTsimple) then
  write(6,*) 'Flexible simple model selected!'
  NPHY = 1
else if (Model_ID==EFTsimIRON) then
  write(6,*) 'Flexible simple model with Iron selected!'
  do_IRON = .TRUE.
  NPHY = 1
else if (Model_ID==EFTcont) then
  write(6,*) 'Flexible continous model selected!'
  NPHY = 1
else if (Model_ID==NPZDcont) then
  write(6,*) 'NPZD continous model selected!'
  NPHY    = 1
  DO_IRON = .TRUE.
else if (Model_ID==EFT2sp .OR. Model_ID==NPZD2sp) then
  write(6,*) 'Two species phytoplankton model selected!'
  NPHY = 2
else if (Model_ID==NPPZDD .OR. Model_ID==EFTPPDD) then
  write(6,*) 'Two phytoplankton two detritus model selected!'
  NPHY = 2
endif

allocate(iPHY(NPHY))
allocate(iCHL(NPHY))
allocate(oPHY(NPHY))
allocate(oCHL(NPHY))
allocate(omuNet(NPHY))
allocate(oLno3(NPHY))
allocate(oSI(NPHY))
allocate(oGraz(NPHY))
allocate(oD_PHY(NPHY))
allocate(oD_CHL(NPHY))
allocate(otheta(NPHY))
allocate(oQN(NPHY))
if (N2fix) allocate(oQp(NPHY))

do i=1,NPHY
   iPHY(i) = i + iNO3
enddo 

iZOO = iPHY(NPHY)+1
iDET = iZOO+1
if (Model_ID==Geiderdisc .or. Model_ID==Geidersimple .or. Model_ID==GeidsimIRON) then
   do i=1,NPHY
      iCHL(i) = i + iDET
   enddo 
   NVAR = iCHL(NPHY)
else if(Model_ID==EFTcont .or. Model_ID==NPZDcont) then
   iPMU = iDET  + 1
   iVAR = iPMU  + 1
   NVAR = iVAR 
   if (Model_ID == NPZDcont) then
      ifer = iVAR + 1
      NVAR = ifer
   endif
else if(Model_ID==NPPZDD .or. Model_ID==EFTPPDD) then
   iDET2= iDET + 1
   NVAR = iDET2
else if(Model_ID==NPZDN2) then
   iDETp = iDET + 1
   iPO4  = iDETp+ 1
   iDIA  = iPO4 + 1
   NVAR  = iDIA
else
   NVAR = iDET
endif

allocate(Vars(NVAR,nlev))
Vars(:,:)=0d0

if (Model_ID==Geiderdisc .or. Model_ID==Geidersimple .or. Model_ID==GeidsimIRON) then
    NVsinkterms = 1 + NPHY * 2  ! Include phyto and Chl
else if (Model_ID==NPPZDD .or. Model_ID==EFTPPDD .or. Model_ID==NPZDN2) then
    NVsinkterms = 2 + NPHY
else
    NVsinkterms = 1 + NPHY
endif

allocate(Windex(NVsinkterms))
do i=1,NPHY
   Windex(i)=iPHY(i)
   if (Model_ID==Geiderdisc .or. Model_ID==Geidersimple .or. Model_ID==GeidsimIRON) then
      Windex(i+NPHY)=iCHL(i)
   endif
enddo
if (Model_ID==NPPZDD .or. Model_ID==EFTPPDD) then
   Windex(NVsinkterms-1)=iDET
   Windex(NVsinkterms)  =iDET2
elseif (Model_ID==NPZDN2) then
   Windex(NVsinkterms-1)=iDET
   Windex(NVsinkterms)  =iDETp
else
   Windex(NVsinkterms)=iDET
endif

! Output array matrices (the order must be consistent with Vars)
do i=1,NPHY
   oPHY(i)=i+oNO3
enddo
oZOO =oPHY(NPHY)+1
oDET =oZOO+1
if(Model_ID==EFTcont .or. Model_ID==NPZDcont) then
   oPMU=oDET+1
   oVAR=oPMU+1
   if (do_IRON) then
      ofer = oVAR+1
      do i=1,NPHY
         oCHL(i) = i + ofer
      enddo 
   else
      do i=1,NPHY
         oCHL(i) = i + oVAR
      enddo 
   endif
else if (Model_ID == NPPZDD .or. Model_ID==EFTPPDD) then
   oDET2=oDET+1
   do i=1,NPHY
      oCHL(i) = i + oDET2
   enddo 
else if (Model_ID == NPZDN2) then
   oDETp=oDET+1
   oPO4 =oDETp+1
   oDIA =oPO4 +1
   oPOP =oDIA +1
   oDIAu=oPOP +1
   do i=1,NPHY
      oCHL(i) = i + oDIAu
   enddo 
else
  do i=1,NPHY
     oCHL(i) = i + oDET
  enddo 
endif

! The above must match with i** indeces   
oPHYt=oCHL(NPHY)+1
oCHLt=oPHYt+1

if (Model_ID==Geiderdisc .or. Model_ID==NPZDdisc .or. Model_ID==EFTdisc&
.or.Model_ID==EFTcont .or. Model_ID==NPZDcont) then
   do i=1,4
      oCHLs(i)=oCHLt+i
   enddo
endif

!allocate(ow_p(NPHY))

if (Model_ID==Geiderdisc .or. Model_ID==NPZDdisc .or. Model_ID==EFTdisc&
.or.Model_ID==EFTcont    .or. Model_ID==NPZDcont) then
   do i=1,NPHY
      omuNet(i)= oCHLs(4) + i
   enddo
else
   do i=1,NPHY
      omuNet(i)= oCHLt + i
   enddo
end if

do i=1,NPHY
   oGraz(i) = omuNet(NPHY) + i
enddo

do i=1,NPHY
   oLno3(i) = oGraz(NPHY)  + i
enddo

do i=1,NPHY
   oSI(i)=oLno3(NPHY)+i
enddo
do i=1,NPHY
   oQN(i)=oSI(NPHY)+i
enddo

if (N2fix) then
   do i=1,NPHY
      oQp(i)=oQN(NPHY)+i
   enddo
   do i=1,NPHY
      otheta(i)=oQp(NPHY)+i
   enddo
else
   do i=1,NPHY
      otheta(i)=oQN(NPHY)+i
   enddo
endif
oZ2N=otheta(NPHY)+1
oD2N=oZ2N+1
oPPt=oD2N+1
oPON=oPPt+1
oPAR_ =oPON+1
omuAvg=oPAR_+1
oD_NO3=omuAvg+1

do i=1,NPHY
   oD_PHY(i)=oD_NO3+i
enddo

oD_ZOO=oD_PHY(NPHY)+1
oD_DET=oD_ZOO+1
if (Model_ID==Geiderdisc .or. Model_ID==Geidersimple .or. Model_ID==GeidsimIRON) then
   do i=1,NPHY
      oD_CHL(i)=oD_DET+1
   enddo
   Nout=oD_CHL(NPHY)
else if(Model_ID==NPPZDD .or. Model_ID==EFTPPDD) then
   oD_DET2=oD_DET+1
   Nout=oD_DET2
else if(Model_ID==EFTcont .or. Model_ID==NPZDcont) then
   oD_PMU=oD_DET+1
   oD_VAR=oD_PMU+1
   oD_fer=oD_VAR+1
   od2mu =oD_fer+1
   odmudl=od2mu +1
   od3mu =odmudl+1
   od4mu =od3mu +1
   Nout  =od4mu
else if(Model_ID==NPZDN2) then
   oD_DETp=oD_DET +1
   oD_PO4 =oD_DETp+1
   oD_DIA =oD_PO4+1
   Nout   =oD_DIA
else
   Nout   =oD_DET
endif
allocate(Varout(Nout,nlev))
IF (AllocateStatus /= 0) STOP "*** Error in allocating Varout ***"
allocate(Labelout(Nout+ ow ))

Labelout(oTemp  )='Temp'
Labelout(oPAR   )='PAR '
Labelout(oAks   )='Aks '
Labelout(oDust  )='Dust'
Labelout(ow     )='w   '
Labelout(oNO3+ow)='NO3 '
do i=1,NPHY
   if (i < 10) then
       format_string = "(A3,I1)"
   else
       format_string = "(A3,I2)"
   endif
   write(Labelout(oPHY(i)  +ow),  format_string) 'PHY',i
   write(Labelout(oSI(i)   +ow),  format_string) 'SI_',i
   write(Labelout(oQN(i)   +ow),  format_string) 'QN_',i
   if(N2fix) write(Labelout(oQP(i)+ow), format_string) 'QP_',i
   write(Labelout(oLno3(i) +ow),  format_string) 'Lno',i
   write(Labelout(otheta(i)+ow),  format_string) 'The',i
   write(Labelout(oCHL(i)  +ow),  format_string) 'CHL',i
enddo

Labelout(oZOO +ow)='ZOO'
Labelout(oDET +ow)='DET'
if (Model_ID==NPPZDD  .or. Model_ID==EFTPPDD) Labelout(oDET2 +ow)='DET2'
!
if (Model_ID==NPZDN2) then
    Labelout(oDETp +ow)='DETp'
    Labelout(oPO4  +ow)='DIP'  ! Consistent with data file
    Labelout(oPOP  +ow)='POP'
    Labelout(oDIA  +ow)='DIA'
    Labelout(oDIAu +ow)='uDIA'
    Labelout(oD_DETp+ow)='DDETp'
    Labelout(oD_PO4 +ow)='D_DIP'
    Labelout(oD_DIA +ow)='D_DIA'
endif
if (Model_ID==EFTcont .or. Model_ID==NPZDcont) then
    Labelout(oPMU  +ow)='PMU'
    Labelout(oVAR  +ow)='VAR'
    Labelout(odmudl+ow)='dmudl'
    Labelout(od2mu +ow)='d2mu '
    Labelout(od3mu +ow)='d3mu'
    Labelout(od4mu +ow)='d4mu'
    Labelout(oD_PMU+ow)='D_PMU'
    Labelout(oD_VAR+ow)='D_VAR'
    if (do_IRON) then
       Labelout(ofer  +ow)='Fer'
       Labelout(oD_fer+ow)='D_Fe'
    endif
endif
Labelout(oPHYt+ow)='PHY_T'
Labelout(oCHLt+ow)='CHL_T'

do i=1,NPHY
   if (i < 10) then
       format_string = "(A3,I1)"
   else
       format_string = "(A3,I2)"
   endif

   write(Labelout(omuNet(i) + ow), format_string) 'muN',i
   write(Labelout(oGraz(i)  + ow), format_string) 'Gra',i
   write(Labelout(oD_PHY(i) + ow), format_string) 'D_P',i
   if (Model_ID==GeidsimIRON.or.Model_ID==Geiderdisc .or. Model_ID==Geidersimple) then
      write(Labelout(oD_CHL(i) + ow), format_string) 'DCH',i
   endif
enddo

Labelout(oZ2N  + ow)='Z2N'
Labelout(oD2N  + ow)='D2N'
Labelout(oPPt  + ow)='NPP_T'
Labelout(oPON  + ow)='PON'
Labelout(oPAR_ + ow)='PAR_'
Labelout(omuAvg+ ow)='muAvg'
Labelout(oD_NO3+ ow)='D_NO3'
Labelout(oD_ZOO+ ow)='D_ZOO'
Labelout(oD_DET+ ow)='D_DET'
if(Model_ID==NPPZDD.or.Model_ID==EFTPPDD) Labelout(oD_DET2+ow)='DDET2'
if(Model_ID==Geiderdisc.or.Model_ID==NPZDdisc .or.Model_ID==EFTdisc .or.&
   Model_ID==EFTcont .or. Model_ID==NPZDCONT) then
   do i=1,4
      if (i < 10) then
          format_string = "(A4,I1)"
      else
          format_string = "(A4,I2)"
      endif
      write(Labelout(oCHLs(i) +ow), format_string) 'CHLs',i
   enddo
endif

do i = 1, Nout+ow
   write(6,*) 'Labelout(',i,') = ',trim(Labelout(i))
enddo

! Initialize parameters
! Indices for parameters that will be used in MCMC                 
! For EFT models, the affinity approach not used for now
! Need to have tradeoffs for maximal growth rate (mu0) and Kn
! Common parameters:
imu0    =  1
igmax   =  imu0+1
if (Model_ID == EFT2sp .OR. Model_ID==EFTPPDD .OR. Model_ID==NPZD2sp .OR. Model_ID==NPPZDD) then
    imu0B   =  igmax + 1  ! The ratio of mu0 of the second species to the first
    iaI0B   =  imu0B + 1  ! The ratio of aI0 of the second species to the first
    if (Model_ID==NPZD2sp .OR. Model_ID==NPPZDD) then
       ibI0B=  iaI0B  + 1
       imz  =  ibI0B  + 1
    else
       iaI0 =  iaI0B  + 1
       imz  =  iaI0   + 1
    endif
else if(Model_ID==NPZDN2) then
      iKPHY =  igmax  + 1  
       imz  =  iKPHY  + 1
else if(Model_ID==NPZDcont .or. Model_ID==NPZDdisc .or. Model_ID==NPZDFix) then
      imz   =  igmax  + 1
else
    iaI0    =  igmax  + 1
    imz     =  iaI0   + 1
endif

if (nutrient_uptake .eq. 1) then
   iKN     =  imz    + 1
   if (Model_ID==NPZD2sp .or. Model_ID==NPPZDD) then
      iKN2 =  iKN    + 1
      iQ0N =  iKN2   + 1
   else if (Model_ID==NPZDN2) then
     ! iKP  =  iKN    + 1
     iKPnif=  iKN    + 1
     iLnifp=  iKPnif + 1
     iRDN_N=  iLnifp + 1
     iRDN_P=  iRDN_N + 1
      iQ0N =  iRDN_P + 1
   else
      iQ0N =  iKN    + 1
   endif
elseif (nutrient_uptake.eq.2) then
   if (Model_ID==NPZD2sp .or. Model_ID==NPZDFix .or. Model_ID==NPPZDD .or. Model_ID==NPZDcont .or. Model_ID==NPZDN2) then
      write(6,*) 'We do not use affinity-based equation for NPZD model!'
      stop
   endif
   iA0N    =  imz   + 1
   if (Model_ID==EFT2sp .or. Model_ID == EFTPPDD) then
      iA0N2 =iA0N+1   ! The ratio of A0N of the second species to the first  
     ialphaG=iA0N2+1
      iQ0N  =ialphaG+1
   else   
      iQ0N = iA0N + 1
   endif
endif

if (Model_ID == NPPZDD .or. Model_ID == EFTPPDD) then
   itau   = iQ0N  +1
   iwDET2 = itau  +1
   iwDET  = iwDET2+1
else
   iwDET  = iQ0N  +1
endif

if (Model_ID==NPZDFix .or.Model_ID==NPZD2sp    .or.Model_ID==NPPZDD  &
.or.Model_ID==NPZDdisc.or.Model_ID==NPZDFixIRON.or.Model_ID==NPZDcont&
.or.Model_ID==NPZDN2) then
   iaI0_C  =  iwDET    + 1
   if (Model_ID==NPZDFix .or. Model_ID==NPZDN2) then
      NPar = iaI0_C
   else if (Model_ID ==NPZD2sp .or. Model_ID == NPPZDD) then
      iRL2 = iaI0_C + 1
    ialphaG= iRL2   + 1
      NPar = ialphaG
   else if (Model_ID == NPZDFixIRON) then
      iKFe = iaI0_C  + 1
      NPar = iKFe
   else
      ialphamu=iaI0_C+1
      ibetamu =ialphamu+1
      ialphaI =ibetamu+1
      !ialphaG =ialphaI +1
      if (nutrient_uptake.eq.1) then
         ialphaKN=ialphaI+1
         if (Model_ID == NPZDcont) then
           
           iVTR   =ialphaKN+1
           !igb    =iVTR    +1
           if (do_IRON) then
             iKFe =iVTR    +1
          ialphaFe=iKFe    +1
          idustsol=ialphaFe+1
             NPar =idustsol
           else
             NPar =iVTR
           endif
         else
           NPar   =ialphaKN
         endif
      elseif(nutrient_uptake.eq.2) then
         ialphaA  =ialphaG+1
           NPar   =ialphaA
      endif
   endif

else if(Model_ID==Geiderdisc.or.Model_ID==EFTdisc.or.Model_ID==EFTcont) &
then
   ialphaI     =iwDET   + 1
   ialphaG     =ialphaI + 1
   ialphamu    =ialphaG + 1
   if (nutrient_uptake.eq.1) then
       ialphaKN   =ialphamu+1
           NPar   =ialphaKN
   elseif(nutrient_uptake.eq.2) then
       ialphaA     =ialphamu+1
            NPar   =ialphaA
   endif
else if (Model_ID==GeidsimIRON .or. Model_ID==EFTsimIRON) then
  iKFe = iQ0N + 1
  NPar = iKFe
else if (Model_ID==EFT2sp .or. Model_ID==EFTPPDD) then
  iRL2 = iwDET + 1 ! The grazing preference on the second species (lower grazing impact)
  NPar = iRL2
else
  NPar = iwDET
endif

write(6,'(I2,1x,A20)') NPar,'parameters in total to be estimated.'
allocate(params(NPar))
allocate(ParamLabel(NPar))

ParamLabel(imu0) = 'mu0hat '
ParamLabel(imz)  = 'mz'
ParamLabel(igmax)= 'gmax'

if (Model_ID .eq. NPZDN2) then
   params(imz)   = 0.15*16d0
else
   params(imz)   = 0.15
endif
params(igmax)    = 1.0
if (Model_ID == NPZDdisc .or. Model_ID == NPZDFix .or. Model_ID==NPZDcont &
.or.Model_ID == NPZDN2) then
   params(imu0)     = 1.2
else
   params(imu0)     = 2.5d0
endif

if(Model_ID==EFT2sp .or. Model_ID==EFTPPDD .or.  Model_ID==NPZD2sp .or. Model_ID==NPPZDD) then
   ParamLabel(imu0B)='mu0B'
   params(imu0B)=0.3
endif

if(Model_ID==NPZDdisc.or.Model_ID==Geiderdisc.or.  &
   Model_ID==EFTdisc .or.Model_ID==EFTcont   .or.  &
   Model_ID==NPZDCONT) then

   ParamLabel(ialphamu) = 'alphamu'
   ParamLabel(ialphaI)  = 'alphaI'

   params(ialphamu) = 0.25
   params(ialphaI ) = 0.08

   !if (Model_ID .ne. NPZDCONT) then
   !   ParamLabel(ialphaQ) = 'alphaQ'
   !   params(ialphaQ)     = -0.05
   !endif
   if (nutrient_uptake.eq.1) then

     ParamLabel(ialphaKN)= 'alphaKN'
         params(ialphaKN)= 0.27

   else if (nutrient_uptake.eq.2) then
     ParamLabel(ialphaA) = 'alphaA'
     params(ialphaA)     = -0.3
   endif
endif
if (Model_ID.eq.NPZD2sp .OR. Model_ID.eq.NPPZDD) then
   ParamLabel(ibI0B)='bI0B'
      params(ibI0B) = -8d0
endif

if (nutrient_uptake .eq. 1) then
   ParamLabel(iKN)  = 'KN'
      params(iKN )  = 5d-1
  if (Model_ID .eq. NPZDN2) then
   !ParamLabel(iKP)  = 'KP'
  ParamLabel(iKPnif)= 'KPnif'
   !   params(iKP)   = .05
      params(iKPnif)= 1D-3
  ParamLabel(iLnifp)= 'Lnifp'
      params(iLnifp)= 0.17*16.
  ParamLabel(iKPHY) = 'KPHY'
      params(iKPHY) = .5/16d0
  ParamLabel(iRDN_N)= 'RDNn'
      params(iRDN_N)= 0.05
  ParamLabel(iRDN_P)= 'RDNp'
      params(iRDN_P)= 0.1
  endif
  if (Model_ID .eq. NPZD2sp .or. Model_ID.eq.NPPZDD) then
   ParamLabel(iKN2) = 'KN2'
      params(iKN2)  = 1D0
  endif
else if (nutrient_uptake .eq. 2) then
   ParamLabel(iA0N    ) = 'A0N    '
   params(iA0N)         = 5d0
  if (Model_ID .eq. EFT2sp .or. Model_ID .eq. EFTPPDD) then
     ParamLabel(iA0N2)='A0N2'   
     params(iA0N) = 7d1
     params(iA0N2)= -1d0  ! It is a ratio after log10 transformation
  endif
endif

if (Model_ID .eq. EFT2sp .or. Model_ID .eq. EFTPPDD .or. Model_ID .eq. NPZD2sp .or. Model_ID.eq.NPPZDD) then
   ParamLabel(iaI0B)='aI0B'
   ParamLabel(iRL2)='RL2'
   params(iaI0B)=0.3  ! A ratio after log transformation
   params(iRL2) =0.5
endif

if(Model_ID.eq.NPZD2sp .OR. Model_ID.eq.EFTdisc .OR. Model_ID .eq. EFT2sp   &
  .OR. Model_ID.eq.NPPZDD .OR. Model_ID .eq. EFTPPDD) then
   ParamLabel(ialphaG)='alphaG'
   params(ialphaG)=-5d0
endif
!
ParamLabel(iwDET) = 'wDET   '
params(iwDET)     = 6d-1  ! wDET = 10**params(iwDET)
!
if (Model_ID == NPPZDD .or. Model_ID == EFTPPDD) then
   ParamLabel(iwDET2)= 'wDET2'
   params(iwDET2)=1d0

   ParamLabel(itau)= 'Tau'
   params(itau)=5D-3
endif
!
if(Model_ID==NPZDdisc.or.Model_ID==NPZD2sp.or.Model_ID==NPPZDD.or.Model_ID==NPZDFix .or. Model_ID==NPZDFixIRON .or. Model_ID==NPZDcont.or. Model_ID==NPZDN2) then
  ParamLabel(iaI0_C)='aI0_C'
  params(iaI0_C)    =0.055
  if (Model_ID == NPZDcont) then
     !ParamLabel(igb)='gb'
     !params(igb)    =-1D-12
     ParamLabel(iVTR)='VTR'
     params(iVTR)    =0.08
     ParamLabel(ibetamu)='betamu'
     params(ibetamu)    =-0.025
     if (do_IRON) then
            ParamLabel(iKFe)='KFe'
        ParamLabel(ialphaFe)='alphaFe'
        ParamLabel(idustsol)='dustsol'
                params(iKFe)=0.08
            params(ialphaFe)=0.27
            params(idustsol)=0.02
     endif
  endif
endif

if(Model_ID==GeidsimIRON .or. Model_ID==EFTsimIRON .or. Model_ID==NPZDFixIRON) then
  ParamLabel(iKFe)  ='KFe'
  params(iKFe)      =0.08
endif

ParamLabel(iQ0N   ) = 'Q0N    '
   params(iQ0N)     = 0.06

if (Model_ID .ne. NPZDFix .and. Model_ID .ne. NPZDcont .and. Model_ID .ne. NPZDdisc .and. Model_ID .ne. NPZD2sp .and. Model_ID .ne. NPPZDD .and. Model_ID .ne. NPZDN2) then
ParamLabel(iaI0    ) = 'aI0'
   params(iaI0   )   = 0.2      ! aI0_Chl, Chl-specific P-I slope
endif

if (Model_ID==NPZDdisc.or.Model_ID==Geiderdisc.or.  &
    Model_ID==EFTdisc) then
    call assign_PMU
endif
end subroutine choose_model
