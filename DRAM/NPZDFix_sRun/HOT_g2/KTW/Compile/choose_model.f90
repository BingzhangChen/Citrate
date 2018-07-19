subroutine choose_model
use bio_MOD
implicit none
integer             :: i
namelist /Model/    Stn, Model_ID, nutrient_uptake, grazing_formulation
character(len=100)  :: format_string

!  open the namelist file and read station name.
open(namlst,file='Model.nml',status='old',action='read')
read(namlst,nml=Model)
close(namlst)

if (Model_ID==NPZDdisc) then
  write(6,*) 'Inflexible NPZD discrete model selected!'
  NPHY = 20
else if (Model_ID==NPZDdiscFe) then
  write(6,*) 'Inflexible NPZD discrete model with Iron selected!'
  NPHY    = 20
  do_IRON = .TRUE.
else if (Model_ID==Geiderdisc) then
  write(6,*) 'Geider discrete model selected!'
  NPHY = 20
else if (Model_ID==EFTdisc) then
  write(6,*) 'Flexible discrete model selected!'
  NPHY = 20
else if (Model_ID==EFTdiscFe) then
  write(6,*) 'Flexible discrete model selected!'
  NPHY = 20
  do_IRON = .TRUE.
else if (Model_ID==NPZDFix) then
  write(6,*) 'Inflexible NPZD model selected!'
  NPHY = 1
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
else if(Model_ID==EFTcont) then
   iPMU = iDET      + 1
   iVAR = iPMU      + 1
   NVAR = iVAR 
else
   NVAR = iDET
endif

allocate(Vars(NVAR,nlev))
Vars(:,:)=0d0

if (Model_ID==Geiderdisc .or. Model_ID==Geidersimple .or. Model_ID==GeidsimIRON) then
    NVsinkterms = 1 + NPHY * 2  ! Include phyto and Chl
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
Windex(NVsinkterms)=iDET

! Output array matrices (the order must be consistent with Vars)
do i=1,NPHY
   oPHY(i)=i+oNO3
enddo
oZOO =oPHY(NPHY)+1
oDET =oZOO+1
if(Model_ID==EFTcont) then
   oPMU=oDET+1
   oVAR=oPMU+1
   do i=1,NPHY
      oCHL(i) = i + oVAR
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
.or.Model_ID==EFTcont) then
   do i=1,4
      oCHLs(i)=oCHLt+i
   enddo
endif

!allocate(ow_p(NPHY))

if (Model_ID==Geiderdisc .or. Model_ID==NPZDdisc .or. Model_ID==EFTdisc&
.or.Model_ID==EFTcont) then
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

do i=1,NPHY
   otheta(i)=oQN(NPHY)+i
enddo
oZ2N=otheta(NPHY)+1
oD2N=oZ2N+1
oPPt=oD2N+1
oD_NO3=oPPt +1

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
else if(Model_ID==EFTcont) then
   oD_PMU=oD_DET+1
   oD_VAR=oD_PMU+1
   Nout=oD_VAR
else
   Nout=oD_DET
endif
allocate(Varout(Nout,nlev))
IF (AllocateStatus /= 0) STOP "*** Error in allocating Varout ***"
allocate(Labelout(Nout+ ow ))

Labelout(oTemp  )='Temp'
Labelout(oPAR   )='PAR '
Labelout(oAks   )='Aks '
Labelout(oDFe   )='DFe '
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
   write(Labelout(oLno3(i) +ow),  format_string) 'Lno',i
   write(Labelout(otheta(i)+ow),  format_string) 'The',i
   write(Labelout(oCHL(i)  +ow),  format_string) 'CHL',i
enddo

Labelout(oZOO +ow)='ZOO'
Labelout(oDET +ow)='DET'
if (Model_ID==EFTcont) then
    Labelout(oPMU+ow)='PMU'
    Labelout(oVAR+ow)='VAR'
    Labelout(oD_PMU+ow)='D_PMU'
    Labelout(oD_VAR+ow)='D_VAR'
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
Labelout(oD_NO3+ ow)='D_NO3'
Labelout(oD_ZOO+ ow)='D_ZOO'
Labelout(oD_DET+ ow)='D_DET'

if(Model_ID==Geiderdisc.or.Model_ID==NPZDdisc .or.Model_ID==EFTdisc .or.&
   Model_ID==EFTcont) then
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
iEp     =  1
iEz     =  iEp    + 1
imu0    =  iEz    + 1
iaI0    =  imu0   + 1
iQ0N    =  iaI0   + 1

if (nutrient_uptake .eq. 1) then
   iKN     =  iQ0N   + 1
   igmax   =  iKN    + 1
elseif (nutrient_uptake.eq.2) then
   iA0N    =  iQ0N   + 1
   igmax   =  iA0N   + 1
endif

ikp     =  igmax  + 1 
iwDET   =  ikp    + 1
irdN    =  iwDET  + 1
imz     =  irdN   + 1

if (Model_ID==NPZDFix.or.Model_ID==NPZDdisc.or.Model_ID==NPZDFixIRON) then
   itheta  =  imz    + 1
   if (Model_ID==NPZDFix) then
      NPar = itheta
   else if(Model_ID == NPZDFixIRON) then
      iKFe = itheta  + 1
      NPar = iKFe
   else
      ialphamu=itheta+1
      if (nutrient_uptake.eq.1) then
         ialphaKN=ialphamu+1
         if(kill_the_winner) then
           ialphaG=ialphaKN+1
           NPar   =ialphaG
         else
           NPar   =ialphaKN
         endif
      elseif(nutrient_uptake.eq.2) then
         ialphaA  =ialphamu+1
         if(kill_the_winner) then
           ialphaG=ialphaA+1
           NPar   =ialphaG
         else
           NPar   =ialphaA
         endif
      endif
   endif

else if(Model_ID==Geiderdisc.or.Model_ID==EFTdisc.or.Model_ID==EFTcont) &
then
   ialphamu    =imz + 1
   if (nutrient_uptake.eq.1) then
       ialphaKN    =ialphamu+1
       if(kill_the_winner) then
            ialphaG=ialphaKN+1
            NPar   =ialphaG
       else
            NPar   =ialphaKN
       endif
   elseif(nutrient_uptake.eq.2) then
       ialphaA     =ialphamu+1
       if(kill_the_winner) then
            ialphaG=ialphaA+1
            NPar   =ialphaG
       else
            NPar   =ialphaA
       endif
   endif
else if (Model_ID==GeidsimIRON .or. Model_ID==EFTsimIRON) then
  iKFe = imz + 1
  NPar = iKFe
else
  NPar = imz
endif

write(6,'(I2,1x,A20)') NPar,'parameters in total to be estimated.'
allocate(params(NPar))
allocate(ParamLabel(NPar))

ParamLabel(imu0) = 'mu0hat '

if (Model_ID == NPZDdisc) then
   ! Update parameters on Nov. 9 2016
   params(imu0)     = 2.5D0

elseif (Model_ID == EFTdisc) then

   params(imu0)     = 5d0

else

   params(imu0)     = 2.5d0

endif

if(Model_ID==NPZDdisc.or.Model_ID==Geiderdisc.or.  &
   Model_ID==EFTdisc .or.Model_ID==EFTcont) then

   ParamLabel(ialphamu) = 'alphamu'

   if (Model_ID .eq. NPZDdisc) then
       params(ialphamu) = 0.08
   else
       params(ialphamu) = 0.1
   endif

   if (nutrient_uptake.eq.1) then

     ParamLabel(ialphaKN)= 'alphaKN'

     if (Model_ID .eq. NPZDdisc) then
         params(ialphaKN)= 0.26
     else
         params(ialphaKN)= 0.3
     endif

   else if (nutrient_uptake.eq.2) then

     ParamLabel(ialphaA) = 'alphaA'
     params(ialphaA)     = -0.3

   endif
endif

if (nutrient_uptake .eq. 1) then
   ParamLabel(iKN     ) = 'KN     '
   if (Model_ID .eq. NPZDdisc) then
      params(iKN    )   = 0.11d0
   else
      params(iKN    )   = 1d0
   endif

else if (nutrient_uptake .eq. 2) then
   ParamLabel(iA0N    ) = 'A0N    '
   params(iA0N)         = 1d0
endif

ParamLabel(igmax   ) = 'gmax   '
params(igmax  ) = 1d0

ParamLabel(iwDET   ) = 'wDET   '
params(iwDET  ) = 2D0

ParamLabel(ikp     ) = 'kp     '
params(ikp    ) = 5d-1

ParamLabel(irdn    ) = 'rdn    '
params(irdn   ) = 0.1d0

ParamLabel(imz     ) = 'mz     '
params(imz         ) = 0.15D0

ParamLabel(iEp     ) = 'Ep     '
params(iEp         ) = 0.5

ParamLabel(iEz     ) = 'Ez     '
params(iEz         ) = 0.6

if(Model_ID==NPZDdisc.or.Model_ID==NPZDFix .or. Model_ID==NPZDFixIRON) then
  ParamLabel(itheta)='theta'
  params(itheta)    =0.24
endif

if(Model_ID==GeidsimIRON .or. Model_ID==EFTsimIRON .or. Model_ID==NPZDFixIRON) then
  ParamLabel(iKFe)  ='KFe'
  params(iKFe)      =0.08
endif

ParamLabel(iQ0N    ) = 'Q0N    '

if (Model_ID == NPZDdisc .or. Model_ID == NPZDFix .or. Model_ID == NPZDFixIRON) then
   params(iQ0N ) = 0.15
else
   params(iQ0N ) = 0.076
endif

ParamLabel(iaI0    ) = 'aI0    '

if (Model_ID .eq. NPZDdisc) then
   params(iaI0   )   = 0.01
else
   params(iaI0   )   = 0.2
endif

if((Model_ID==NPZDdisc.or.Model_ID==Geiderdisc.or. &
    Model_ID==EFTdisc .or.Model_ID==EFTcont).and.kill_the_winner) then
   ParamLabel(ialphaG ) = 'alphaG '
   params(ialphaG)=1.1
endif

if (Model_ID==NPZDdisc.or.Model_ID==Geiderdisc.or.                  &
    Model_ID==EFTdisc) then
    call assign_PMU
endif
do i = 1, NPar
   write(6,*) 'ParamLabel(',i,') = ',trim(ParamLabel(i))
enddo
end subroutine choose_model
