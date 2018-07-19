subroutine choose_model
use bio_MOD
implicit none
integer             :: i
namelist /Model/    Stn, Model_ID, nutrient_uptake, grazing_formulation

!  open the namelist file and read station name.
open(namlst,file='Model.nml',status='old',action='read')
read(namlst,nml=Model)
close(namlst)

SELECT CASE(model_ID)
  case(NPZDFix)
    write(6,*) 'Inflexible NPZD model selected!'
    NPHY = 1
    allocate(iPHY(NPHY))
    do i=1,NPHY
       iPHY(i)=i+iNO3
    enddo 

    iZOO=iPHY(NPHY)+1
    iDET=iZOO+1
    NVAR=iDET

    allocate(Vars(NVAR,nlev))

    NVsinkterms = 1 + NPHY
    allocate(Windex(NVsinkterms))
    do i=1,NPHY
       Windex(i)=iPHY(i)
    enddo
    Windex(NVsinkterms)=iDET

    ! Output array matrices
    allocate(oPHY(NPHY))
    do i=1,NPHY
       oPHY(i)=i+oNO3
    enddo
    oZOO =oPHY(NPHY)+1
    oDET =oZOO+1
    oCHLt=oDET+1
    ! The above must match with i** indeces   
 
    allocate(omuNet(NPHY))
    allocate(oGraz(NPHY))
    !allocate(ow_p(NPHY))
    allocate(oD_PHY(NPHY))
    allocate(oSI(NPHY))
    allocate(oLno3(NPHY))

    do i=1,NPHY
       omuNet(i)=oDET + i
    enddo
    oSI(1)   = omuNet(1)+1
    oLno3(1) = oSI(1) +1
    do i=1,NPHY
       oGraz(i)=oLno3(NPHY)+i
    enddo

    oZ2N=oGraz(NPHY)+1
    oD2N=oZ2N+1

    oPPt  =oD2N+1
    oD_NO3=oPPt +1

    do i=1,NPHY
       oD_PHY(i)=oD_NO3+i
    enddo

    oD_ZOO=oD_PHY(NPHY)+1
    oD_DET=oD_ZOO+1
    Nout  =oD_DET

    allocate(Varout(  Nout,nlev))
    allocate(Labelout(Nout+ ow ))

    Labelout(oTemp  )='Temp '
    Labelout(oPAR   )='PAR  '
    Labelout(oAks   )='Aks  '
    Labelout(ow     )='w    '
    Labelout(oNO3+ow)='NO3  '

    do i=1,NPHY
       write(Labelout(oPHY(i)+ow), "(A3,I2)") 'PHY',i
       write(Labelout(oSI(i)+ow),  "(A3,I2)") 'SI ',i
       write(Labelout(oLno3(i)+ow),"(A3,I2)") 'LNO',i
    enddo

    Labelout(oZOO +ow)='ZOO  '
    Labelout(oDET +ow)='DET  '
    Labelout(oCHLt+ow)='CHL_T'
    do i=1,NPHY
       write(Labelout(omuNet(i) + ow), "(A3,I2)") 'muN',i
       write(Labelout(oGraz(i)  + ow), "(A3,I2)") 'Gra',i
       write(Labelout(oD_PHY(i) + ow), "(A3,I2)") 'D_P',i
    enddo
  
    Labelout(oZ2N  + ow)='Z2N  '
    Labelout(oD2N  + ow)='D2N  '
    Labelout(oPPt  + ow)='NPP_T'
    Labelout(oD_NO3+ ow)='D_NO3'
    Labelout(oD_ZOO+ ow)='D_ZOO'
    Labelout(oD_DET+ ow)='D_DET'
    Labelout(oD_CHL+ ow)='D_CHL'

    do i = 1, Nout+ow
       write(6,*) 'Labelout(',i,') = ',trim(Labelout(i))
    enddo
    ! Initialize parameters
    ! Indices for parameters that will be used in MCMC                 
    imu0    =  1
    iaI0    =  imu0   + 1
    iQ0N    =  iaI0   + 1
    iKN     =  iQ0N   + 1
    itheta  =  ikN    + 1
    igmax   =  itheta + 1
    ikp     =  igmax  + 1 
    iwDET   =  ikp    + 1
    irdN    =  iwDET  + 1
    imz     =  irdN   + 1
    iEp     =  imz    + 1
    iEz     =  iEp    + 1
    NPar    =  iEz
    write(6,'(I2,1x,A20)') NPar,'parameters in total to be estimated.'
    allocate(params(NPar))
    allocate(ParamLabel(NPar))

    ParamLabel(imu0  ) = 'mu0hat '
    ParamLabel(iaI0  ) = 'aI0    '
    ParamLabel(iKN   ) = 'KN     '
    ParamLabel(iQ0N  ) = 'Q0N    '
    ParamLabel(itheta) = 'theta  '
    ParamLabel(iwDET ) = 'wDET   '
    ParamLabel(igmax ) = 'gmax   '
    ParamLabel(ikp   ) = 'kp     '
    ParamLabel(irdn  ) = 'rdn    '
    ParamLabel(imz   ) = 'mz     '
    ParamLabel(iEp   ) = 'Ep     '
    ParamLabel(iEz   ) = 'Ez     '

    params(imu0   ) = 2d0
    params(iaI0   ) = 0.5
    params(iKN    ) = 1d0
    params(iEp    ) = 0.4
    params(iEz    ) = 0.6
    params(igmax  ) = 1d0
    params(ikp    ) = 5d-1
    params(iQ0N   ) = 0.15
    params(itheta ) = 0.02
    params(iwDET  ) = 1D0
    params(irdn   ) = 5d-2
    params(imz    ) = 5d-2

  case(NPZDdisc)
    write(6,*) 'Inflexible NPZD discrete model selected!'
    NPHY = 20

    allocate(iPHY(NPHY))
    do i=1,NPHY
       iPHY(i) = i + iNO3
    enddo 

    iZOO = iPHY(NPHY)+1
    iDET = iZOO+1
    NVAR = iDET

    allocate(Vars(NVAR,nlev))

    NVsinkterms = 1 + NPHY
    allocate(Windex(NVsinkterms))
    do i=1,NPHY
       Windex(i)=iPHY(i)
    enddo
    Windex(NVsinkterms)=iDET

    ! Output array matrices
    allocate(oPHY(NPHY))
    do i=1,NPHY
       oPHY(i)=i+oNO3
    enddo
    oZOO =oPHY(NPHY)+1
    oDET =oZOO+1
    oPHYt=oDET+1
    oCHLt=oPHYt+1
    do i=1,4
       oCHLs(i)=oCHLt+i
    enddo

    ! The above must match with i** indeces   
 
    allocate(omuNet(NPHY))
    allocate(oGraz(NPHY))
    !allocate(ow_p(NPHY))
    allocate(oD_PHY(NPHY))

    do i=1,NPHY
       omuNet(i)= oCHLs(4) + i
    enddo

    do i=1,NPHY
       oGraz(i) = omuNet(NPHY) + i
    enddo

    oZ2N=oGraz(NPHY)+1
    oD2N=oZ2N+1

    oPPt  =oD2N+1

    oD_NO3=oPPt +1

    do i=1,NPHY
       oD_PHY(i)=oD_NO3+i
    enddo

    oD_ZOO=oD_PHY(NPHY)+1
    oD_DET=oD_ZOO+1
    Nout  =oD_DET

    allocate(Varout(  Nout,nlev))
    allocate(Labelout(Nout+ ow ))

    Labelout(oTemp  )='Temp '
    Labelout(oPAR   )='PAR  '
    Labelout(oAks   )='Aks  '
    Labelout(ow     )='w    '
    Labelout(oNO3+ow)='NO3  '

    do i=1,NPHY
       write(Labelout(oPHY(i)+ow), "(A3,I2)") 'PHY',i
    enddo

    Labelout(oZOO +ow)='ZOO  '
    Labelout(oDET +ow)='DET  '
    Labelout(oPHYt+ow)='PHY_T'
    Labelout(oCHLt+ow)='CHL_T'

    do i=1,NPHY
       write(Labelout(omuNet(i) + ow), "(A3,I2)") 'muN',i
       write(Labelout(oGraz(i)  + ow), "(A3,I2)") 'Gra',i
       write(Labelout(oD_PHY(i) + ow), "(A3,I2)") 'D_P',i
    enddo
  
    Labelout(oZ2N  + ow)='Z2N  '
    Labelout(oD2N  + ow)='D2N  '
    Labelout(oPPt  + ow)='NPP_T'
    Labelout(oD_NO3+ ow)='D_NO3'
    Labelout(oD_ZOO+ ow)='D_ZOO'
    Labelout(oD_DET+ ow)='D_DET'
    Labelout(oD_CHL+ ow)='D_CHL'

    do i=1,4
       write(Labelout(oCHLs(i) +ow), "(A4,I2)") 'CHLs',i
    enddo

    do i = 1, Nout+ow
       write(6,*) 'Labelout(',i,') = ',trim(Labelout(i))
    enddo
    ! Initialize parameters
    ! Indices for parameters that will be used in MCMC                 
    ! Need to have tradeoffs for maximal growth rate (mu0) and Kn
    imu0    =  1
    ialphamu=  imu0     + 1  ! mu0 increases with size
    iaI0    =  ialphamu + 1
    iQ0N    =  iaI0     + 1
    iKN     =  iQ0N     + 1
    ialphaKN=  iKN      + 1  ! KN increases with size
    itheta  =  ialphaKN + 1
    igmax   =  itheta   + 1
    ikp     =  igmax    + 1 
    iwDET   =  ikp      + 1
    irdN    =  iwDET + 1
    imz     =  irdN     + 1
    iEp     =  imz      + 1
    iEz     =  iEp      + 1
    ialphaG =  iEz      + 1
    NPar    =  ialphaG

    write(6,'(I2,1x,A20)') NPar,'parameters in total to be estimated.'
    allocate(params(NPar))
    allocate(ParamLabel(NPar))

    ParamLabel(imu0  ) = 'mu0hat '
    ParamLabel(ialphaKN) = 'alphaKN'
    ParamLabel(ialphamu) = 'alphamu'
    ParamLabel(iKN   ) = 'KN     '
    ParamLabel(igmax ) = 'gmax   '
    ParamLabel(iwDET ) = 'wDET   '
    ParamLabel(ikp   ) = 'kp     '
    ParamLabel(irdn  ) = 'rdn    '
    ParamLabel(imz   ) = 'mz     '
    ParamLabel(iEp   ) = 'Ep     '
    ParamLabel(iEz   ) = 'Ez     '
    ParamLabel(itheta) = 'theta  '
    ParamLabel(iQ0N  ) = 'Q0N    '
    ParamLabel(iaI0  ) = 'aI0    '
    ParamLabel(ialphaG)= 'alphaG '

    params(imu0   ) = 2d0
    params(ialphamu)= 0.2
    params(ialphaKN)= 0.2
    params(iaI0   ) = 0.5
    params(iKN    ) = 1d0
    params(iEp    ) = 0.4
    params(iEz    ) = 0.6
    params(igmax  ) = 1d0
    params(ikp    ) = 5d-1
    params(iQ0N   ) = 0.15
    params(itheta ) = 0.02
    params(iwDET  ) = 1D0
    params(irdn   ) = 5d-2
    params(imz    ) = 5d-2
    params(ialphaG) = 1.1

    call assign_PMU

  case(Geidersimple)
    write(6,*) 'Geider simple model selected!'
    NPHY = 1
    allocate(iPHY(NPHY))

    do i=1,NPHY
       iPHY(i)=i+iNO3
    enddo 

    iZOO=iPHY(NPHY)+1
    iDET=iZOO+1
    iCHL=iDET+1
    NVAR=iCHL 

    allocate(Vars(NVAR,nlev))

    NVsinkterms = 1+NPHY
    allocate(Windex(NVsinkterms))
    do i=1,NPHY
       Windex(i)=iPHY(i)
    enddo
    Windex(NVsinkterms)=iDET

    ! Output array matrices
    allocate(oPHY(NPHY))
    do i=1,NPHY
       oPHY(i)=i+oNO3
    enddo
    oZOO =oPHY(NPHY)+1
    oDET =oZOO+1
    oCHLt=oDET+1
    ! The above must match with i** indeces   
 
    allocate(omuNet(NPHY))
    allocate(oTheta(NPHY))
    allocate(oGraz(NPHY))
    !allocate(ow_p(NPHY))
    allocate(oD_PHY(NPHY))
    allocate(oSI(NPHY))
    allocate(oLno3(NPHY))

    do i=1,NPHY
       omuNet(i)=oCHLt + i
    enddo
    do i=1,NPHY
       oTheta(i)=omuNet(NPHY) + i
    enddo
    oSI(1)=oTheta(1)+1
    oLno3(1)=oSI(1) +1
    do i=1,NPHY
       oGraz(i)=oLno3(NPHY)+i
    enddo

    oZ2N=oGraz(NPHY)+1
    oD2N=oZ2N+1
    !do i=1,NPHY
    !   ow_p(i)=oD2N+i
    !enddo

    oPPt  =oD2N+1
    oD_NO3=oPPt +1

    do i=1,NPHY
       oD_PHY(i)=oD_NO3+i
    enddo

    oD_ZOO=oD_PHY(NPHY)+1
    oD_DET=oD_ZOO+1
    oD_CHL=oD_DET+1
    Nout  =oD_CHL

    allocate(Varout(  Nout,nlev))
    allocate(Labelout(Nout+ ow ))

    Labelout(oTemp  )='Temp '
    Labelout(oPAR   )='PAR  '
    Labelout(oAks   )='Aks  '
    Labelout(ow     )='w    '
    Labelout(oNO3+ow)='NO3  '

    do i=1,NPHY
       write(Labelout(oPHY(i)+ow), "(A3,I2)") 'PHY',i
       write(Labelout(oSI(i)+ow),  "(A3,I2)") 'SI ',i
       write(Labelout(oLno3(i)+ow),"(A3,I2)") 'LNO',i
    enddo

    Labelout(oZOO +ow)='ZOO  '
    Labelout(oDET +ow)='DET  '
    Labelout(oCHLt+ow)='CHL_T'
    do i=1,NPHY
       write(Labelout(omuNet(i) + ow), "(A3,I2)") 'muN',i
       write(Labelout(oTheta(i) + ow), "(A3,I2)") 'THE',i
       write(Labelout(oGraz(i)  + ow), "(A3,I2)") 'Gra',i
!       write(Labelout(ow_p(i)   + ow), "(A3,I2)") 'w_p',i
       write(Labelout(oD_PHY(i) + ow), "(A3,I2)") 'D_P',i
    enddo
  
    Labelout(oZ2N  + ow)='Z2N  '
    Labelout(oD2N  + ow)='D2N  '
    Labelout(oPPt  + ow)='NPP_T'
    Labelout(oD_NO3+ ow)='D_NO3'
    Labelout(oD_ZOO+ ow)='D_ZOO'
    Labelout(oD_DET+ ow)='D_DET'
    Labelout(oD_CHL+ ow)='D_CHL'

    do i = 1, Nout+ow
       write(6,*) 'Labelout(',i,') = ',Labelout(i)
    enddo
    ! Initialize parameters
    ! Indices for parameters that will be used in MCMC                 
    imu0    =  1
    iaI0    =  imu0   + 1
    iQ0N    =  iaI0   + 1
    iKN     =  iQ0N   + 1
    igmax   =  iKN    + 1
    ikp     =  igmax  + 1 
    iwDET   =  ikp    + 1
    irdN    =  iwDET  + 1
    imz     =  irdN   + 1
    iEp     =  imz    + 1
    iEz     =  iEp    + 1
    NPar    =  iEz
    write(6,'(I2,1x,A20)') NPar,'parameters in total to be estimated.'
    allocate(params(NPar))
    allocate(ParamLabel(NPar))

    ParamLabel(imu0  ) = 'mu0hat '
    ParamLabel(iaI0  ) = 'aI0    '
    ParamLabel(iKN   ) = 'KN     '
    ParamLabel(iQ0N  ) = 'Q0N    '
    ParamLabel(iwDET ) = 'wDET   '
    ParamLabel(igmax ) = 'gmax   '
    ParamLabel(ikp   ) = 'kp     '
    ParamLabel(irdn  ) = 'rdn    '
    ParamLabel(imz   ) = 'mz     '
    ParamLabel(iEp   ) = 'Ep     '
    ParamLabel(iEz   ) = 'Ez     '

    params(imu0   ) = 2d0
    params(iaI0   ) = 0.5
    params(iKN    ) = 1d0
    params(iEp    ) = 0.4
    params(iEz    ) = 0.6
    params(igmax  ) = 1d0
    params(ikp    ) = 5d-1
    params(iQ0N   ) = 0.15
    params(iwDET  ) = 1D0
    params(irdn   ) = 5d-2
    params(imz    ) = 5d-2

   CASE(EFTsimple)

    write(6,*) 'FlexEFT simple model (only ONE PHYTO size class) selected!'
    NPHY = 1
    allocate(iPHY(NPHY))

    do i=1,NPHY
       iPHY(i)=i+iNO3
    enddo 

    iZOO=iPHY(NPHY)+1
    iDET=iZOO+1
    NVAR=iDET 

    allocate(Vars(NVAR,nlev))

    NVsinkterms = 1+NPHY
    allocate(Windex(NVsinkterms))
    do i=1,NPHY
       Windex(i)=iPHY(i)
    enddo
    Windex(NVsinkterms)=iDET

    ! Output array matrices
    allocate(oPHY(NPHY))
    do i=1,NPHY
       oPHY(i)=i+oNO3
    enddo
    oZOO =oPHY(NPHY)+1
    oDET =oZOO+1
    oCHLt=oDET+1

    ! The above must match with i** indeces   
    allocate(omuNet(NPHY))
    allocate(oGraz(NPHY))
    allocate(oTheta(NPHY))
    allocate(oQN(NPHY))
    allocate(oSI(NPHY))
    allocate(oLno3(NPHY))
    !allocate(ow_p(NPHY))
    allocate(oD_PHY(NPHY))

    do i=1,NPHY
       omuNet(i)=oCHLt + i
    enddo

    do i=1,NPHY
       oGraz(i)=omuNet(NPHY)+i
    enddo

    oZ2N=oGraz(NPHY)+1
    oD2N=oZ2N+1
    !do i=1,NPHY
    !   ow_p(i)=oD2N+i
    !enddo
    do i=1,NPHY
       oSI(i)=oD2N+i
    enddo

    do i=1,NPHY
       oLno3(i)=oSI(NPHY)+i
    enddo
    do i=1,NPHY
       oTheta(i)=oLno3(NPHY)+i
    enddo

    do i=1,NPHY
       oQN(i)=oTheta(NPHY)+i
    enddo

    oPPt  =oQN(NPHY)+1
    oD_NO3=oPPt +1

    do i=1,NPHY
       oD_PHY(i)=oD_NO3+i
    enddo

    oD_ZOO=oD_PHY(NPHY)+1
    oD_DET=oD_ZOO+1
    Nout  =oD_DET
    write(6,'(A8,1x, I3, A25)') 'Totally',Nout,'Variables for diagosis!'
    allocate(Varout(  Nout,nlev))
    allocate(Labelout(Nout+ ow ))

    Labelout(oTemp  )='Temp '
    Labelout(oPAR   )='PAR  '
    Labelout(oAks   )='Aks  '
    Labelout(ow     )='w    '
    Labelout(oNO3+ow)='NO3  '

    do i=1,NPHY
       write(Labelout(oPHY(i)+ow), "(A3,I2)") 'PHY',i
    enddo

    Labelout(oZOO +ow)='ZOO  '
    Labelout(oDET +ow)='DET  '
    Labelout(oCHLt+ow)='CHL_T'
    do i=1,NPHY
       write(Labelout(omuNet(i) + ow), "(A3,I2)") 'muN',i
       write(Labelout(oGraz(i)  + ow), "(A3,I2)") 'Gra',i
     !  write(Labelout(ow_p(i)   + ow), "(A3,I2)") 'w_p',i
       write(Labelout(oD_PHY(i) + ow), "(A3,I2)") 'D_P',i
       write(Labelout(oTheta(i) + ow), "(A3,I2)") 'THE',i
       write(Labelout(oQN(i)    + ow), "(A3,I2)") 'QN ',i
       write(Labelout(oSI(i)    + ow), "(A3,I2)") 'SI ',i
       write(Labelout(oLno3(i)  + ow), "(A3,I2)") 'LNO',i
    enddo
  
    Labelout(oZ2N  + ow)='Z2N  '
    Labelout(oD2N  + ow)='D2N  '
    Labelout(oPPt  + ow)='NPP_T'
    Labelout(oD_NO3+ ow)='D_NO3'
    Labelout(oD_ZOO+ ow)='D_ZOO'
    Labelout(oD_DET+ ow)='D_DET'

    do i = 1, Nout+ow
       write(6,'(A9,I2,A4,A5)') 'Labelout(',i,') = ',Labelout(i)
    enddo
    ! Initialize parameters
    ! Indices for parameters that will be used in MCMC                 
    !imu0    =  1
    !iaI0    =  imu0    + 1
    !iV0N    =  iaI0    + 1
    iaI0 = 1

    if (nutrient_uptake .eq. 1) then
       iKN  =  iaI0    + 1
       iQ0N =  iKN     + 1
    else if(nutrient_uptake .eq. 2) then 
       iA0N =  iaI0    + 1
       iQ0N =  iA0N    + 1
    else
       write(6,*) 'Option of nutrient uptake incorrect! Quit!'
       stop
    endif

    igmax   =  iQ0N    + 1
    ikp     =  igmax   + 1 
    iwDET   =  ikp     + 1
    irdN    =  iwDET   + 1
    imz     =  irdN    + 1
    iEp     =  imz     + 1
    iEz     =  iEp     + 1
    izetaN  =  iEz     + 1
    izetaChl=  izetaN  + 1
    NPar    =  izetaChl

    allocate(params(NPar))
    write(6,'(A23,I3,A30)') 'There are totally ',NPar,'parameters to be fitted.'

    allocate(ParamLabel(NPar))
    !ParamLabel(imu0 ) = 'mu0hat '
    ParamLabel(iaI0 )   = 'aI0    '
    ParamLabel(iEp  )   = 'Ep     '
    ParamLabel(iEz  )   = 'Ez     '
    ParamLabel(izetaN)  = 'zetaN  '
    ParamLabel(izetaChl)= 'zetaChl'

    if (nutrient_uptake .eq. 1) then
        ParamLabel(iKN  ) = 'KN     '
    else if(nutrient_uptake .eq. 2) then
        ParamLabel(iA0N ) = 'A0N    '
    else
        write(6,*) 'Option of nutrient uptake incorrect! Quit!'
        stop
    endif

    !ParamLabel(iV0N )  = 'V0N    '
    ParamLabel(iQ0N )  = 'Q0N    '
    ParamLabel(iwDET)  = 'wDET   '
    ParamLabel(igmax)  = 'gmax   '
    ParamLabel(ikp  )  = 'kp     '
    ParamLabel(irdn )  = 'rdn    '
    ParamLabel(imz  )  = 'mz     '
    !ParamLabel(iPHYini)= 'PHYini '
    !ParamLabel(iZOOini)= 'ZOOini '
    !ParamLabel(iDETini)= 'DETini '

    !params(imu0)       = 5d0
    params(iaI0)       = 1d0
    params(iEp )       = 0.4
    params(iEz )       = 0.6
    params(izetaN)     = 0.6
    params(izetaChl)   = 0.8
    !params(iV0N)       = 5d0
    params(igmax)      = 1d0
    params(ikp)        = 1d0
    if (nutrient_uptake .eq. 1) then
       params(iKN)     = 0.1
    elseif (nutrient_uptake .eq. 2) then
       params(iA0N )   = 1d0
    else
       write(6,*) 'Option of nutrient uptake incorrect! Quit!'
       stop
    endif
    params(iQ0N )   = 0.05
    params(iwDET)   = 1d0
    params(irdn )   = 0.05
    params(imz  )   = 0.05

  case(EFTdiscrete)

    WRITE(6,*) 'EFTdiscrete model (20 PHY size classes) selected!'
    NPHY=20
    ! Assign variable array indices 
    ! and not including CHL in the variable lists 
    allocate(iPHY(NPHY))

    do i=1,NPHY
       iPHY(i)=i+iNO3
    enddo 

    iZOO=iPHY(NPHY)+1
    iDET=iZOO+1
    NVAR=iDET 

    allocate(Vars(NVAR,nlev),  STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

    NVsinkterms = 1+NPHY
    allocate(Windex(NPHY+1),  STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

    do i=1,NPHY
       Windex(i)=iPHY(i)
    enddo
    Windex(NVsinkterms)=iDET

    ! Output array matrices
    allocate(oPHY(NPHY))
    do i=1,NPHY
       oPHY(i)=i+oNO3
    enddo
    oZOO=oPHY(NPHY)+1
    oDET=oZOO+1
    oFER=oDET+1

    allocate(oTheta(NPHY))
    allocate(oQN(NPHY))
    allocate(omuNet(NPHY))
    allocate(oGraz(NPHY))
    !allocate(ow_p(NPHY))
    allocate(oSI(NPHY))
    allocate(oLno3(NPHY))
    allocate(oD_PHY(NPHY))
    allocate(oTheHat(NPHY))

    do i=1,NPHY
       oTheta(i)=oFER+i
    enddo

    do i=1,NPHY
       oQN(i)=oTheta(NPHY)+i
    enddo

    do i=1,NPHY
       omuNet(i)=oQN(NPHY)+i
    enddo

    do i=1,NPHY
       oGraz(i)=omuNet(NPHY)+i
    enddo

    oZ2N=oGraz(NPHY)+1

    oD2N=oZ2N+1
    !do i=1,NPHY
    !   ow_p(i)=oD2N+i
    !enddo

    do i=1,NPHY
       oSI(i)=oD2N+i
    enddo

    do i=1,NPHY
       oLno3(i)=oSI(NPHY)+i
    enddo

    do i=1,NPHY
       oTheHat(i)=oLno3(NPHY)+i
    enddo

    oPHYt =oTheHat(NPHY)+1
    oCHLt =oPHYt+1

    do i=1,4
       oCHLs(i)=oCHLt+i
    enddo
    oPPt  =oCHLs(4)+1
    oPMU  =oPPt+1
    oVAR  =oPMU+1
    oD_NO3=oVAR+1

    do i=1,NPHY
       oD_PHY(i)=oD_NO3+i
    enddo

    oD_ZOO=oD_PHY(NPHY)+1
    oD_DET=oD_ZOO+1
    Nout  =oD_DET

    allocate(Varout(    Nout,nlev),  STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Problem in allocating Varout ***"

    allocate(Labelout(  Nout+ ow ),  STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Problem in allocating Labelout ***"

    Labelout(oTemp  )='Temp '
    Labelout(oPAR   )='PAR  '
    Labelout(oAks   )='Aks  '
    Labelout(ow     )='w    '
    Labelout(oNO3+ow)='NO3  '
    do i=1,NPHY
       write(Labelout(oPHY(i)+ow), "(A3,I2)") 'PHY',i
    enddo
    Labelout(oZOO+ow)='ZOO  '
    Labelout(oDET+ow)='DET  '
    Labelout(oFER+ow)='FER  '
    do i=1,NPHY
       write(Labelout(oTheta(i)  + ow), "(A3,I2)") 'THE',i
       write(Labelout(oQN(i)     + ow), "(A3,I2)") 'QN ',i
       write(Labelout(omuNet(i)  + ow), "(A3,I2)") 'muN',i
       write(Labelout(oGraz(i)   + ow), "(A3,I2)") 'Gra',i
!       write(Labelout(ow_p(i)    + ow), "(A3,I2)") 'w_p',i
       write(Labelout(oSI(i)     + ow), "(A3,I2)") 'SI ',i
       write(Labelout(oLno3(i)   + ow), "(A3,I2)") 'LNO',i
       write(Labelout(oTheHat(i) + ow), "(A3,I2)") 'THA',i
       write(Labelout(oD_PHY(i)  + ow), "(A3,I2)") 'D_P',i
    enddo
    do i=1,4
       write(Labelout(oCHLs(i)   + ow), "(A4,I2)") 'CHLs',i
    enddo

    Labelout(oZ2N  + ow)='Z2N  '
    Labelout(oD2N  + ow)='D2N  '
    Labelout(oPHYt + ow)='PHY_T'
    Labelout(oCHLt + ow)='CHL_T'
    Labelout(oPPt  + ow)='NPP_T'
    Labelout(oPMU  + ow)='PMU  '
    Labelout(oVAR  + ow)='VAR  '
    Labelout(oD_NO3+ ow)='D_NO3'
    Labelout(oD_ZOO+ ow)='D_ZOO'
    Labelout(oD_DET+ ow)='D_DET'
    do i = 1, Nout+ow
       write(6,'(A9,I4,A4,A5)') 'Labelout(',i,') = ',Labelout(i)
    enddo

    ! Initialize parameters
    ! Indices for parameters that will be used in MCMC                 
    !imu0    =  1
    !iaI0    =  imu0    + 1
    !ialphamu=  1
    !ibetamu =  ialphamu + 1
    iaI0    =  1
    ialphaI =  iaI0     + 1
    !iV0N    =  ialphaI + 1
    !ialphaV =  iV0N    + 1
    select case(nutrient_uptake)
    case(1)
      iKN     =  ialphaI + 1
      ialphaK =  iKN     + 1
      iQ0N    =  ialphaK + 1
    case(2)
      iA0N    =  ialphaI + 1
      ialphaA =  iA0N    + 1
      iQ0N    =  ialphaA + 1
    case default
      write(6,*) 'Option of nutrient uptake incorrect! Quit!'
      stop
    end select

    ialphaQ =  iQ0N    + 1
    igmax   =  ialphaQ + 1
    ikp     =  igmax   + 1 
    ialphaG =  ikp     + 1
    iwDET   =  ialphaG + 1
    irdN    =  iwDET   + 1
    imz     =  irdN    + 1
    iEp     =  imz     + 1
    iEz     =  iEp     + 1
    !izetaN  =  iEz     + 1
    !izetaChl=  izetaN  + 1
    NPar    =  iEz

    allocate(params(NPar))
    allocate(ParamLabel(NPar))

    !ParamLabel(ialphamu) = 'alphamu'
    !ParamLabel(ibetamu ) = 'betamu'
    ParamLabel(iaI0    ) = 'aI0    '
    ParamLabel(ialphaI ) = 'alphaI '
    if (nutrient_uptake .eq. 1) then
       ParamLabel(iKN    )  = 'K0N    '
       ParamLabel(ialphaK ) = 'alphaK '
    elseif (nutrient_uptake .eq. 2) then
       ParamLabel(iA0N    ) = 'A0N    '
       ParamLabel(ialphaA ) = 'alphaA '
    else
       print *, 'Option of nutrient uptake incorrect! Quit!'
       stop
    endif
!    ParamLabel(iV0N    ) = 'V0N    '
!    ParamLabel(ialphaV ) = 'alphaV '
    ParamLabel(iQ0N    ) = 'Q0N    '
    ParamLabel(ialphaQ ) = 'alphaQ '
    ParamLabel(ialphaG ) = 'alphaG '
    ParamLabel(iwDET   ) = 'wDET   '
    ParamLabel(igmax   ) = 'gmax   '
    ParamLabel(ikp     ) = 'kp     '
    ParamLabel(irdn    ) = 'rdn    '
    ParamLabel(imz     ) = 'mz     '
    ParamLabel(iEp     ) = 'Ep '
    ParamLabel(iEz     ) = 'Ez '
    !ParamLabel(izetaN  ) = 'zetaN '
    !ParamLabel(izetaChl) = 'zetaChl '

    !params(ialphamu)=0d0
    !params(ibetamu) =0d0
    params(iaI0)    =0.5
  !  params(iV0N)    =5d0
    params(igmax)   =1d0
    params(ikp     )=5d-1
    params(ialphaI) =-0.13
    if (nutrient_uptake .eq. 1) then
       params(iKN)     =1d0
       params(ialphaK) =0.27
    elseif (nutrient_uptake .eq. 2) then
       params(iA0N   ) =1D0
       params(ialphaA) =-0.3d0
    endif
  !  params(ialphaV) =1D-6      ! To avoid zero
    params(ialphaG) =1.1d0
    params(iQ0N   ) =4d-2
    params(ialphaQ) =-0.17d0
    params(iwDET  ) =1D0
    params(irdn   ) =5d-2
    params(imz    ) =5d-2
    params(iEp    ) =.4
    params(iEz    ) =.6
    !params(izetaN ) =.6
    !params(izetaChl)=.8

    call assign_PMU

    CASE(EFTcont)

    NPHY=1
    ! Assign variable array indices 
    ! and not including CHL in the state variable lists 
    allocate(iPHY(NPHY))

    do i=1,NPHY
       iPHY(i)=i+iNO3
    enddo 

    iZOO = iPHY(NPHY)+ 1
    iDET = iZOO      + 1
    iPMU = iDET      + 1
    iVAR = iPMU      + 1
    NVAR = iVAR 

    allocate(Vars(NVAR,nlev))
    NVsinkterms = NPHY+1

    allocate(Windex(NVsinkterms))
    do i=1,NPHY
       Windex(i)=iPHY(i)
    enddo

    Windex(NVsinkterms)=iDET

    ! Output array matrices
    allocate(oPHY(NPHY))
    do i=1,NPHY
       oPHY(i)=i+oNO3
    enddo

    oZOO=oPHY(NPHY)+1
    oDET=oZOO+1
    oPMU=oDET+1
    oVAR=oPMU+1
    oFER=oVAR+1

    allocate(oTheta(NPHY))
    allocate(oQN(NPHY))
    allocate(omuNet(NPHY))
    allocate(oGraz(NPHY))
   ! allocate(ow_p(NPHY))
    allocate(oD_PHY(NPHY))

    do i=1,NPHY
       oTheta(i)=oFER+i
    enddo

    do i=1,NPHY
       oQN(i)=oTheta(NPHY)+i
    enddo

    do i=1,NPHY
       omuNet(i)=oQN(NPHY)+i
    enddo

    do i=1,NPHY
       oGraz(i)=omuNet(NPHY)+i
    enddo

    oCHLt=oGraz(NPHY)+1
    do i=1,4
       oCHLs(i)=oCHLt+i
    enddo
    oPPt=oCHLs(4)+1
    oZ2N=oPPt+1
    oD2N=oZ2N+1
    !do i=1,NPHY
    !   ow_p(i)=oD2N+i
    !enddo

    oD_NO3=oD2N+1
    do i=1,NPHY
       oD_PHY(i)=oD_NO3+1
    enddo
    oD_ZOO=oD_PHY(NPHY)+1
    oD_DET=oD_ZOO+1
    oD_PMU=oD_DET+1
    oD_VAR=oD_PMU+1
    odmudl=oD_VAR+1
    od2mu =odmudl+1
    od2gdl=od2mu +1
    Nout  =od2gdl

    allocate(Varout(Nout,nlev), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Problem in allocating Varout ***"
    allocate(Labelout(Nout+ow))

    Labelout(oTemp ) ='Temp '
    Labelout(oPAR  ) ='PAR  '
    Labelout(oAks  ) ='Aks  '
    Labelout(ow    ) ='w    '
    Labelout(oNO3+ow)='NO3  '
    do i=1,NPHY
       write(Labelout(oPHY(i)+ow), "(A3,I2)") 'PHY',i
    enddo
    Labelout(oZOO+ow)='ZOO  '
    Labelout(oDET+ow)='DET  '
    Labelout(oPMU+ow)='PMU  '
    Labelout(oVAR+ow)='VAR  '
    Labelout(oFER+ow)='FER  '
    do i=1,NPHY
       write(Labelout(oTheta(i)+ow), "(A3,I2)") 'The',i
       write(Labelout(oQN(i)   +ow), "(A3,I2)") 'QN ',i
       write(Labelout(omuNet(i)+ow), "(A3,I2)") 'muN',i
       write(Labelout(oGraz(i) +ow), "(A3,I2)") 'Gra',i
!       write(Labelout(ow_p(i)  +ow), "(A3,I2)") 'w_p',i
       write(Labelout(oD_PHY(i)+ow), "(A3,I2)") 'D_P',i
    enddo
    do i=1,4
       write(Labelout(oCHLs(i) +ow), "(A4,I2)") 'CHLs',i
    enddo

    Labelout(oCHLt +ow)='CHL_T'
    Labelout(oPPt  +ow)='NPP_T'
    Labelout(oZ2N  +ow)='Z2N  '
    Labelout(oD2N  +ow)='D2N  '
    Labelout(oD_NO3+ow)='D_NO3'
    Labelout(oD_ZOO+ow)='D_ZOO'
    Labelout(oD_DET+ow)='D_DET'
    Labelout(oD_PMU+ow)='D_PMU'
    Labelout(oD_VAR+ow)='D_VAR'
    Labelout(odmudl+ow)='dmudl'
    Labelout(od2mu +ow)='d2mu '
    Labelout(od2gdl+ow)='d2gdl'

    do i = 1, Nout+ow
       write(6,*) 'Labelout(',i,') = ',Labelout(i)
    enddo
    ! Initialize parameters
    !===================================================
    ! Define indices for parameters that will be used in MCMC                 
    !===================================================
!    ialphamu=  1
!    ibetamu =  ialphamu+ 1
    iaI0    =  1
    ialphaI =  iaI0    + 1
    !iV0N    =  ialphaI + 1
    !ialphaV =  iV0N    + 1

    select case(nutrient_uptake)
    case(1)
      iKN     =  ialphaI + 1
      ialphaK =  iKN     + 1
      iQ0N    =  ialphaK + 1
    case(2)
      iA0N    =  ialphaI + 1
      ialphaA =  iA0N    + 1
      iQ0N    =  ialphaA + 1
    case default
      write(6,*) 'Option of nutrient uptake incorrect! Quit!'
      stop
    end select

    ialphaQ =  iQ0N    + 1
    igmax   =  ialphaQ + 1
    ikp     =  igmax   + 1 
    ialphaG =  ikp     + 1
    iwDET   =  ialphaG + 1
    irdN    =  iwDET   + 1
    imz     =  irdN    + 1
    iPenfac =  imz     + 1
    ithetamin= iPenfac + 1
    iQNmin  =  ithetamin+1
    iEp     =  iQNmin   + 1
    iEz     =  iEp      + 1
    !izetaN  =  iEz      + 1
    !izetaChl=  izetaN   + 1
    NPar    =  iEz

    allocate(params(NPar))
    allocate(ParamLabel(NPar))
    !===================================================
    ! Give parameter labels
    !===================================================

!    ParamLabel(ialphamu) = 'alphamu'
!    ParamLabel(ibetamu) = 'betamu'
    ParamLabel(iaI0    ) = 'aI0    '
    ParamLabel(ialphaI ) = 'alphaI '
    if (nutrient_uptake .eq. 1) then
        ParamLabel(iKN     ) = 'K0N    '
        ParamLabel(ialphaK ) = 'alphaK '
    elseif (nutrient_uptake .eq. 2) then
        ParamLabel(iA0N    ) = 'A0N    '
        ParamLabel(ialphaA ) = 'alphaA '
    endif
!    ParamLabel(iV0N    ) = 'V0N    '
!    ParamLabel(ialphaV ) = 'alphaV '
    ParamLabel(iQ0N    ) = 'Q0N    '
    ParamLabel(ialphaQ ) = 'alphaQ '
    ParamLabel(ialphaG ) = 'alphaG '
    ParamLabel(iwDET   ) = 'wDET   '
    ParamLabel(igmax   ) = 'gmax   '
    ParamLabel(ikp     ) = 'kp     '
    ParamLabel(irdn    ) = 'rdn    '
    ParamLabel(imz     ) = 'mz     '
    ParamLabel(iPenfac ) = 'Penfac '
    ParamLabel(ithetamin)= 'thetmin'
    ParamLabel(iQNmin )  = 'QNmin '
    ParamLabel(iEp     ) = 'Ep '
    ParamLabel(iEz     ) = 'Ez '
   ! ParamLabel(izetaN  ) = 'zetaN '
   ! ParamLabel(izetaChl) = 'zetaChl '

!    params(ialphamu)=0d0
!    params(ibetamu)=0d0
    params(iaI0    )=0.25d0
!    params(iV0N    )=5d0
    params(igmax   )=1d0
    params(ialphaI )=-0.13
    if (nutrient_uptake .eq. 1) then
       params(iKN     )=1d0
       params(ialphaK )=0.27
    elseif (nutrient_uptake .eq. 2) then
       params(iA0N    )=4D1
       params(ialphaA )=-3d-1
    else
       write(6,*) 'Option of nutrient uptake incorrect!'
       stop
    endif
!    params(ialphaV)=2d-1
    params(ialphaG)=1.1d0
    params(iQ0N   )=1d-1
    params(ialphaQ)=-0.17d0
    params(iwDET  )=1D0
    params(ikp    )=1D0
    params(irdn   )=5d-2
    params(imz    )=5d-2
    params(iPenfac)=1d0
    params(ithetamin)=0.05
    params(iQNmin)=1d-2
    params(iEp)   =0.4
    params(iEz   )=0.6
   ! params(izetaN)=0.6
   ! params(izetaChl)=0.8

  CASE DEFAULT

    write(6,*) 'Error: Incorrect option for biological models!'
    stop

ENDSELECT
end subroutine choose_model
