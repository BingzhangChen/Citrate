SUBROUTINE FLEXEFT_CONT
use bio_MOD
implicit none
!INPUT PARAMETERS:
real :: tC,par_
!LOCAL VARIABLES of phytoplankton:
real :: PMUPHY,VARPHY,NO31,PHY1,ZOO1
real :: PMUPHY1,VARPHY1  !State variables
real :: larg,Ep,Ez,zetaN, zetaChl 
real :: PMU,VAR,PMU1,VAR1,X,B,dXdl,dBdl,KFe,Fe,dmu0hatdl, d2mu0hatdl2
real :: V0hat,Kn,Lmin,A0N,A0N0, A0hat,dA0hatdl,d2A0hatdl2,fA
real :: VNhat,dVNhatdl,d2VNhatdl2 ! Nutrient uptake variables
real :: mu0hat,muIhat,mu0hatSI,dmu0hatSIdl,d2mu0hatSIdl2
real :: dmuIhatdl,d2muIhatdl2! Growth rate variables
real :: fV,dfVdl,d2fVdl2  !fV
real :: ZINT,dZINdl,d2ZINdl2 !ZINT
real :: aI,SI,dSIdl,d2SIdl2         !Light dependent component
real :: RMchl,Theta,ThetaHat,dThetaHatdl,d2ThetaHatdl2 !Chl related variables
real :: QN,Qs,dQNdl,d2QNdl2  ! cell quota related variables
real :: larg1,dmu0hat_aIdl,d2mu0hat_aIdl2,dlargdl
real :: d2largdl2,W_Y,dWYYdl,daI_mu0hatdl,d2aI_mu0hatdl2
real :: dV0hatdl,d2V0hatdl2,muNet,dmuNetdl,d2muNetdl2
real :: ThetaAvg,QNavg, dgdlbar,d2gdl2bar, alphaA
real :: Q0N,alphaQ,aI0,alphaI,Lref,Penfac
real :: alphaK, K0N
real :: NO3, PHY, ZOO, DET, DET1
!Declarations of zooplankton:
integer :: k,j

real :: pp_PN, pp_NP, PPpn
real :: pCHL(4) = 0d0
logical, parameter :: do_IRON = .false.
real,    external  :: pnorm, WAPR
!-----------------------------------------------------------------------
Ep     = params(iEp)
alphamu= 0d0
betamu = 0d0
Ez     = params(iEz)
zetaN  = 0.6
zetaChl= 0.8
!zetaN  = params(izetaN)
!zetaChl= params(izetaChl)

DO k = nlev, 1, -1   
   ! Retrieve current (local) state variable values.
   tC     = Temp(k)
   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif

   NO3    = max(Vars(iNO3,k)   , 1D-10)
   PHY    = max(Vars(iPHY(1),k), 1D-10)
   PHY    = min(PHY,             1D+10)
   ZOO    = max(Vars(iZOO,k)   , 1D-10)
   DET    = max(Vars(iDET,k)   , 1D-10)
   PMUPHY = max(Vars(iPMU,k)   , 1D-10)
   VARPHY = max(Vars(iVAR,k)   , 1D-10)
!All rates have been multiplied by dtdays to get the real rate correponding to the actual time step
   tf_p   = TEMPBOL(Ep,tC)
   PMU    = min(PMUPHY/PHY, PMUmax)
   PMU    = PMU - log(1d1)  ! Correct PMU to the real one
   VAR    = min(VARPHY/PHY, VARmax)

! Fe related:
   if (do_IRON) then
!      Fe   = Vars(iFER,k)

      !Dissolved Fe concentration
      Fe   = max(Fe,Femin)  

      !The half saturation constant of Fe at average size
      KFe  = ScaleTrait(PMU, K0Fe,alphaFe)
   endif
       
   alphaQ  = params(ialphaQ)
   Q0N     = params(iQ0N)
   Qs      = ScaleTrait(PMU, Q0N, alphaQ)/2d0
   mu0hat  = dtdays*tf_p*mu0*exp(alphamu*PMU + betamu*PMU**2)
   X       = 0d0
   if (do_IRON) then
     mu0hat  = mu0hat * Fe/(Fe + KFe)
     X       = alphaFe*KFe/(Fe + KFe) 
   endif
   
   dmu0hatdl  = mu0hat*(alphamu + 2D0*betamu*PMU - X)
   d2mu0hatdl2= dmu0hatdl*(alphamu + 2D0*betamu*PMU-X) + mu0hat*2D0*betamu
   
   if (do_IRON) d2mu0hatdl2 = d2mu0hatdl2 - mu0hat*alphaFe*Fe*X/(Fe+KFe)
   
   ! Iron limits nitrogen uptake

   !V0N    = params(iV0N)
   !alphaV = params(ialphaV)
   V0hat  = ScaleTrait(PMU, dtdays*tf_p*V0N, alphaV)
     
   if (do_IRON) V0hat = V0hat*Fe/(Fe + KFe)
     
   dV0hatdl   = V0hat*(alphaV- X)
   d2V0hatdl2 = dV0hatdl*(alphaV - X)
     
   if (do_Iron) d2V0hatdl2 = d2V0hatdl2 - V0hat*alphaFe*Fe*X/(Fe+KFe) 
  
   ! Initial slope of P-I curve
   aI0    = params(iaI0)
   alphaI = params(ialphaI)
   aI     = ScaleTrait(PMU, dtdays*aI0, alphaI)

   ! Cost of photosynthesis
   RMchl  = tf_p*RMchl0*dtdays

   ! Threshold irradiance and RMchl is set temperature dependent
   I_zero = zetaChl*RMchl/aI  
  
   !Define VNhat: the saturation function of ambient nutrient concentration
   SELECTCASE(nutrient_uptake)  
   ! case 1: Classic Michaelis Menton 
   case(1)
   ! Potential maximal nutrient-limited uptake
   ! Half-saturation constant of nitrate uptake
    alphaK = params(ialphaK)
    K0N    = params(iKN)
    Kn     = ScaleTrait(PMU,K0N, alphaK) 
    VNhat  = V0hat*NO3/(NO3 + Kn)
   
    dVNhatdl  = -VNhat*alphaK*Kn/(NO3+Kn) + NO3/(NO3 + Kn)*dV0hatdl
   
     d2VNhatdl2 = -alphaK*(VNhat * alphaK * NO3/(NO3+Kn)**2*Kn          &
   + Kn/(NO3 + Kn)*dVNhatdl)+ NO3/(NO3+ Kn)*d2V0hatdl2                  &
   - dV0hatdl*NO3/(NO3 + Kn)**2*alphaK*Kn
  
     ! case 2: optimal uptake based on Pahlow (2005) and Smith et al. (2009)
   case(2)

     Lref  = 0.5       !params(iLref) 
     Penfac= params(iPenfac)
     alphaA= params(ialphaA)

     Lmin  = log(Lref**3/6d0*pi) + log(1d0+Penfac)/( Penfac*alphaA) 

     A0N0  = params(iA0N)
     A0N   = dtdays*tf_p*A0N0
     A0hat = PenAff(PMU, alphaA, Penfac,Lmin)*ScaleTrait(PMU,A0N,alphaA)
     A0hat = max(A0hat,1D-9*A0N)  ! Maintain positivity
  
     dA0hatdl = alphaA*A0hat-A0N*exp(PMU*alphaA)*Penfac*alphaA           &
             * exp(Penfac*alphaA *(PMU-Lmin))
    
     d2A0hatdl2 = alphaA*dA0hatdl                                        &
             - Penfac*alphaA*exp(alphaA*((1d0+Penfac)*PMU-Penfac*Lmin))  &
             * (dA0hatdl + A0hat*alphaA*(1d0+Penfac))  
      
             !Define fA
      fA    = min(1D0/( 1D0 + sqrt(A0hat * NO3/V0hat) ), 0.9999D0 )
      VNhat = (1D0-fA)*V0hat*fA*A0hat*NO3/((1D0-fA)*V0hat+fA*A0hat*NO3) 
      VNhat = max(V0hat,1D-9*V0hat)
      
      !X: temporary variable
      X    = V0hat/A0hat + 2d0*sqrt(V0hat*NO3/A0hat) + NO3       
    
      !B: d(V0/A0)dl
      B    = dV0hatdl/A0hat - V0hat/A0hat**2*dA0hatdl
    
      dXdl = B*(1d0 + sqrt(NO3*A0hat/V0hat))
    
      dBdl = d2V0hatdl2/A0hat - dV0hatdl*dA0hatdl/A0hat**2              &
       - (V0hat/A0hat**2*d2A0hatdl2                                     &
       +  dA0hatdl*(dV0hatdl/A0hat**2 - 2d0*V0hat*dA0hatdl/A0hat**3))
      
      dVNhatdl = NO3*(dV0hatdl/X-V0hat/X**2*B*(1d0+ sqrt(NO3*A0hat/V0hat) ) )
    
      d2VNhatdl2 = NO3*(d2V0hatdl2/X - dV0hatdl*dXdl/X**2               &
       - (V0hat/X**2*B*(-sqrt(NO3)/2d0*(A0hat/V0hat)                    &
       * sqrt(A0hat/V0hat) * B)                                         &
       + B*(1d0 + sqrt(NO3*A0hat/V0hat) )                               &
       * (dV0hatdl/X**2 - 2d0*V0hat*dXdl/X**3)                          &
       + V0hat/X**2*(1d0+sqrt(NO3*A0hat/V0hat))*dBdl))
    
   CASE DEFAULT
      write(6,*) 'Error: Incorrect option for nutrient uptake!'
      STOP
   ENDSELECT  

! Calculate thetahat (optimal g Chl/mol C for the chloroplast under nutrient replete conditions)
! Only calculate within the euphotic zone, otherwise many numerical problems.

      if( par_ .gt. I_zero ) then
        
        larg1 = exp(1d0 + min(aI*par_/(mu0hat*zetaChl),6d2))
        larg  = (1d0 + RMchl/mu0hat)*larg1   
        dmu0hat_aIdl   = (dmu0hatdl - alphaI*mu0hat)/aI
        d2mu0hat_aIdl2 = d2mu0hatdl2/aI -alphaI/aI*dmu0hatdl-aI*dmu0hat_aIdl
        daI_mu0hatdl   = -(aI/mu0hat)**2*dmu0hat_aIdl

        d2aI_mu0hatdl2 = -((aI/mu0hat)**2*d2mu0hat_aIdl2                 &
       - 2d0/((mu0hat/aI)*(mu0hat/aI)*(mu0hat/aI))                       &
       * (dmu0hat_aIdl)**2    )
   
        dlargdl = -RMchl*larg1/(mu0hat**2)*dmu0hatdl                     &
       + (1d0+RMchl/mu0hat)*larg1 * par_/zetaChl*daI_mu0hatdl
        
        d2largdl2 = -RMchl*(larg1/(mu0hat*mu0hat)*d2mu0hatdl2            &
       + larg1*par_/zetaChl*daI_mu0hatdl/(mu0hat*mu0hat)*dmu0hatdl       &
       + larg1*dmu0hatdl*(-2d0/(mu0hat*mu0hat*mu0hat)*dmu0hatdl))        &
       + par_/zetaChl*((1.+RMchl/mu0hat)*larg1*d2aI_mu0hatdl2            &
       + (1.+RMchl/mu0hat)*larg1*par_/zetaChl*daI_mu0hatdl*daI_mu0hatdl  &
       + RMchl*(-dmu0hatdl/(mu0hat*mu0hat))*larg1*daI_mu0hatdl)
   
       W_Y      = WAPR(larg)
       ThetaHat = 1d0/zetaChl + (1d0- W_Y)*mu0hat/(aI * par_)
   
       dThetaHatdl = 1d0/par_*(-W_Y/larg/(1.+W_Y)*dlargdl*mu0hat/aI      &
     +  (1d0-W_Y)*dmu0hat_aIdl)
       
       dWYYdl = dlargdl*(-W_Y**2/(larg*larg)/((1d0+W_Y)*(1.+W_Y)*(1.+W_Y))  &
       -  W_Y/larg**2/(1d0+W_Y) + W_Y/(larg*(1d0+W_Y))**2)
   
       d2ThetaHatdl2 = 1d0/par_*(-(W_Y/larg/(1d0+W_Y)*dlargdl*dmu0hat_aIdl  &
       +  W_Y/larg/(1d0+W_Y)*d2largdl2*mu0hat/aI   &
       +  dWYYdl*dlargdl*mu0hat/aI)               &
       -  W_Y/larg/(1d0+W_Y)*dlargdl * dmu0hat_aIdl &
       +  (1d0-W_Y)*d2mu0hat_aIdl2)
   
        SI = 1d0 - max(exp(-aI*par_*ThetaHat/mu0hat),0d0)
   
        dSIdl = ( (alphaI- dmu0hatdl/mu0hat)   &
        * ThetaHat + dThetaHatdl) * (1d0-SI)*aI*par_/mu0hat    !confirmed
   
       d2SIdl2 = par_*(- dSIdl*aI/mu0hat*(ThetaHat*alphaI     &
      - ThetaHat/mu0hat*dmu0hatdl + dThetaHatdl) + (1d0-SI)   &
      * (ThetaHat*alphaI- ThetaHat/mu0hat*dmu0hatdl + dThetaHatdl) &
      * daI_mu0hatdl + (1d0-SI)*aI/mu0hat*(                 &
      - (d2mu0hatdl2/mu0hat - dmu0hatdl**2/mu0hat**2)*ThetaHat  &
      + (alphaI-dmu0hatdl/mu0hat)*dThetaHatdl + d2ThetaHatdl2)  )
   
       ! Light dependent growth rate 
       ! (needs to take into account the cost of dark and light-dependent chl maintenance)
       mu0hatSI = mu0hat*SI  ! Gross specific carbon uptake (photosynthesis)
       muIhat   = mu0hatSI-(mu0hatSI+RMchl)*zetaChl*ThetaHat ! Net specific carbon uptake
       muIhat   = max(muIhat,1D-10*mu0hatSI)
   
       dmu0hatSIdl   = SI*dmu0hatdl + mu0hat*dSIdl
       d2mu0hatSIdl2 = d2mu0hatdl2*SI+2d0*dmu0hatdl*dSIdl+mu0hat*d2SIdl2 !Correct
       dmuIhatdl = (1D0-zetaChl*ThetaHat)*dmu0hatSIdl - dThetaHatdl*zetaChl*(mu0hatSI+RMchl) !Correct
    
       d2muIhatdl2=d2mu0hatSIdl2                                      &
       -zetaChl*(ThetaHat*d2mu0hatSIdl2+2d0*dThetaHatdl*dmu0hatSIdl   &
       +mu0hatSI*d2ThetaHatdl2)-zetaChl*RMchl*d2ThetaHatdl2   !Correct
    
       ZINT   = Qs*(muIhat/VNhat + zetaN)
    
       dZINdl = Qs*(dmuIhatdl/VNhat - muIhat*dVNhatdl/VNhat**2)+alphaQ*ZINT    
       
       d2ZINdl2 = Qs/VNhat*((alphaQ-dVNhatdl/VNhat)*dmuIhatdl & 
       + d2muIhatdl2) - Qs/VNhat**2*(muIhat*d2VNhatdl2        &
       + dVNhatdl*(dmuIhatdl+alphaQ*muIhat-2d0*muIhat          &
       / VNhat*dVNhatdl)) + alphaQ*dZINdl
       
       fV = (-1d0 + sqrt(1d0 + 1d0/ZINT))*Qs*muIhat/VNhat
       fV = max(fV, 0.01)
    !
       else
    ! Under the conditions of no light:
          ThetaHat      = 0.01  !  a small positive value 
          dThetaHatdl   = 0d0
          d2ThetaHatdl2 = 0d0
          ZINT          = Qs*zetaN
          dZINdl        = alphaQ*ZINT
          d2ZINdl2      = alphaQ*dZINdl
          fV            = 0.01
          muIhat        = 0d0
          dmuIhatdl     = 0d0
          d2muIhatdl2   = 0d0
       endif
    
          ! Optimal nutrient quota:
       QN = (1d0+ sqrt(1d0+1d0/ZINT))*Qs
    
       dQNdl  = alphaQ*QN-dZINdl*Qs/(2d0*ZINT*sqrt(ZINT*(1d0+ZINT))) !confirmed  
    
       d2QNdl2 = alphaQ*dQNdl - Qs/(2d0*ZINT*sqrt(ZINT*(ZINT+1d0)))  &
        *(d2ZINdl2+alphaQ*dZINdl-(2d0*ZINT+1.5d0)/(ZINT*(ZINT+1d0))*dZINdl**2)      ! Confirmed
    
       dfVdl = alphaQ*Qs*(1d0/QN+2d0*zetaN)-(zetaN+Qs/QN**2)*dQNdl  !Confirmed
    
       d2fVdl2 = (alphaQ**2)*Qs*(1d0/QN + 2d0*zetaN)                    &
         -  2d0*alphaQ*Qs*dQNdl/(QN*QN) + 2d0*(dQNdl**2)*Qs/(QN*QN*QN)    &
         -     (zetaN + Qs/(QN*QN)) * d2QNdl2  ! Confirmed
    
       if (par_ .gt. I_zero) then
        !Net growth rate (d-1) of phytoplankton at the average size
          !muNet = muIhat*X - zetaN*fV*VNhat
          muNet = muIhat*(1D0 - 2D0 * Qs/QN)

!Here the derivative of muNet includes respiratory costs of both N Assim and Chl maintenance       
          dmuNetdl = muIhat*(Qs/(QN**2)*dQNdl                           &
      -  alphaQ*Qs/QN-dfVdl) +  X*dmuIhatdl                             & 
      -  zetaN*(fV*dVNhatdl+VNhat*dfVdl)

         d2muNetdl2 = (Qs/(QN*QN))*dQNdl*dmuIhatdl                      &
     + muIhat*(Qs/QN**2)*(dQNdl*(alphaQ - 2d0*dQNdl/QN) + d2QNdl2)      & 
     - (alphaQ*(Qs/QN*(dmuIhatdl + alphaQ*muIhat)                       &
     - muIhat*Qs/(QN**2)*dQNdl))                                        &
     - (muIhat*d2fVdl2+dfVdl*dmuIhatdl)                                 &
     + dmuIhatdl*(Qs/(QN**2)*dQNdl-alphaQ*Qs/QN-dfVdl)  & 
     + X*d2muIhatdl2                   &
     - zetaN*(fV*d2VNhatdl2 + 2d0*dfVdl*dVNhatdl + VNhat*d2fVdl2)  !dC

       else

        muNet      = 0d0
        dmuNetdl   = 0d0
        d2muNetdl2 = 0d0
       endif
          X  = 1d0-fV-Qs/QN
!  chl:C ratio [ g chl / mol C ] of the whole cell, at the mean cell size 
       Theta = ThetaHat*X
! Calculate the mean chl:C ratio of phytoplankton (averaged over all sizes) 
    ! using the second derivative of the cell quota w.r.t. log size. (Why negative?)
       ThetaAvg = max(Theta + (VAR/2d0)*(d2ThetaHatdl2*X                &     
       - ThetaHat*(d2fVdl2 + Qs*(dQNdl**2*(2d0/QN)-d2QNdl2)/(QN*QN))),  &
       params(ithetamin))

      if (ThetaAvg .le. 0.01 .and. k .lt. nlev) then
        ThetaAvg = Varout(oTheta(1),k+1)
     endif
            
    ! Calculate the mean N:C ratio of phytoplankton (averaged over all sizes) 
    ! using the second derivative of the cell quota w.r.t. log size. (How to derive?)
     QNAvg = max(QN/(1d0+((2d0/QN)*dQNdl**2 - d2QNdl2)*VAR/(2d0*QN)),   &
        params(iQNmin))  
!=============================================================
  !  ! Calculate the community sinking rate of phytoplankton
  !  ! Phytoplankton sinking rate at the average size
  !  w_p    = ScaleTrait(PMU,abs(w_p0),alphaW)  !Positive

  !  ! Constrain the sinking rate not too large for large cells (Smayda 1970)
  !  w_p    = min(w_p, Wmax)
  !  dwdl   = alphaW*w_p
  !  d2wdl2 = alphaW*dwdl
 
  !! Sinking rate (m) of phytoplankton
  !  w_pAvg = w_p+0.5*VAR*d2wdl2
  !! sinking rate must be negative!
  !  w_pAvg = -min(w_pAvg, Wmax)
    
!=============================================================
!! ZOOplankton section:
   tf_z = TEMPBOL(Ez,tC)
  
!   ! The grazing dependence on total prey
   gbar = grazing(grazing_formulation,params(ikp),PHY)

   !Zooplankton per capita total ingestion rate
    INGES = tf_z*dtdays*params(igmax)*gbar

   !Zooplankton excretion rate (-> DOM)
   RES = INGES*(1d0-GGE-unass)

   !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass
    
! Grazing rate on the mean size of PHY (specific to N-based Phy biomass, unit: d-1) (Eq. 12)
   gbar      = INGES*ZOO/PHY*sqrt(params(ialphaG))
   dgdlbar   = 0d0
   d2gdl2bar = (1d0-params(ialphaG))/VAR*gbar
!!End of zooplankton section
!=============================================================
!! Solve ODE functions:
   Zmort = ZOO*ZOO*dtdays*params(imz)*tf_z  !Mortality term for ZOO

   if (d2gdl2bar .lt. 0d0) then
      VAR1 = VAR*(1d0-VAR*d2gdl2bar)
   else                                       
      VAR1 = VAR/(1d0+VAR*d2gdl2bar)
   endif

   ! Update PMU and VAR:    
   if (d2muNetdl2 .gt. 0d0) then
! self%dVARdt=VAR*VAR*(self%d2muNetdl2- self%d2gdl2bar-self%d2wdl2)  !Eq. 19 
      VAR1 =VAR1*(1d0+VAR*d2muNetdl2)
   else    
      VAR1 =VAR1/(1d0-VAR*d2muNetdl2)
   endif
!   VAR1    =VAR1/(1d0+VAR*d2wdl2)
   
   !Eq. 18, Update PMU:
   !   self%dPMUdt = self%p%VAR*(self%dmuNetdl - self%dgdlbar - self%dwdl)
   PMU = PMU + log(1D1)  ! Restore to positive values
   if (dmuNetdl .gt. 0d0) then
      PMU1 = PMU  + VAR * dmuNetdl
   else
      PMU1 = PMU/(1d0-VAR/PMU*dmuNetdl)
   endif
   
!   PMU1 = PMU1/(1d0 + VAR/PMU*dwdl)
   !Contrain the PMU and VAR: 
   PMU1 = min(max(PMU1,0.01),PMUmax)
   VAR1 = min(max(VAR1,0.01),VARmax)
          
    ! For production/destruction matrix:
      
   pp_ND=dtdays*params(irDN)*DET*tf_z   
   pp_NZ=ZOO*RES        
   pp_DZ=ZOO*EGES+Zmort 
   pp_ZP=ZOO*INGES      
   PPpn =PHY*(muNet + 0.5*VAR*d2muNetdl2)
   pp_PN=max(0.5*( PPpn + abs(PPpn)), 0d0)
   pp_NP=max(0.5*(-PPpn + abs(PPpn)), 0d0)
   DET1 = (DET + pp_DZ)/(1d0 + dtdays*params(irDN)*tf_z)
   DET1 = max(DET1,1D-16)
   NO31 = (NO3+pp_ND+pp_NZ+pp_NP)/(1d0+pp_PN/NO3)
  
   !if (NO31 .le. 1D-10) then
   !   write(6,*) 'NO3=',NO3
   !   write(6,*) 'muNet=',muNet
   !   write(6,*) 'd2muNetdl2=',d2muNetdl2
   !   write(6,*) 'VAR=',VAR
   !   write(6,*) 'pp_PN=',pp_PN
   !   stop
   !endif
   PHY1 = (PHY+pp_PN)/(1d0+(pp_ZP+pp_NP )/PHY)
   ZOO1 = (ZOO+pp_ZP)/(1d0+ EGES+ZOO*dtdays*params(imz)*tf_z+RES)
   
   PMUPHY1 =(PMUPHY+PHY*PMU1+PMU*PHY1)/(1d0+ 2d0*PHY*PMU/PMUPHY)
   VARPHY1 =(VARPHY+PHY*VAR1+VAR*PHY1)/(1d0+ 2d0*PHY*VAR/VARPHY)
   
   Varout(oNO3,k)   = NO31
   Varout(oPHY(1),k)= PHY1
   Varout(oZOO,k)   = ZOO1

   ! Calculate total CHL:
   Varout(oCHLt,k)  = PHY1 / QNAvg * ThetaAvg
   Varout(oPPt,k)   = pp_PN/ QNAvg / dtdays
   Varout(oDET,k)   = DET1
   Varout(oPMU,k)   = PMUPHY1
   Varout(oVAR,k)   = VARPHY1

   if(do_IRON) Varout(oFER,k) = Fe

   Varout(oTheta(1),k) = ThetaAvg
   Varout(oQN(1)   ,k) = QNAvg
   Varout(omuNet(1),k) = PP_PN/dtdays/PHY
   Varout(oGraz(1) ,k) = PP_ZP/dtdays/PHY
   Varout(oZ2N     ,k) = PP_NZ/dtdays
   Varout(oD2N     ,k) = PP_ND/dtdays
   Varout(odmudl   ,k) = dmuNetdl/dtdays
   Varout(od2mu    ,k) = d2muNetdl2/dtdays
   Varout(od2gdl   ,k) = d2gdl2bar/dtdays
   !Varout(ow_p(1)  ,k) = w_pAvg  !Should not corrected by dtdays

   ! Calculate size-fractionated Chl
   ! Calculate the proportions of different phyto size fractions
   ! Here the PMU1 is already the normal PMU + log(10) (to avoid negative values)
   pCHL(4) = pnorm(PMU_1 , PMU1, sqrt(VAR1))
   pCHL(3) = pnorm(PMU_3 , PMU1, sqrt(VAR1))
   pCHL(2) = pnorm(PMU_10, PMU1, sqrt(VAR1))
   pCHL(1) = 1d0    -pCHL(2)
   pCHL(2) = pCHL(2)-pCHL(3)
   pCHL(3) = pCHL(3)-pCHL(4)
   do j = 1, 4
      !Varout(oCHLs(j),k) = pCHL(j) * Varout(oCHLt,k) ! Absolute concentration
      Varout(oCHLs(j),k) = pCHL(j) ! Percentage of one size class
   enddo

ENDDO
END SUBROUTINE FLEXEFT_CONT 
