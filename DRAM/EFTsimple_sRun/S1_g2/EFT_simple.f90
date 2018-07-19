SUBROUTINE FlexEFT_simple
use bio_MOD
implicit none
integer :: k
real    :: Qs,QN, mu0hat,V0hat, aI, RMchl
real    :: VNhat
real    :: A0hat,fA,larg1,larg
real    :: W_Y,ThetaHat,SI,mu0hatSI,muIhat,par_
real    :: Lno3, Ep, Ez, zetaN, zetaChl
real    :: NO3, PHY, ZOO, DET
real    :: pp_PN
real    :: ZINT,fV,muNet
real, external :: WAPR

Ep     = params(iEp)
Ez     = params(iEz)
zetaN  = params(izetaN)
zetaChl= params(izetaChl)

DO k=nlev,1,-1    ! k cannot be reversed
   DET  = Vars(iDET,k)
   NO3  = Vars(iNO3,k)
   PHY  = Vars(iPHY(1),k)
   ZOO  = Vars(iZOO,k)  !Zooplankton biomass

   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif

! Phytoplankton section:
   tf_p = TEMPBOL(Ep,Temp(k))
   Qs   = params(iQ0N)/2d0
 mu0hat = dtdays*tf_p*mu0
  V0hat = dtdays*tf_p*V0N
  ! Initial slope of P-I curve
  aI    = dtdays*params(iaI0)
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
      VNhat  = V0hat*NO3/(NO3 + params(iKN))

 ! case 2: optimal uptake based on Pahlow (2005) and Smith et al. (2009)
    case(2)
      A0hat = dtdays*tf_p*params(iA0N)
      !Define fA
      fA    = 1D0/(1D0 + sqrt(A0hat*NO3/V0hat) ) 
      VNhat = (1D0-fA)*V0hat*fA*A0hat*NO3/((1D0-fA)*V0hat+fA*A0hat*NO3) 
    case default
     write(6,*) 'Error: Incorrect option for nutrient uptake!'
     stop
  ENDSELECT  

! Calculate thetahat (optimal g Chl/mol C for the chloroplast under nutrient replete conditions)
! Only calculate within the euphotic zone, otherwise many numerical problems.
  if( par_ .gt. (I_zero+.1) ) then
    larg1 = exp(1d0 + aI*par_/(mu0hat*zetaChl))
    larg  = (1d0 + RMchl/mu0hat)*larg1   
     W_Y  = WAPR(larg)
   ThetaHat = 1d0/zetaChl + (1d0- W_Y)*mu0hat/(aI * par_)

    ! Effect of light limitation
    SI = 1d0 - max(exp(-aI * par_ * ThetaHat/mu0hat),0d0)
    Varout(oSI(1),k) = SI
! Light dependent growth rate 
! (needs to take into account the cost of dark and light-dependent chl maintenance)
    mu0hatSI = mu0hat*SI !Gross specific carbon uptake (photosynthesis)
    muIhat   = mu0hatSI-(mu0hatSI+RMchl)*zetaChl*ThetaHat ! Net specific carbon uptake
    muIhat   = max(muIhat,0d0)
   
    ZINT   = Qs*(muIhat/VNhat + zetaN)
           
    fV = (-1d0 + sqrt(1d0 + 1d0/ZINT))*Qs*muIhat/VNhat
    fV = max(fV,0.01)
   
    else
    ! Under the conditions of no light:
    !   ThetaHat = 1d-2  !  a small positive value 
       ZINT     = Qs*zetaN
       muIhat   = 0d0
       Varout(oSI(1),k) = 0d0
    endif
      ! Nutrient limitation index:
     Lno3 =1d0/(1d0 + sqrt(1D0 +1D0/ZINT)) 

     Varout(oLno3(1),k) = Lno3
     ! Optimal nutrient quota:
     QN   = Qs/Lno3
     Varout(oQN(1),k)   = QN
    
     if (par_ .gt. (I_zero+.1)) then
   !Net growth rate (d-1) of phytoplankton at the average size
       muNet = muIhat*(1d0-2d0*Lno3)  !Droop model
   !Chl:C ratio [ g chl / mol C ] of the whole cell
       Varout(oTheta(1),k) = ThetaHat*(1D0-fV-Qs/QN)

     else
       muNet = 0d0
      !Chl:C ratio [ g chl / mol C ] of the whole cell
       Varout(oTheta(1),k) =  Varout(oTheta(1),k+1)
     endif

   !Save net growth rate
    Varout(omuNet(1),k) = muNet/dtdays
    
   !Phytoplankton sinking rate at the average size
   ! Varout(ow_p(1),  k) = abs(w_p0)  !Positive
!---------------------------------------------------------------
!! ZOOplankton section:
   tf_z = TEMPBOL(Ez,Temp(k))

 ! The grazing dependence on total prey (dimensionless)
   gbar = grazing(grazing_formulation,params(ikp),PHY)
 !Zooplankton total ingestion rate
   INGES = tf_z*dtdays*params(igmax)*gbar

 !Zooplankton excretion rate (-> DOM)
   RES  = INGES*(1d0-GGE-unass)

 !ZOOPLANKTON EGESTION (-> POM)
   EGES = INGES*unass
    
! Grazing rate on PHY(specific to N-based Phy biomass, unit: d-1)

   ! Calculate the specific grazing rate for PHY
   Varout(oGraz(1),k) = INGES*ZOO/PHY/dtdays
   
   Varout(oPHY(1),k)  = PHY*(1d0+muNet)/(1d0 + INGES*ZOO/PHY)
!!End of zooplankton section
!=============================================================
!! Solve ODE functions:
  Zmort = ZOO*ZOO*dtdays*params(imz)*tf_z  !Mortality term for ZOO
 
  ! For production/destruction matrix:
  pp_PN=PHY*muNet
  pp_ND=dtdays* params(irDN)*DET*tf_z   
  pp_NZ=ZOO*RES        
  pp_DZ=ZOO*EGES+Zmort 
  pp_ZP=ZOO*INGES      
  
  Varout(oDET,k) = (DET+pp_DZ)/(1d0 + dtdays * params(irDN)*tf_z)
  Varout(oNO3,k) = (NO3+pp_ND+pp_NZ)/(1d0+pp_PN/NO3)

  Varout(oZOO,k) = (ZOO+pp_ZP)/(1d0   &
                 + EGES+ZOO*dtdays*params(imz)*tf_z+RES)
    
  Varout(oZ2N,k) = pp_NZ/dtdays
  Varout(oD2N,k) = pp_ND/dtdays
  Varout(oPPt,k) = pp_PN/dtdays/QN
  Varout(oCHLt,k)= Vars(iPHY(1),k)/Varout(oQN(1),k)*Varout(oTheta(1),k)
ENDDO
End subroutine FlexEFT_simple

