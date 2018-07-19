subroutine EFT_size(PMU, NO3, tC, par_, muNet, QN, Theta, LSI, Lno3, ThetaHat)
use bio_MOD
implicit none
!INPUT PARAMETERS:
real, intent(in)    :: PMU, NO3, tC, par_

!Output parameters:
real, intent(out)   :: muNet, QN, Theta, LSI, Lno3, ThetaHat
!LOCAL VARIABLES of phytoplankton:
real :: X
real :: larg   !Environmental variables
real :: zetaN, zetaChl
real :: KFe,Fe
real :: V0hat,Kn,A0N,A0hat,fA,muIhat
real :: VNhat
real :: mu0hat,mu0hatSI,RMchl
real :: fV, SI
real :: ZINT
real :: aI
real :: Qs
real :: larg1, W_Y
logical, parameter :: do_IRON = .false.
real,    external  :: WAPR
!-----------------------------------------------------------------------
!Warning: All rates should be multiplied by dtdays to get the real rate

alphamu = params(ialphamu)
betamu  = params(ibetamu)
tf_p    = TEMPBOL(params(iEp),tC)
zetaChl = 0.8 !params(izetaChl)
zetaN   = 0.6 !params(izetaN)

! Fe related:

if (do_IRON) then
   Fe   = max(Fe,Femin)  !Dissolved Fe concentration
   KFe  = ScaleTrait(PMU, K0Fe,alphaFe) !The half saturation constant of Fe at average size
endif
    
Qs      = ScaleTrait(PMU, params(iQ0N), params(ialphaQ))/2d0
mu0hat  = dtdays*tf_p*mu0 * exp(alphamu*PMU + betamu*PMU*PMU)

! X: a scratch variable
X       = 0d0

if (do_IRON) then
  mu0hat= mu0hat * Fe/(Fe + KFe)
  X     = alphaFe*KFe/(Fe + KFe) 
endif
  
  ! Iron limits nitrogen uptake
  V0hat=ScaleTrait(PMU, dtdays*tf_p*V0N, alphaV)
  
  if (do_IRON) V0hat = V0hat*Fe/(Fe + KFe)
  
  ! Initial slope of P-I curve
  aI = ScaleTrait(PMU, dtdays*params(iaI0), params(ialphaI))
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
      Kn     = ScaleTrait(PMU,params(iKN), params(ialphaK)) 
      VNhat  = V0hat*NO3/(NO3 + Kn)
  
     ! case 2: optimal uptake based on Pahlow (2005) and Smith et al. (2009)
    case(2)
  
    !  Lmin  = log((params(iLref)**3)/6d0*pi) &
    !  + log(1D0 +params(iPenfac))/( params(iPenfac)*params(ialphaA)) 

      A0N   = dtdays*tf_p*params(iA0N)
      A0hat = ScaleTrait(PMU, A0N, params(ialphaA))

    !  A0hat = PenAff(PMU, params(ialphaA), params(iPenfac),Lmin) &
    !   * ScaleTrait(PMU, A0N, params(ialphaA))

      A0hat = max(A0hat,1D-9*A0N)
 
      !Define fA
      fA = 1D0/( 1D0 + sqrt(A0hat*NO3/V0hat) ) 
      VNhat = (1D0-fA)*V0hat*fA*A0hat*NO3/((1D0-fA)*V0hat + fA*A0hat*NO3) 
      
    case default
     write(6,*) 'Error: Incorrect option for nutrient uptake!'
     stop
  ENDSELECT  

! Calculate thetahat (optimal g Chl/mol C for the chloroplast under nutrient replete conditions)
! Only calculate within the euphotic zone, otherwise many numerical problems.

  if(par_ .gt. (I_zero+.1)) then    !To avoid large C:Chl ratios
    
    larg1   = exp(1d0 + min(aI*par_/(mu0hat*zetaChl),6d2))
    larg    = (1d0 + RMchl/mu0hat)*larg1   
   W_Y      = WAPR(larg)
   ThetaHat = 1d0/params(izetaChl) + (1d0- W_Y)*mu0hat/(aI * par_)

    ! Effect of light limitation
    SI = 1d0 - max(exp(-aI*par_*ThetaHat/mu0hat),0d0)

! Light dependent growth rate 
! (needs to take into account the cost of dark and light-dependent chl maintenance)
    mu0hatSI = mu0hat*SI  ! Gross specific carbon uptake (photosynthesis)
    muIhat   = mu0hatSI-(mu0hatSI+RMchl)*zetaChl*ThetaHat ! Net specific carbon uptake
    muIhat   = max(muIhat,1D-9*mu0hat)
    LSI      = 1D0-SI
    ZINT     = Qs*(muIhat/VNhat + zetaN)
    fV       = (-1d0 + sqrt(1d0 + 1d0/ZINT))*Qs*muIhat/VNhat
    fV       = max(fV,0.01)
    
  else
    ! Under the conditions of no light:

       ThetaHat      = 1d-2  !  a small positive value 
       ZINT          = Qs*zetaN
       fV            = 1d-2
       muIhat        = 0d0
       LSI           = 1D0
    endif
     ! Nutrient limitation index:
     Lno3 =1d0/(1d0 + sqrt(1D0 +1D0/ZINT)) 
     ! Optimal nutrient quota:
     QN = Qs/Lno3
    
     if (par_ .gt. I_zero) then
!Net growth rate (d-1) of phytoplankton at the average size
      ! muNet = muIhat*(1d0-fV-Qs/QN) - zetaN*fV*VNhat
       muNet = muIhat*(1d0-2d0*Lno3)
     else
       muNet = 0d0
     endif
! Chl:C ratio [ g chl / mol C ] of the whole cell
   Theta  = ThetaHat*(1D0-fV-Qs/QN)
   
! Phytoplankton sinking rate at the average size 
! (should not normalize to dtdays)
   !w_p    = ScaleTrait(PMU,abs(w_p0),alphaW)  !Positive
! Constrain the sinking rate not too large for large cells (Smayda 1970)
   !wPHY   = min(w_p, Wmax)
end subroutine EFT_size

