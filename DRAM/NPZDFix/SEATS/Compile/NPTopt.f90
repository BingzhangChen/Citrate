Module NPTopt_MOD
implicit none

CONTAINS

SUBROUTINE BIOLOGY
USE PARAM_MOD
IMPLICIT NONE
integer :: k, i

!INPUT PARAMETERS:
real :: par_

!Local variables
real :: KN, mu0, aI0, Mm, Emort, p, NO3, NPPt = 0.
real, parameter :: Ez = 0.65, phi = 0.5, PHYmax = 200. 

KN   = exp(params(iKN))
mu0  = exp(params(imu0))
aI0  = exp(params(iaI0))
Mm   = exp(params(iMm))

DO k = nlev, 1, -1   
   ! Check whether in the MLD or not
   if (k < N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif

   Varout(oPAR_,k) = par_
   NO3    = Vars(iNO3,    k)

!=============================================================
!! Solve ODE functions:
!All rates have been multiplied by dtdays to get the real rate correponding to the actual time step
   Tchl      = 0.
   mu_avg    = 0.
   mumax_avg = 0.
   PHYt      = sum(Vars(iPHY(:,k)))
   NPPt      = 0.
   do j = 1, NG
      do i = 1,NTAXA

         p = Vars(iPHY((j-1)*NG+i), k)

         IF (mod(it, nsave) .EQ. 1) THEN
           !Use nonmixed light to calculate NPP

           call growth(temp(k), NO3, PAR(k), mu0, Topt(i), KN,         &
                     aI0, mu, QN, theta, mumax)
           NPPt = NPPt + p * mu/QN * 12.
         ENDIF

         !Phytoplankton growth rate
         call growth(temp(k), NO3, par_, mu0, Topt(i), KN,         &
                     aI0, mu, QN, theta, mumax)

         Tchl     = Tchl      + p / QN * theta
         mu_avg   = mu_avg    + mu*p
         mumax_avg= mumax_avg + mumax*p

         !Maintain diversity using Record et al. (2014)
         Emort    = exp(TK(temp(k))*Ez)
         mort     = Emort * Mm * p**(1.+phi) * PHYt**(1.-phi)

         p = p + dtdays*(mu*p - mort)
         Vars(iPHY((j-1)*NG+i), k) = min(max(p, 0.), PHYmax) 

         !Compute total phyto. uptake - mortality
         total_uptake = total_uptake + (mu*p - mort)
       enddo
   enddo
   mu_avg    = mu_avg/PHYt
   mumax_avg = mumax_avg/PHYt

   ! Save some model outputs:
   Varout(oCHLt, k) = Tchl 
   Varout(oNPP,  k) = NPPt

   NO3 = NO3 + dtdays*(-total_uptake) 
   NO3 = min(max(NO3, 0.), PHYmax)
   Varout(oNO3,  k) = NO3

ENDDO
END SUBROUTINE BIOLOGY

!phytoplankton growth rate, Chl:C and N:C function
subroutine growth(t, N, par,r0,To, KN, aI0, gr, QN, theta, um)
IMPLICIT NONE
real, intent(in)  :: t     !Environmental Temperature (ºC)
real, intent(in)  :: N     !Nitrate (mmol m-3)
real, intent(in)  :: par   !light (mol E m-2 d-1)
real, intent(in)  :: r0    !Normalized growth rate for mumax (d-1)
real, intent(in)  :: To    !Optimal temperature (ºC)
real, intent(in)  :: KN    !Half-saturation constant  (mmol m-3)
real, intent(in)  :: aI0   !Slope of P_I curve (W m-2)-1 d-1
real, intent(out) :: gr    !Specific growth rate (d-1)
real, intent(out) :: QN    !Phyto. N:C ratio
real, intent(out) :: theta !Phyto. Chl:C ratio (gChl molC-1)
real, intent(out) :: um    !Phyto. maximal growth rate (excluding N and light limitation)

real, parameter   :: Q0N      = 0.05  !Minimal N:C ratio
real, parameter   :: Qm       = 0.15  !Maximal N:C ratio
real, parameter   :: thetamin = 0.02  !Minimal Chl:C ratio (gChl molC-1)
real, parameter   :: thetamax = 0.63  !Maximal Chl:C ratio (gChl molC-1)
real, parameter   :: Ea       = 0.81  !Activation energy (eV)
real, parameter   :: dED      = 1.74  !Eh - Ea (eV)

real :: par_, SI, fN


um   = Johnson(t, mu0(To, r0), Ea, dED, To)
fN   = N / (N + KN)
if (um < 0.01) then
  SI = 1.
else
  SI = 1. - exp(-aI0*par_/um)
endif
QN   = Q0N/(1.-(1.-Q0N/Qm)*fN)
gr   = um*SI*fN  !Growth rate
theta= thetamin + gr/par_/aI0 * (thetamax-thetamin)   !Unit: gChl/molC

return
end subroutine growth

real function Johnson(tC, mu0, Ea, dED, To) result(y)
implicit none
real,   intent(in)  :: tC, mu0, Ea, dED, To
real,   parameter   :: kb   = 8.62D-5
real,   parameter   :: T0   = 273.15D0
real,   parameter   :: Tref = 15D0
real                :: ED, a, b

if (dED .le. 0d0) then
  stop "dED must be positive!"
else
  ED = dED + Ea
   a = Ea/dED
   b = exp(ED/kb*(1D0/(To+T0)-1D0/(tC+T0)))
   b = a*b
   a = Ea/kb*(1D0/(Tref+T0) - 1D0/(tC+T0))
   y = mu0*exp(a)/(1D0+b)   
endif
return
End function Johnson

pure real function TK(tC)
implicit none
!DESCRIPTION:
!The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
! tC: in situ temperature
! Tr: reference temperature
!
!INPUT PARAMETERS:
real, intent (in) :: tC
! boltzman constant constant [ eV /K ]
real, parameter   :: kb = 8.62d-5, Tr = 15.0

TK = -(1./kb)*(1./(273.15 + tC) - 1./(273.15 + Tr))
return 
end function TK

real function mu0(x, mu0p)
implicit none
real, intent(in) :: x           !Topt in ºC
real, intent(in) :: mu0p        !Normalized growth rate
real, parameter  :: E0 = -0.41  !Nominal AE of growth rate normalizations
mu0 =  mu0p * exp(TK(x) * E0) 
end function mu0

END Module NPTopt_MOD
