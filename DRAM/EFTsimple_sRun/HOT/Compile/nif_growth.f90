subroutine nif_growth(PO4,PAR,temp,tau,mu,theta)
use bio_MOD, only : params, iKPnif, TEMPBOL,Ep
implicit none
real, intent(in)  :: PO4, PAR, temp,tau
real, intent(out) :: mu, theta
real, parameter   :: tau_crit = 0.062 !Critical wind stress (N m-2)
real, parameter   :: tem_crit = 24.75 !Critical temperature (ºC)
real, parameter   :: aI0      = 0.01  !Initial slope of the P-I curve
real, parameter   :: mumax    = 0.3   !Maximal growth rate of diazotrophs at 15 ºC
real, parameter   :: thetamin = 0.02, thetamax = 0.47

! Amplification factor
real :: ampl
real :: tf
real :: SI  ! Light limitation index
real :: Kpo4nih  = 0.05  !Half-saturation constant for PO4 uptake by nif

tf      = TEMPBOL(Ep,temp)
Kpo4nih = params(iKPnif)
! This subroutine follows Fennel et al. (2002) to incorporate PO4
! Diatroph growth rates depends on temperature, light, and wind stress
ampl = tanh(2.*(temp - tem_crit)) + 2.
!if (tau .le. tau_crit) then
   ampl = ampl/3.
!else
!   ampl = ampl/6.
!endif

!The light limitation index (SI)
! Unit of aI0: (W m-2)-1 d-1
SI    = 1d0 - exp(-aI0*PAR/mumax)
mu    = ampl*SI*PO4/(PO4+Kpo4nih)*mumax*tf

!if (mu < 0.) then
!   write(6,*) 'ampl = ', ampl
!   write(6,*) 'SI   = ', SI
!   write(6,*) 'LnoP = ', PO4/(PO4+Kpo4nih)
!   stop
!endif
theta = thetamin+mu/PAR/aI0*(thetamax-thetamin)   !Unit: gChl/molC
end subroutine
