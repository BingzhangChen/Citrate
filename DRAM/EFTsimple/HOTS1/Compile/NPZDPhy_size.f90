subroutine NPZDPhy_size(PMU, NO3, tC, par_, muNet, SI, Lno3, QN, theta)
use bio_MOD
implicit none

!INPUT PARAMETERS:
real, intent(in)    :: PMU, NO3, tC, par_

!Output parameters:
real, intent(out)   :: muNet, SI, Lno3

!Local variable:
real                :: Kn, mu0hat, aI, cff, cff1
real                :: Qmax, fN, QN, theta
real, parameter     :: alphaK = 0.27
real, parameter     :: Qmin   = 0.06
real, parameter     :: thetamin = 0.02

tf_p    = TEMPBOL(Ep,tC)
mu0hat  = dtdays*tf_p * params(imu0) * exp(alphamu*PMU + betamu * PMU**2)

Kn      = params(iKN) * exp( alphaK * PMU)

! Effect of light limitation
aI      = params(iaI0) * exp(params(ialphaI)*PMU)
SI      = 1d0 - max(exp(-params(iaI0) * par_ ),0d0)
Lno3    = NO3/(NO3+Kn)
muNet   = mu0hat * Lno3 * SI

Qmax    = 3.*Qmin
cff1    = 1d0-Qmin/Qmax
cff     = 1d0-cff1*fN
!N:C ratio
QN      = Qmin/cff

!Chl:C ratio 
cff     = (thetamax - thetamin)/PAR_
theta   = thetamin+muNet/aI*cff   !Unit: gChl/molC

return
end subroutine NPZDPhy_size

