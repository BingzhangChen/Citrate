subroutine GeiderPhy_size(PMU, NO3, tC, par_,theta, muNet, SI, Lno3)
use bio_MOD
implicit none

!INPUT PARAMETERS:
real, intent(in)    :: PMU, NO3, tC, par_, theta

!Output parameters:
real, intent(out)   :: muNet, SI, Lno3

!Local variable:
real                :: Kn, mu0hat

alphamu = params(ialphamu)
tf_p    = TEMPBOL(Ep,tC)

mu0hat  = dtdays*tf_p * params(imu0) * exp(alphamu*PMU)

Kn      = params(iKN) * exp( params(ialphaKN) * PMU)

!The light limitation index
SI = 1d0 - exp(-params(iaI0) * dtdays *par_ * theta/mu0hat)

!The nutrient limitation index
Lno3    = NO3/(NO3+Kn)
muNet   = mu0hat * Lno3 * SI
return
end subroutine GeiderPhy_size

