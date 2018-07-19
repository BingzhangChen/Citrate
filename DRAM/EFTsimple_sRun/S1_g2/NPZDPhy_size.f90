subroutine NPZDPhy_size(PMU, NO3, tC, par_, muNet)
use bio_MOD
implicit none

!INPUT PARAMETERS:
real, intent(in)    :: PMU, NO3, tC, par_

!Output parameters:
real, intent(out)   :: muNet

!Local variable:
real                :: Kn, SI, mu0hat

alphamu = params(ialphamu)
tf_p    = TEMPBOL(params(iEp),tC)
mu0hat  = dtdays*tf_p * params(imu0) * exp(alphamu*PMU)

Kn      = params(iKN) * exp( params(ialphaKN) * PMU)

! Effect of light limitation
SI      = 1d0 - max(exp(-params(iaI0) * par_ ),0d0)

muNet   = mu0hat * NO3 / (NO3 + Kn) * SI
return
end subroutine NPZDPhy_size

