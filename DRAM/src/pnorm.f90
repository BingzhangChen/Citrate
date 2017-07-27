real function pnorm(x,mean,sd)
! This function calculates the area of the tail left to x 
! of the Gaussian curve with mean and sd
implicit none
real, intent(in) :: x, mean, sd
real             :: x_
real, parameter  :: inv_sqrt2 = 1d0/sqrt(2d0)

x_    = (x-mean)/sd
pnorm = 0.5*(1d0+erf(x_*inv_sqrt2))
return
end function pnorm

