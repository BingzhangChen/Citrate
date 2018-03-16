subroutine MONOD(Temp,PAR,NO3,PO4,mu0,Qnmin,Qpmin,aI0,bI0,KN,KP,DFe,KFe,muNet,QN,QP,theta,SI,Lno3)
! This subroutine calculate nondiatroph phytoplankton growth rates, QN and Theta based on the simplest MONOD equation
use bio_MOD, only: TEMPBOL, Ep, do_IRON, N2fix,thetamax
real, intent(in)  :: Temp, PAR, NO3,PO4,mu0,Qnmin,Qpmin,aI0,bI0,KN,KP,DFe,KFe
real, intent(out) :: QN, QP, theta, SI, Lno3,muNet
real              :: rmax_T,Qnmax,Qpmax,Lpo4
real, parameter   :: thetamin = 0.02

  ! Maximal N:C ratio
   Qnmax=3.*Qnmin
   Qpmax=3.*Qpmin

  ! The maximal growth rate (rmax_T) under temperature tC 
   rmax_T = mu0*TEMPBOL(Ep,Temp)

  !The light limitation index (SI)
  ! Unit of aI0_C: (W m-2)-1, different from other models (aI0_Chl)
  ! Include photoinhibition (Platt et al. J Phycol 1980)
  SI = (1d0 - exp(-aI0*PAR/rmax_T))*exp(-bI0*PAR/rmax_T)

  ! Growth rate (d-1) is the product function of temperature, light, and nutrient   
  Lno3  = NO3/(NO3 + KN)
  if (N2fix) then
     Lpo4 = PO4/(PO4 + KP)
     Lno3 = min(Lno3, Lpo4)
  endif

  if (DO_IRON) then
     Lpo4  = DFe/(DFe + KFe)
     Lno3  = min(Lno3, Lpo4)
  endif
  muNet = rmax_T*Lno3*SI
  QN    = Qnmin/(1d0-(1d0-Qnmin/Qnmax)*Lno3)
  if (N2fix) QP=Qpmin/(1d0-(1d0-Qpmin/Qpmax)*(PO4/(PO4+KP)))
  
  theta  = thetamin+muNet/PAR/aI0*(thetamax-thetamin)   !Unit: gChl/molC
  return
end subroutine

