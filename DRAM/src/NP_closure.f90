SUBROUTINE NP_CLOSURE
    ! This NP_closure model has 5 tracers: <N>, <P>, <N'>^2, <P'>^2, <N'P'>
    ! Governing functions follow Mandal et al. JPR (2016)
use bio_MOD
use MOD_1D, only: it, nsave
implicit none

integer :: k,j,i
!INPUT PARAMETERS:
real :: tC,par_
!LOCAL VARIABLES of phytoplankton:
real :: PHY1
real :: NO3, PHY
real :: mu0,aI0
real :: QN,Qnmax,Qnmin  ! cell quota related variables
real :: muNet
real :: alphaI,SI,Lno3
real :: muNet1, SI1, Lno31, QN1
real :: rmax_T ! a scratch variable to temporally store phyto. growth rate
real :: theta
real :: Ptot,  CHLt,KN, NPPt

!-----------------------------------------------------------------------
DO k = nlev, 1, -1   

   ! Retrieve current (local) state variable values.
   tC     = Temp(k)
   ! Check whether in the MLD or not
   if (k .lt. N_MLD) then
      par_ = PAR(k)
   else
      par_ = PARavg
   endif
   Varout(oPAR_,k)=par_
   NO3    = Vars(iNO3,    k)
   PHY    = Vars(iPHY(1), k)
   VPHY   = Vars(iVPHY,   k)
   VNO3   = Vars(iVNO3,   k)
   COVNP  = Vars(iCOVNP,  k)

   call PHY_NPZDCONT(NO3,par_,tC,Fe,PMU,muNet,dmuNetdl,d2muNetdl2,d3mudl3,d4mudl4,SI,Lno3, theta, &
    QN,dQNdL,d2QNdL2,dthetadL,d2thetadl2)

   Varout(oTheta(1),k)= theta
   Varout(oQN(1)   ,k)= QN   

   Varout(oSI(1),k)  =SI
   Varout(oLno3(1),k)=Lno3

!=============================================================
! Calculate the community sinking rate of phytoplankton
! Phytoplankton sinking rate at the average size
    w_p    = ScaleTrait(PMU,abs(w_p0),alphaW)  !Positive
! 
! Sinking rate (m) of phytoplankton
    w_pAvg = w_p+0.5*VAR*d2wdl2
! sinking rate must be negative!
    w_pAvg = -min(w_pAvg, Wmax)
!=============================================================
     NPPt          = PHY*(muNet1/QN1+0.5*VAR*d2muNet_QNdl2) 

     Varout(oCHLt,k) = CHLt
     Varout(oPPt, k) = NPPt*12d0 ! Correct the unit to ug C/L
   Endif
!=============================================================
!! Solve ODE functions:
   Zmort2 = MES*MES*mz2*tf_z  !Mortality term for MES

   !Eq. 18, Update PMU:
   PMU = PMU + PMU0  ! Restore to positive values
   PMU1= PMU + dtdays*(VAR*(dmuNetdl-dgdlbar1-dgdlbar2+VTR*d3mudl3) - 3d0*VTR*dmuNetdl)

!All rates have been multiplied by dtdays to get the real rate correponding to the actual time step
   PP_ND= dtdays*RDN*DET*tf_z   
   PP_NZ= (MIC*RES1+MES*RES2)*dtdays     !Sum of MIC and MES to DIN   
   PP_DZ= (MIC*EGES1+MES*EGES2+Zmort2)*dtdays  !Sum of MIC and MES to detritus
   PP_PN= PHY*(muNet+0.5*VAR*(d2muNetdl2+VTR*d4mudl4)-1.5*VTR*d2muNetdl2)

!Update tracers:
   NO3  = NO3 + PP_ND + PP_NZ - PP_PN
   PHY1 = PHY + PP_PN - PP_MIC_P - PP_MES_P
   MIC  = MIC + PP_MIC_P*GGE - PP_MES_MIC
   MES  = MES + (MES*INGES2*GGE - Zmort2)*dtdays
   
   Varout(oNO3,k)      = NO3
   Varout(oPHY(1),k)   = PHY1
   Varout(oPHYt,  k)   = PHY1
   Varout(oZOO,k)      = MIC
   Varout(oZOO2,k)     = MES
   Varout(oDET,k)      = DET1
   Varout(oDETFe,k)    = DETFe
   Varout(ofer,k)      = Fe
   Varout(oPMU,k)      = PMUPHY
   Varout(oVAR,k)      = VARPHY
   Varout(omuNet(1),k) = muNet               !Growth rate of mean size
   Varout(omuAvg,   k) = PP_PN/dtdays/PHY    !Avg. growth rate of the total community
ENDDO
END SUBROUTINE NP_CLOSURE 

! The subroutine only for phytoplankton
subroutine PHY_NPCLOSURE(NO3,PAR_,Temp_,PHY, VNO3, COVNP, muNet,SI,theta,QN, Snp, Chl)
use bio_MOD, only : TEMPBOL, params, mu_Edwards2015 
use bio_MOD, only : iIopt, imu0, iaI0_C, iKN, Ep
use bio_MOD, only : thetamax, thetamin
implicit none
real, intent(in)  :: NO3, PAR_,Temp_,PHY, VNO3, COVNP 
real, intent(out) :: muNet, theta, QN, SI,  Snp, Chl

! muNet: mean growth rate at <N>
real :: mu0hat, mu0hatSI
real :: alphaI, cff,cff1
real :: KN, tf
real, parameter :: Qmin = 0.06, Qmax=0.18

alphaI   = params(ialphaI)
tf       = TEMPBOL(Ep,Temp_)

!The temperature and light dependent component
mu0hat   = tf*params(imu0)
SI       = mu_Edwards2015(PAR_, params(iIopt),mu0hat, alphaI) 
mu0hatSI = mu0hat*SI

KN = params(iKN)
fN = NO3/(NO3 + Kn)  !Nitrogen limitation index

! Phytoplankton growth rate at the mean size:
muNet = mu0hat*SI*fN

! Snp: ensemble mean production
Snp   = PHY*muNet + mu0hatSI*(Kn*COVNP/(Kn+NO3)**2 - Kn*PHY*VNO3/(Kn+NO3)**3)

!N:C ratio at <N>
cff1  = 1d0-Qmin/Qmax
cff   = 1d0-cff1*fN
QN    = Qmin/cff
Q     = 1./QN
!Chl:C ratio at <N>
cff   = (thetamax - thetamin)/PAR_
theta = thetamin+muNet/alphaI*cff   !Unit: gChl/molC

!Chl:N ratio at <N>
eta   = theta/QN
dYdN  = Kn/(NO3 + Kn)**2
dEta_dN  = dYdN*(theta*(1./Qmax - 1./Qmin) + Q*mu0hatSI/alphaI*cff)
d2ChldN2 = -2.*PHY*Kn/(NO3 + Kn)**3 * (thetamin*(1./Qmax - 1./Qmin) + Q*mu0hatSI/alphaI*cff)

Chl = PHY*eta + .5*(2.*COVNP*dEta_dN + VNO3*d2ChldN2)

return
end subroutine

