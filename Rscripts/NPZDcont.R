library(Rcpp)
sourceCpp("~/Working/FlexEFT1D/Rscripts/lambert.cpp")

#Mechaelis-Mention functions and derivatives
MM = function(N, K0, alphaK, L){
   # Half saturation constant for growth at avg. size
   Kn = ScaleTrait(L, K0, alphaK)
   fN = N/(N + Kn) # Nutrient limitation index at avg. size
   dfNdl = -alphaK*Kn*N/(N+Kn)**2
 d2fNdl2 = -alphaK**2*N*Kn*(1/(N+Kn)**2 - 2*Kn/(N+Kn)**3)
 d3fNdl3 = alphaK**3*N*Kn*(2.*N*Kn-(Kn-N)**2)/(Kn+N)**4
 d4fNdl4 = alphaK**4*N*Kn*(11.*Kn*N*(N-Kn)+Kn**3-N**3)/(N+Kn)**5  #Correct
 return(list(fN=fN,dfNdl=dfNdl,d2fNdl2=d2fNdl2,d3fNdl3=d3fNdl3,d4fNdl4=d4fNdl4))
}

#Test NPZDCONT code:
#NPZDCONT(Ep=0.5,mu0=1.2,Temp_=23.44, PAR_=44.7733, NO3=0.028543,
#         PMU=-2.03,
#          Fe=0.0665 )

NPZDCONT <- function(Ep=0.5,mu0=1.5, alphamu=0.3,betamu=-0.03,
                  Temp_=15, PAR_=100,thetamin =0.02,thetamax=0.47,
                    NO3=1, PMU=0, aI0=.04, alphaI=0.07, Q0N=.06, 
                    K0N=.5, alphaK=0.27, dtdays=1, DO_IRON=TRUE,
                   K0Fe=0.04, alphaFe=0.27, Fe=0.1){

tf     = TEMPBOL(Ep,Temp_)
mu0hat = dtdays*tf*mu0*exp(alphamu*PMU + betamu*PMU**2)

dmu0hatdl = mu0hat*(alphamu + 2 * betamu * PMU)

d2mu0hatdl2=2.*mu0hat*betamu+mu0hat*(alphamu+2*betamu*PMU)**2
d3mu0hatdl3=(2.*betamu+(alphamu+2.*betamu*PMU)**2)*dmu0hatdl+4.*betamu*mu0hat*(alphamu+2.*betamu*PMU)
d4mu0hatdl4=dmu0hatdl*8.*betamu*(alphamu+2.*betamu*PMU)+(2.*betamu+(alphamu+2.*betamu*PMU)**2)*d2mu0hatdl2+8.*betamu**2*mu0hat

# Initial slope of P-I curve
aI=ScaleTrait(PMU, aI0, alphaI) * dtdays

#The light limitation index (SI)
cff=exp(-aI*PAR_/mu0hat)
SI =1-cff

daI_mu0hatdl=aI*(alphaI-alphamu-2*betamu*PMU)/mu0hat

dSIdl   =cff*PAR_*daI_mu0hatdl
mu0hatSI=mu0hat*SI

dmu0hatSIdl    = mu0hatSI*(alphamu+2*betamu*PMU) + mu0hat*dSIdl
dmu0hat_aIdl   = (dmu0hatdl - alphaI*mu0hat)/aI
d2mu0hat_aIdl2 = mu0hat/aI*2*betamu+(alphamu-alphaI+2*betamu*PMU)*dmu0hat_aIdl

daI_mu0hat2dl  = daI_mu0hatdl/mu0hat  - aI/mu0hat**3*dmu0hatdl #Correct

d2aI_mu0hatdl2 = alphaI*daI_mu0hatdl - aI/mu0hat**2*d2mu0hatdl2 - daI_mu0hat2dl*dmu0hatdl

#d2aI_mu0hatdl2 = -2*betamu*aI/mu0hat-(2*betamu*PMU+(alphamu-alphaI))*daI_mu0hatdl #(same with above)
  
d3aI_mu0hatdl3 = d2aI_mu0hatdl2*(alphaI-alphamu-2.*betamu*PMU)-4.*betamu*daI_mu0hatdl

d2SIdl2 = PAR_*((1.-SI)*d2aI_mu0hatdl2 - daI_mu0hatdl * dSIdl)

d2mu0hatSIdl2 = (alphamu+2*betamu*PMU)*dmu0hatSIdl + 2*betamu*mu0hatSI + mu0hat*(alphamu+2*betamu*PMU)*dSIdl + mu0hat*d2SIdl2  #Correct

d3SIdl3 = PAR_*(-2.*dSIdl*d2aI_mu0hatdl2 - d2SIdl2*daI_mu0hatdl + (1.-SI)*d3aI_mu0hatdl3)  #Correct

daI_mu0hat3dl  = daI_mu0hat2dl/mu0hat - aI*dmu0hatdl/mu0hat**4 #Correct

d2aI_mu0hat2dl2= d2aI_mu0hatdl2/mu0hat - daI_mu0hatdl/mu0hat**2*dmu0hatdl - daI_mu0hat3dl * dmu0hatdl - aI/mu0hat**3*d2mu0hatdl2 #Correct 

d3aI_mu0hatdl3 = alphaI*d2aI_mu0hatdl2 - aI/mu0hat**2*d3mu0hatdl3 -2.*daI_mu0hat2dl*d2mu0hatdl2 -d2aI_mu0hat2dl2*dmu0hatdl #Correct

daI_mu0hat4dl  = daI_mu0hat3dl/mu0hat - aI/mu0hat**5*dmu0hatdl  #Correct

d2aI_mu0hat3dl2= d2aI_mu0hat2dl2/mu0hat - daI_mu0hat2dl/mu0hat**2*dmu0hatdl - aI/mu0hat**4*d2mu0hatdl2 - daI_mu0hat4dl*dmu0hatdl   #Correct

d3aI_mu0hat2dl3= (d3aI_mu0hatdl3*mu0hat - daI_mu0hatdl*d2mu0hatdl2)/mu0hat**2 -2./mu0hat**3*(d2aI_mu0hatdl2*mu0hat - daI_mu0hatdl*dmu0hatdl)*dmu0hatdl - (aI/mu0hat**3*d3mu0hatdl3+2.*daI_mu0hat3dl*d2mu0hatdl2+d2aI_mu0hat3dl2*dmu0hatdl)  #Correct

d4aI_mu0hatdl4 = alphaI*d3aI_mu0hatdl3 - (aI/mu0hat**2*d4mu0hatdl4 + 3.* daI_mu0hat2dl*d3mu0hatdl3+3.*d2aI_mu0hat2dl2 * d2mu0hatdl2)- d3aI_mu0hat2dl3*dmu0hatdl  #Correct

d4SIdl4 = PAR_*((1.-SI)*d4aI_mu0hatdl4 - 3.*dSIdl*d3aI_mu0hatdl3 - 3.*d2aI_mu0hatdl2*d2SIdl2-daI_mu0hatdl*d3SIdl3)   #Correct

d3muIhatdl3 = d2mu0hatSIdl2*(alphamu + 2.*betamu*PMU) + dmu0hatSIdl*4.*betamu + d2mu0hatdl2*dSIdl + 2.*dmu0hatdl*d2SIdl2 + mu0hat*d3SIdl3  #Correct

d4muIhatdl4 = mu0hat*d4SIdl4 + 4.*dmu0hatdl*d3SIdl3 + 6.*d2mu0hatdl2*d2SIdl2 + 4.*dSIdl*d3mu0hatdl3+SI*d4mu0hatdl4  #Correct

cff    = MM(NO3,K0N,alphaK,PMU)
Kn     = ScaleTrait(PMU, K0N, alphaK)
fN     = cff$fN
dfNdl  = cff$dfNdl
d2fNdl2= cff$d2fNdl2
d3fNdl3= cff$d3fNdl3
d4fNdl4= cff$d4fNdl4
# Phytoplankton growth rate at the mean size:
muNet = mu0hat*SI*fN
dmudl = dmu0hatSIdl*fN - mu0hatSI*alphaK*Kn*NO3/(NO3 + Kn)**2
d2mudl2=2*dmu0hatSIdl*dfNdl+d2mu0hatSIdl2*fN+mu0hatSI*d2fNdl2 
d3mudl3=3.*(d2mu0hatSIdl2*dfNdl+dmu0hatSIdl*d2fNdl2) +fN*d3muIhatdl3 + mu0hatSI*d3fNdl3 #Correct
d4mudl4=4.*d3muIhatdl3*dfNdl+6.*d2mu0hatSIdl2*d2fNdl2+4.*dmu0hatSIdl*d3fNdl3 + fN*d4muIhatdl4 + mu0hatSI*d4fNdl4  #Correct

#Add iron limitation:
if (DO_IRON){
   cff     = MM(Fe,K0Fe,alphaFe,PMU)
   KFe     = ScaleTrait(PMU, K0Fe, alphaFe)
   fFe     = cff$fN
   dfFedl  = cff$dfNdl
   d2fFedl2= cff$d2fNdl2
   d3fFedl3= cff$d3fNdl3
   d4fFedl4= cff$d4fNdl4
   muNet_  = muNet*fFe
   dmudl_  = muNet*dfFedl + dmudl*fFe
   d2mudl2_= muNet*d2fFedl2 + 2.*dmudl*dfFedl + fFe*d2mudl2
   d3mudl3_= muNet*d3fFedl3 + 3.*(dmudl*d2fFedl2+d2mudl2*dfFedl)+fFe*d3mudl3
   d4mudl4 = muNet*d4fFedl4 + 4.*dmudl*d3fFedl3 + 6.*d2mudl2*d2fFedl2 + 4.*d3mudl3*dfFedl + fFe*d4mudl4  
   muNet   = muNet_
   dmudl   = dmudl_
   d2mudl2 = d2mudl2_
   d3mudl3 = d3mudl3_
}  #All correct

Qmin=Q0N
Qmax=3*Qmin

#N:C ratio at avg. size
 cff1   = 1.-Qmin/Qmax
 QN     = Qmin/(1.-cff1*fN)


cff     = 1.-cff1*fN
dcffdl  = -cff1*dfNdl
dQNdL   = cff1/cff**2 * dfNdl * Qmin
d2QNdL2 = Qmin*cff1*(d2fNdl2/cff**2 - 2./cff**3*dfNdl*dcffdl)  #Correct

cff     = (thetamax - thetamin)/PAR_

theta   = thetamin+muNet/aI*cff   #Unit: gChl/molC

dthetadl = cff * dY_Xdl(muNet, aI, dmudl, aI*alphaI)  #Correct

d2thetadl2 = cff*d2Y_Xdl2(muNet,aI,dmudl,aI*alphaI, d2mudl2, aI*alphaI**2) #Correct

 return(list(muNet=muNet, dmudl=dmudl, d2mu=d2mudl2, d3mu=d3mudl3,
            d4mu  =d4mudl4,  fN=fN,    fFe =fFe, SI = SI,
             theta=theta,    QN=QN,    dQNdL=dQNdL,
          dthetadl=dthetadl, d2QNdL2=d2QNdL2, d2thetadl2=d2thetadl2))
}

#Test the growth ~ temperature function:
Mu = function(Topt = 300, Temp = 293, mu0 = 5.7, k = -0.11, 
              Ea = 1, Eh = 5, kb = 8.62E-5, T0 = 298.15){
     T_Zero = 273.15
     if(Topt < 100) Topt=Topt+T_Zero
     if(Temp < 100) Temp=Temp+T_Zero
     a =Ea/kb
     b =Eh/kb
     c =Ea/(Eh-Ea)
     Y =1.+c*exp(b*(1./Topt-1./Temp))
     dYdx = (1.-Y)*b/Topt**2
     mu   =mu0*exp(k*(Topt-T0))*exp(a*(1./T0 - 1./Temp))/Y
     dmudl=mu*(k + (Y-1)*b/Topt**2*Y)
     d2mudl2=dmudl*(k+b*(Y-1.)/Topt**2/Y) - b*mu*(Y-1.)*(2.+b/Y/Topt)/Topt**3/Y

     d3mudl3=d2mudl2*(k+b/Topt**2+dmudl/mu+2./Topt) - dmudl**3/mu**2-b*mu*(Y-1)/Topt**2/Y*(b**2*(Y-1)/Topt**4/Y**2 -2./Topt**2 - 2.*b/Topt**3/Y) - k*dmudl*(2./Topt + b/Y/Topt**2)  #Correct

     A1=-2.*b/Topt**3 + d2mudl2/mu - dmudl**2/mu**2 - 2./Topt**2
     dAdx = d3mudl3*(k+b/Topt**2+dmudl/mu+2./Topt) + d2mudl2*A1
     dBdx = -2./mu**3*dmudl**4 + 3./mu**2*dmudl**2*d2mudl2
     C1a  = (2.-Y)/Topt**4/Y**3*dYdx
     C1b  = -(Topt**3*dYdx + Y*3.*Topt**2)/(Topt**3*Y)**2
     C1   = b**2*C1a + 1./Topt**3 - 2.*b*C1b
     D    = d2mudl2 - k*dmudl
     dCdx = b*mu*(Y-1.)/Topt**2/Y*C1 + D*(b**2*(Y-1)/Topt**4/Y**2-2./Topt**2-2.*b/Topt**3/Y)
     cff  = -(2./Topt**3/Y + dYdx/Topt**2/Y**2)
     dDdx = k*(dmudl*(-2./Topt**2 + b*cff) + (2./Topt+b/Y/Topt**2)*d2mudl2)
     d4mudl4 = dAdx - (dBdx+dCdx+dDdx)
     return(list(mu=mu,dmudl=dmudl,d2mudl2=d2mudl2,d3mudl3=d3mudl3, d4mudl4 = d4mudl4))
}

#Three trait function:
mu3 = function(x, Z, I=200, NO3=1, t=15, w=10, K0N=0.2, alphaK=.27, L=0, k=-.47, a0=.149, b=.21,mu0=.21, alphamu=.2, betamu=0, Ep=.0633, T0=15){
    #x: ln optimal light (umol photons m-2 s-1)
    #Z: optimal temperature (ºC)
    #Light term
    umax = mu0*exp(b*x)
    alpha= a0*exp(k*x)
    Iopt = exp(x)
    H  = 1./(I/(alpha*Iopt**2) + 1./umax -2./(alpha*Iopt) + 1./(I*alpha))
    #Note here suboptimal is optimal (even if I < Iopt, because umax is higher, mu is also higher)

    #Size (nutrient) term (L: lnvolume)
    hat    = exp(alphamu*L + betamu*L**2)
    dhatdl = hat*(alphamu + 2.*betamu * L)
    d2hatdl2 = 2.*hat*betamu+hat*(alphamu+2.*betamu*L)**2
    cff    = MM(NO3,K0N,alphaK,L)
    Kn     = ScaleTrait(L, K0N, alphaK)
    fN     = cff$fN
    dfNdl  = cff$dfNdl
    d2fNdl2= cff$d2fNdl2

    G      = hat * fN
    #Calculate the derivatives of G
    dGdL   = hat*dfNdl + dhatdl*fN
    d2GdL2 = hat*d2fNdl2 + 2.*dfNdl*dhatdl + fN*d2hatdl2

    #Calculate the derivatives of H
    Y      = 1./H
    dYdx   = -(k+2.)*I/(alpha*Iopt**2) - b/umax +2.*(k+1.)/(alpha*Iopt) - k/(I*alpha)
    d2Ydx2 = (k+2.)**2*I/(alpha*Iopt**2) + b**2/umax -2.*(k+1.)**2/(alpha*Iopt) + k**2/(I*alpha)
    dHdx   = - H*H*dYdx
    d2Hdx2 = (2*H**3)*(dYdx)**2 - d2Ydx2*H**2

 #Temperature term
    h       = exp(Ep*(t - T0))
    P       = 1.-((t-Z)/w)**2
    f       = h * P
    dPdZ    = 2.*(t-Z)/w**2
    d2PdZ2  = -2./w**2
    dfdZ    = h*dPdZ 
    d2fdZ2  = h*d2PdZ2 

    Mu      = H * G * f
    dMudx   = G*f*dHdx
    d2Mudx2 = G*f*d2Hdx2
    dMudL   = f*H*dGdL
    d2MudL2 = f*H*d2GdL2
    dMudZ   = H*G*dfdZ
    d2MudZ2 = H*G*d2fdZ2
    return(list(mu=Mu, dmudx=dMudx, d2mudx=d2Mudx2, dmudL=dMudL, d2mudL2=d2MudL2, dmudZ=dMudZ, d2mudZ2=d2MudZ2))
}

#Check the shape of mu ~ size
#mu = function(x, mu0hat = 1.2, alphamu = 0.01, betamu = -0.01) mu0hat*exp(alphamu*x+betamu*x**2)
#curve(mu, from = -3, to = 20)
#
#rhol =function(x,a=0,b=-0.2,c=0) exp(a*x*x+b*x+c) 
#curve(rhol, from =-3, to =20)
#
#
#Lmin <- log(pi/6*0.6**3)
#Lmax <- log(pi/6*40**3)
#L    <- seq(Lmin, Lmax, length.out = 40)
#
#NO3  <- c(0.1, 1, 4)
#
#pdffile <- '~/Working/FlexEFT1D/DRAM_0.9/NPZDcont/growth_derivatives.pdf'
#pdf(pdffile, width = 3, height = 8, paper= 'a4')
# op <- par(font.lab = 1,
#             family ="serif",
#             mar    = c(4,4,1.5,2),
#             mgp    = c(2.3,1,0),
#             mfrow  = c(3,1),
#             oma    = c(4,4,0,0)) 
#   
#
#for (i in 1:length(NO3)){
#   
#   mu=sapply(1:length(L), function(l)NPZDCONT(PMU=L[l], NO3=NO3[i])$muNet) 
#    if (i == 1){
#       plot(L, mu, type = 'l', ylim = c(0, 3), 
#            xlab = bquote('Ln volume ('*µm^3*')'), 
#            ylab = bquote('µ ( '*d^-1*')'),
#             col = 1 + i)
#    } else{
#       points(L, mu, type = 'l', col = 1 + i)
#    }
#}
#legend('topleft', lty = 1, legend = paste0('NO3=',NO3),col=2:(length(NO3)+1))
#
#for (i in length(NO3):1){
#   
#    dmudl <- sapply(1:length(L),function(l)NPZDCONT(PMU=L[l],NO3=NO3[i])$dmudl) 
#    if (i == length(NO3)){
#       plot(L, dmudl, type = 'l', ylim = c(-0.2,0.2),
#            xlab = bquote('Ln volume ('*µm^3*')'), 
#            ylab = bquote('dµ/dL ( '*'(d ln'*µm^3*')'^-1*')'),
#             col = 1 + i)
#    } else{
#       points(L, dmudl, type = 'l', col = 1 + i)
#    }
#}
#abline(0,0)
#
##d2µdl2
#for (i in length(NO3):1){
#   
#    d2udl2 <- sapply(1:length(L), function(l)NPZDCONT(PMU=L[l], NO3=NO3[i])$d2mu) 
#    if (i == length(NO3)){
#       plot(L, d2udl2, type = 'l', ylim = c(-0.05,0.05),
#            xlab = bquote('Ln volume ('*µm^3*')'), 
#            ylab = bquote('d2µ/dL2 ( '*'(d ln'*µm^3*')'^-2*')'),
#             col = 1 + i)
#    } else{
#       points(L, d2udl2, type = 'l', col = 1 + i)
#    }
#}
#abline(0,0)
#
#dev.off()

