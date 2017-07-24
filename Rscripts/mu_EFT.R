#The Pahlow model that returns to mu, theta, and QN
library(Rcpp)
sourceCpp("~/Working/FlexEFT1D/Rscripts/lambert.cpp")

#Parameters very important!
mu_fix <-function(NO3=50,par_  = 100, Temp = 15, mu0=1.2, KN0 = 3.486, 
                  Q0N  =.06,Ep = 0.6, aI0C= 0.066, thetamax = 0.42){
 #Try to be as consistent with flexible model as possible
    thetamin <- 0.02
    rmax  <- mu0 * TEMPBOL(Ep, Temp) 
    SI    <- 1 - exp(-aI0C * par_/rmax)
    gr    <- rmax * NO3/(NO3 + KN0) * SI
    Qmax  <- 3*Q0N
    QN    <- Q0N/(1-(1-Q0N/Qmax)*NO3/(NO3+KN0))
  theta   <- thetamin+gr/par_/aI0C*(thetamax-thetamin)   #Unit: gChl/molC
    return(list(mu=gr, QN=QN, Theta=theta))
}

#Test the effect of aI0Chl on theta:
#f1 <- function(I)mu_fix(par_=I)$Theta
#curve(f1, from=0.1,to=200,ylim=c(0,0.6))
#f2 <- function(I)mu_fix(par_=I, aI0Chl=.2,aI0C=0.1)$Theta
#curve(f2, from=0.1,to=200,add=T,col=2)


mu_EFT <- function(temp = 15,  par_=400, NO3=100, aI = 6.103E-01, Q0N = 0.07561,
                   A0N  = 9.139E-02, Ep = 0.5, mu0= 2.541,   V0N= 2.541,
                   RMchl0 = .1,zetaChl=.8,  zetaN=0.6,  KN0 = 1,
                   nutrient_uptake=2,KFe=.2, DFe=0.4, do_IRON=F){

# Phytoplankton section:
   V0N  = mu0
   tf_p = TEMPBOL(Ep,temp)
   Qs   = Q0N/2
 mu0hat = tf_p*mu0
  V0hat = tf_p*V0N

  if (do_IRON){
      mu0hat = mu0hat * DFe /(DFe + KFe)
      V0hat  =  V0hat * DFe /(DFe + KFe)
      aI     =     aI * DFe /(DFe + KFe)  #this does not affect thetaHat too much
  }
  # Cost of photosynthesis
  RMchl  = tf_p*RMchl0
  # Threshold irradiance and RMchl is set temperature dependent
  I_zero = zetaChl*RMchl/aI  
  
#Define VNhat: the saturation function of ambient nutrient concentration
  if (nutrient_uptake == 1){  
  # case 1: Classic Michaelis Menton 
  # Potential maximal nutrient-limited uptake
      # Half-saturation constant of nitrate uptake
      VNhat  = V0hat*NO3/(NO3 + KN0)

 # case 2: optimal uptake based on Pahlow (2005) and Smith et al. (2009)
  }else if(nutrient_uptake == 2){
      A0hat = tf_p*A0N
      #Define fA
      fA    = 1/(1 + sqrt(A0hat*NO3/V0hat) ) 
      VNhat = (1-fA)*V0hat*fA*A0hat*NO3/((1-fA)*V0hat+fA*A0hat*NO3) 
  }else{
      stop('Error: Incorrect option for nutrient uptake')
  }

# Calculate thetahat (optimal g Chl/mol C for the chloroplast under nutrient replete conditions)
# Only calculate within the euphotic zone, otherwise many numerical problems.
  if( par_ > I_zero ) {
    larg1 = exp(1 + aI*par_/(mu0hat*zetaChl))
    larg  = (1 + RMchl/mu0hat)*larg1   
    
     W_Y  = WAPR(larg,0,0)

   ThetaHat = 1/zetaChl + (1- W_Y)*mu0hat/(aI * par_)

    # Effect of light limitation: the same aI is different for flexible and inflexible models.
    # More light limited for EFT model if using the same aI

    SI = 1 - max(exp(-aI * par_ * ThetaHat/mu0hat),0.)
# Light dependent growth rate 
# (needs to take into account the cost of dark and light-dependent chl maintenance)
    mu0hatSI = mu0hat*SI #Gross specific carbon uptake (photosynthesis)
    muIhat   = mu0hatSI-(mu0hatSI+RMchl)*zetaChl*ThetaHat # Net specific carbon uptake
    muIhat   = max(muIhat,0.)
   
    ZINT   = Qs*(muIhat/VNhat + zetaN)
           
    fV = (-1 + sqrt(1 + 1/ZINT))*Qs*muIhat/VNhat
    
    }else{
    # Under the conditions of no light:
       ThetaHat         = .02  #  thetamin
       ZINT             = Qs*zetaN
       fV               = 1E-2
       muIhat           = 0.
    }
     # Nutrient limitation index:
     Lno3 =1/(1 + sqrt(1 +1/ZINT)) 
     # Optimal nutrient quota:
     QN   = Qs/Lno3
    
     if (par_ > I_zero){
   #Net growth rate (d-1) of phytoplankton at the average size
       muNet = muIhat*(1-2.*Lno3)  #Droop model
     }else{
       muNet = 0.
     }

   #Save net growth rate
    muNet = muNet

   #Chl:C ratio [ g chl / mol C ] of the whole cell
    Theta = ThetaHat*(1-fV-Qs/QN)
    return(list(mu=muNet,Theta=Theta,QN=QN))
}

#
#A function that generates the concentrations within each size class
dB <- function(B, mean, var, mu=seq(1,20,length.out=100),N=10000){
   M  <- length(mu)
   X  <- rnorm(N, mean=mean, sd = sqrt(var))
   NewX <- numeric(M)
   for (i in 1:M){
       #Count the number of X within one interval:
       NewX[i] <- length(which(X >= mu[i] & X < mu[i+1]))
   }
   NewX <- NewX/N*B
   return(NewX)
}

