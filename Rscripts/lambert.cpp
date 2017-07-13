#include <Rcpp.h>
#include <math.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double TEMPBOL(double Ea, double tC){
  //DESCRIPTION:
  //The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
  // tC: in situ temperature
  // Tr: reference temperature
  //
  //INPUT PARAMETERS:
  // boltzman constant constant [ eV /K ]
  double kb = 8.62E-5;
  double Tr = 15.;
  double Y;
  Y = exp(-(Ea/kb)*(1/(273.15 + tC)-1/(273.15 + Tr)));
  return Y;
} 

// [[Rcpp::export]]
double ScaleTrait(double logsize, double star, double alpha){ 

   // Calculate the size-scaled value of a trait
   // for the given log (natural, base e) of cell volume as pi/6*ESD**3 (micrometers). 
   double Y;
   
   Y = star * exp( alpha * logsize );
   
   return Y;
}

// [[Rcpp::export]]
double dY_Xdl(double Y, double X, double dYdl, double dXdl){ 
   double cff;
   cff = dYdl/X - Y/pow(X,2)*dXdl;
   return cff;
}

// [[Rcpp::export]]
double d2Y_Xdl2(double Y, double X, double dYdl, double dXdl, double d2Ydl2, double d2Xdl2){ 
   double cff;
   cff = d2Ydl2/X - 2.*dYdl*dXdl/pow(X,2) - Y/pow(X,2)*d2Xdl2 + 2.*Y/pow(X,3)*pow(dXdl,2);
   return cff;
}

// [[Rcpp::export]]
double ALM(double L, double A, double B){
      // L: cell size (ESD)
      // A: Rate at 1 um3
      // B: allometric exponent
      double Vol;
      double pi=3.1415926535;
      double mu;
      
      Vol = log(pi/6.*pow(L,3) );
      mu  = A * exp(Vol * B);
      return mu;
} 

// [[Rcpp::export]]
double WAPR(double X, double NB, double L) {
  int NN     = 6;
  int NBITS  = 23;
  int NITER  = 1;
 
  double EM    = -0.367879441171442;    // -EXP(-1)
  double EM9   = -1.234098040866796E-4;       // -EXP(-9)
  double C13   = 1.E0/3.E0;
  double C23;
  double EM2;
  double E12;
  double TB;
  double TB2;
  double X0    = 0.0350769390096679055; // TB**(1/6)*0.5E0
  double X1    =  0.302120119432784731; //  #(1 - 17*TB**(2/7))*EM
  double AN22  =  3.6E2/83.E0;
  double AN11  =  135./83.E0;
  double AN3   =  8.E0/3.E0;
  double AN4   =  135.E0/83.E0;
  double AN5   =  166.E0/39.E0;
  double AN6   =  3167.E0/3549.E0;
  double S2;
  double S21;
  double S22;
  double S23;
  double DELX;
  double WAP;
  double RETA;
  double XX;
  double AN2;
  double ZL;
  double ZN;
  double TEMP;
  double TEMP2;
  double ETA;
  double T;
  double TS;


  C23  =  2.E0*C13;
  EM2  =  2.E0/EM;
  E12  =  -EM2;
  TB   =  pow(.5E0,NBITS);
  TB2  =  pow(.5E0, (NBITS/2));
  S2   =  sqrt(2);
  S21  =  2.E0*S2-3.E0;
  S22  =  4.E0-3.E0*S2;
  S23  =  S2-2.E0;

  if (L == 1) {
    DELX = X;
    if (DELX < 0) {
      WAP = 1./0.;
    }
     XX=X+EM;
  } else {
    if(X < EM) {
      WAP = 1./0.;
    }
    if(X == EM){
      WAP = -1.E0;
    }
  XX  =X;
  DELX=XX-EM;
  }

  if (NB == 0) {

  /*  Calculations for Wp */

    if(std::abs(XX) <= X0) {
      WAP  = XX/(1.E0+XX/(1.E0+XX/(2.E0+XX/(.6E0+.34E0*XX))));
    } else if (XX <= X1) {
      RETA=  sqrt(E12*DELX);
      WAP =  RETA/(1.E0+RETA/(3.E0+RETA/(RETA/(AN4+RETA/(RETA*AN6+AN5))+AN3)))-1.E0;
    } else if (XX <= 2.E1) {
     RETA = S2*sqrt(1.E0-XX/EM);
      AN2 = 4.612634277343749E0*sqrt(sqrt(RETA+1.09556884765625E0));
      WAP = RETA/(1.E0+RETA/(3.E0+(S21*AN2+S22)*RETA/(S23*(AN2+RETA))))-1.E0;
    } else{
      ZL  =  log(XX);
      WAP =  log(XX/log(XX/pow(ZL,exp(-1.124491989777808E0/(.4225028202459761E0+ZL)))));
    }
  } else{
  //        Calculations for Wm
    if (XX >= 0.E0) {
      WAP  = -1./0.;
    }
    if (XX<=X1) {
      RETA = sqrt(E12*DELX);
      WAP  = RETA/(RETA/(3.E0+RETA/(RETA/(AN4+RETA/(RETA*AN6-AN5))-AN3))-1.E0)-1.E0;
    }else if(XX <= EM9) {
      ZL  =  log(-XX);
      T   =  -1.E0-ZL;
      TS  =  sqrt(T);
      WAP =  ZL-(2.E0*TS)/(S2+(C13-T/(2.7E2+TS*127.0471381349219E0))*TS);
    }else {
      ZL  =  log(-XX);
      ETA =  2.E0-EM2*XX;
      WAP =  log(XX/log(-XX/((1.E0-.5043921323068457E0*(ZL+1.E0))*(sqrt(ETA)+ETA/3.E0)+1.E0)));
    }
  }  

ZN    = log(XX/WAP)-WAP;
TEMP  = 1.E0+WAP;
TEMP2 = TEMP+C23*ZN;
TEMP2 = 2.E0*TEMP*TEMP2;
  WAP = WAP*(1.E0+(ZN/TEMP)*(TEMP2-ZN)/(TEMP2-2.E0*ZN));
  return WAP;

}

//Mechaelis-Mention functions and derivatives
// [[Rcpp::export]]
double MM(N, K0, alphaK, L){
   double Kn, fN, dfNdl, d2fNdl2, d3fNdl3, d4fNdl4;
   
   // Half saturation constant for growth at avg. size
   Kn    = ScaleTrait(L, K0, alphaK);
   fN    = N/(N + Kn); // Nutrient limitation index at avg. size
   dfNdl = -alphaK*Kn*N/pow((N+Kn),2);
 d2fNdl2 = -pow(alphaK,2)*N*Kn*(1./pow((N+Kn),2) - 2.*Kn/pow((N+Kn),3));
 d3fNdl3 = pow(alphaK,3)*N*Kn*(2.*N*Kn-pow((Kn-N),2))/pow((Kn+N),4);
 d4fNdl4 = pow(alphaK,4)*N*Kn*(11.*Kn*N*(N-Kn)+pow(Kn,3)-pow(N,3))/pow((N+Kn),5); 
 return(list(fN=fN,dfNdl=dfNdl,d2fNdl2=d2fNdl2,d3fNdl3=d3fNdl3,d4fNdl4=d4fNdl4))
}


// [[Rcpp::export]]
double NPZDcontFe(double Temp_, double PAR_, double NO3, double Fe, double PMU)

  double Ep      = 0.5;
  double mu0     = 1.5;
  double alphamu = 0.3;
  double betamu  = -0.03;
  double thetamin= 0.02;
  double thetamax= 0.47;
  double aI0     = 0.04;
  double alphaI  = 0.07;
  double Q0N     = 0.06;
  double K0N     = 0.5;
  double alphaK  = 0.27;
  double K0Fe    = 0.04;
  double alphaFe = 0.27;
  
  double mu0hat;
  double dmu0hatdl;
  double d2mu0hatdl2;
  double d3mu0hatdl3;
  double d4mu0hatdl4;
  double aI;
  double SI;
  double cff;
  double cff1;
  double daI_mu0hatdl;
  double dSIdl, mu0hatSI;
  double dmu0hatSIdl, dmu0hat_aIdl,d2mu0hat_aIdl2;
  double d2aI_mu0hatdl2,d3aI_mu0hatdl3,d2SIdl2,d2mu0hatSIdl2;
  double d3SIdl3,daI_mu0hat3dl,d2aI_mu0hat2dl2,d4SIdl4 ;
  double daI_mu0hat4dl,d2aI_mu0hat3dl2,d3aI_mu0hat2dl3 ;
  double d3muIhatdl3, d4muIhatdl4, d4aI_mu0hatdl4;
  double Kn,dfNdl,d2fNdl2,d3fNdl3, d4fNdl4,tf;
  double fFe, dfFedl, d2fFedl2, d3fFedl3, d4fFedl4;
  double KFe, dcffdl, cff2;


  tf     = TEMPBOL(Ep,Temp_);
  mu0hat = tf*mu0*exp(alphamu*PMU + betamu*pow(PMU,2));
  cff1   = alphamu+2.*betamu*PMU;
  
  dmu0hatdl  = mu0hat*cff1;
  d2mu0hatdl2=2.*mu0hat*betamu+mu0hat*pow(cff1,2);
  d3mu0hatdl3=(2.*betamu+pow(cff1,2))*dmu0hatdl+4.*betamu*mu0hat*cff1;
  d4mu0hatdl4=dmu0hatdl*8.*betamu*cff1+(2.*betamu+pow(cff1,2))*d2mu0hatdl2+8.*pow(betamu,2)*mu0hat;
  
  // Initial slope of P-I curve
  aI=ScaleTrait(PMU, aI0, alphaI);
  
  // The light limitation index (SI)
  cff            = exp(-aI*PAR_/mu0hat);
  SI             = 1.-cff;
  daI_mu0hatdl   = aI*(alphaI-cff1)/mu0hat;
  dSIdl          = cff*PAR_*daI_mu0hatdl;
  mu0hatSI       = mu0hat*SI;
  dmu0hatSIdl    = mu0hatSI*cff1 + mu0hat*dSIdl;
  dmu0hat_aIdl   = (dmu0hatdl - alphaI*mu0hat)/aI;
  d2mu0hat_aIdl2 = mu0hat/aI*2.*betamu+(alphamu-alphaI+2.*betamu*PMU)*dmu0hat_aIdl;
  daI_mu0hat2dl  = daI_mu0hatdl/mu0hat - aI/pow(mu0hat,3)*dmu0hatdl;
  d2aI_mu0hatdl2 = alphaI*daI_mu0hatdl - aI/pow(mu0hat,2)*d2mu0hatdl2 - daI_mu0hat2dl*dmu0hatdl;
  
  d3aI_mu0hatdl3 = d2aI_mu0hatdl2*(alphaI-cff1)-4.*betamu*daI_mu0hatdl;
  
  d2SIdl2 = PAR_*((1.-SI)*d2aI_mu0hatdl2 - daI_mu0hatdl * dSIdl);
  
  d2mu0hatSIdl2 = (alphamu+2*betamu*PMU)*dmu0hatSIdl + 2*betamu*mu0hatSI + mu0hat*(alphamu+2*betamu*PMU)*dSIdl + mu0hat*d2SIdl2  #Correct
  
  d3SIdl3 = PAR_*(-2.*dSIdl*d2aI_mu0hatdl2 - d2SIdl2*daI_mu0hatdl + (1.-SI)*d3aI_mu0hatdl3)  #Correct
  
  daI_mu0hat3dl  = daI_mu0hat2dl/mu0hat - aI*dmu0hatdl/mu0hat**4 #Correct
  
  d2aI_mu0hat2dl2= d2aI_mu0hatdl2/mu0hat - daI_mu0hatdl/mu0hat**2*dmu0hatdl - daI_mu0hat3dl * dmu0hatdl - aI/mu0hat**3*d2mu0hatdl2 #Correct 
  
  d3aI_mu0hatdl3 = alphaI*d2aI_mu0hatdl2 - aI/mu0hat**2*d3mu0hatdl3 -2.*daI_mu0hat2dl*d2mu0hatdl2 -d2aI_mu0hat2dl2*dmu0hatdl #Correct
  
  daI_mu0hat4dl  = daI_mu0hat3dl/mu0hat - aI/mu0hat**5*dmu0hatdl  #Correct
  
  d2aI_mu0hat3dl2= d2aI_mu0hat2dl2/mu0hat - daI_mu0hat2dl/mu0hat**2*dmu0hatdl - aI/mu0hat**4*d2mu0hatdl2 - daI_mu0hat4dl*dmu0hatdl; 
  
  d3aI_mu0hat2dl3= (d3aI_mu0hatdl3*mu0hat - daI_mu0hatdl*d2mu0hatdl2)/mu0hat**2 -2./mu0hat**3*(d2aI_mu0hatdl2*mu0hat - daI_mu0hatdl*dmu0hatdl)*dmu0hatdl - (aI/mu0hat**3*d3mu0hatdl3+2.*daI_mu0hat3dl*d2mu0hatdl2+d2aI_mu0hat3dl2*dmu0hatdl)  #Correct
  
  d4aI_mu0hatdl4 = alphaI*d3aI_mu0hatdl3 - (aI/pow(mu0hat,2)*d4mu0hatdl4 + 3.* daI_mu0hat2dl*d3mu0hatdl3+3.*d2aI_mu0hat2dl2 * d2mu0hatdl2)- d3aI_mu0hat2dl3*dmu0hatdl;
  
  d4SIdl4 = PAR_*((1.-SI)*d4aI_mu0hatdl4 - 3.*dSIdl*d3aI_mu0hatdl3 - 3.*d2aI_mu0hatdl2*d2SIdl2-daI_mu0hatdl*d3SIdl3);
  
  d3muIhatdl3 = d2mu0hatSIdl2*(alphamu + 2.*betamu*PMU) + dmu0hatSIdl*4.*betamu + d2mu0hatdl2*dSIdl + 2.*dmu0hatdl*d2SIdl2 + mu0hat*d3SIdl3;  //Correct
  
  d4muIhatdl4 = mu0hat*d4SIdl4 + 4.*dmu0hatdl*d3SIdl3 + 6.*d2mu0hatdl2*d2SIdl2 + 4.*dSIdl*d3mu0hatdl3+SI*d4mu0hatdl4;  //Correct
  
  cff    = MM(NO3,K0N,alphaK,PMU)
  Kn     = ScaleTrait(PMU, K0N, alphaK)
  fN     = cff$fN
  dfNdl  = cff$dfNdl
  d2fNdl2= cff$d2fNdl2
  d3fNdl3= cff$d3fNdl3
  d4fNdl4= cff$d4fNdl4
  // Phytoplankton growth rate at the mean size:
  muNet = mu0hat*SI*fN
  dmudl = dmu0hatSIdl*fN - mu0hatSI*alphaK*Kn*NO3/pow((NO3 + Kn),2);
  d2mudl2=2*dmu0hatSIdl*dfNdl+d2mu0hatSIdl2*fN+mu0hatSI*d2fNdl2; 
  d3mudl3=3.*(d2mu0hatSIdl2*dfNdl+dmu0hatSIdl*d2fNdl2) +fN*d3muIhatdl3 + mu0hatSI*d3fNdl3;
  d4mudl4=4.*d3muIhatdl3*dfNdl+6.*d2mu0hatSIdl2*d2fNdl2+4.*dmu0hatSIdl*d3fNdl3 + fN*d4muIhatdl4 + mu0hatSI*d4fNdl4;
  
  //Add iron limitation:
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

//N:C ratio at avg. size
   cff1   = 1.-Qmin/Qmax
   QN     = Qmin/(1.-cff1*fN)
  cff     = 1.-cff1*fN
  dcffdl  = -cff1*dfNdl
  dQNdL   = cff1/cff**2 * dfNdl * Qmin
  d2QNdL2 = Qmin*cff1*(d2fNdl2/cff**2 - 2./cff**3*dfNdl*dcffdl)  #Correct

cff     = (thetamax - thetamin)/PAR_

theta   = thetamin+muNet/aI*cff   #Unit: gChl/molC

dthetadl = cff * dY_Xdl(muNet, aI, dmudl, aI*alphaI)  #Correct

d2thetadl2 = cff*d2Y_Xdl2(muNet,aI,dmudl,aI*alphaI, d2mudl2, aI*alphaI**2) //Correct


//return Rcpp::List::create(Rcpp::Named("vec") = someVector,
//                          Rcpp::Named("lst") = someList,
//                          Rcpp::Named("vec2") = someOtherVector);
