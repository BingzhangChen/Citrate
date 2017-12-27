#include <Rcpp.h>
#include <math.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP TEMPBOL(double Ea, SEXP tC_){
  //DESCRIPTION:
  //The temperature dependence of plankton rates are fomulated according to the Arrhenuis equation. 
  // tC: in situ temperature
  // Tr: reference temperature
  //
  //INPUT PARAMETERS:
  // boltzman constant constant [ eV /K ]
  double kb = 8.62E-5;
  double Tr = 15.;
  NumericVector tC(clone(tC_));

  int N = tC.size();
  NumericVector Y(N);

  //for (int i = 0; i < N; i++) {
   // Y[i] = exp(-(Ea/kb)*(1./(273.15 + tC[i])-1./(273.15 + Tr)));
    Y = exp(-(Ea/kb)*(1./(273.15 + tC)-1./(273.15 + Tr)));
 // }
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

//return Rcpp::List::create(Rcpp::Named("vec") = someVector,
//                          Rcpp::Named("lst") = someList,
//                          Rcpp::Named("vec2") = someOtherVector);
