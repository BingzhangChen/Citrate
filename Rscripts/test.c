#include <stdio.h>
#define DT 60.    // time step 60 seconds
#define D2SEC 86400.    // How many seconds in one day
#define NDAYS 36500     // How many years to run

int main() 
{
    // A simple NPZ model
    double N, P, Z;
    double dtday;
    int    Nsteps, i;

    Nsteps = NDAYS * int(D2SEC)/int(DT); // The total number of time steps
    dtday  = DT/D2SEC;  // Time of one step in terms of days
    N = 2;
    P = 0.1;
    Z = 0.01;


}

