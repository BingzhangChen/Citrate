        !COMPILER-GENERATED INTERFACE MODULE: Tue Oct 10 17:23:07 2017
        MODULE UPDATECVM__genmod
          INTERFACE 
            SUBROUTINE UPDATECVM(OLDCVM,OLDMEAN,N,PAR,NEWMEAN,NEWCVM)
              USE INTERFACE_MOD, ONLY :                                 &
     &          NPAR,                                                   &
     &          PRIORCVM
              REAL(KIND=8), INTENT(IN) :: OLDCVM(NPAR*(NPAR+1)/2)
              REAL(KIND=8), INTENT(IN) :: OLDMEAN(NPAR)
              REAL(KIND=8), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: PAR(NPAR)
              REAL(KIND=8), INTENT(OUT) :: NEWMEAN(NPAR)
              REAL(KIND=8), INTENT(OUT) :: NEWCVM(NPAR*(NPAR+1)/2)
            END SUBROUTINE UPDATECVM
          END INTERFACE 
        END MODULE UPDATECVM__genmod
