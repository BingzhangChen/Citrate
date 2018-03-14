        !COMPILER-GENERATED INTERFACE MODULE: Tue Oct 10 17:23:07 2017
        MODULE COV__genmod
          INTERFACE 
            SUBROUTINE COV(NR,DAT,CVM)
              USE INTERFACE_MOD, ONLY :                                 &
     &          NPAR
              INTEGER(KIND=4), INTENT(IN) :: NR
              REAL(KIND=8), INTENT(IN) :: DAT(NR,NPAR)
              REAL(KIND=8), INTENT(INOUT) :: CVM(NPAR*(NPAR+1)/2)
            END SUBROUTINE COV
          END INTERFACE 
        END MODULE COV__genmod
