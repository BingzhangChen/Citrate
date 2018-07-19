        !COMPILER-GENERATED INTERFACE MODULE: Thu Apr 13 17:24:30 2017
        MODULE DIFF_CENTER__genmod
          INTERFACE 
            SUBROUTINE DIFF_CENTER(N,DT,CNPAR,POSCONC,H,NUY,YIN,YOUT)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: DT
              REAL(KIND=8), INTENT(IN) :: CNPAR
              INTEGER(KIND=4), INTENT(IN) :: POSCONC
              REAL(KIND=8), INTENT(IN) :: H(N)
              REAL(KIND=8), INTENT(IN) :: NUY(0:N)
              REAL(KIND=8), INTENT(IN) :: YIN(N)
              REAL(KIND=8), INTENT(OUT) :: YOUT(N)
            END SUBROUTINE DIFF_CENTER
          END INTERFACE 
        END MODULE DIFF_CENTER__genmod
