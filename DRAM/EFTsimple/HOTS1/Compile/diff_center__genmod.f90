        !COMPILER-GENERATED INTERFACE MODULE: Tue Oct 10 17:23:00 2017
        MODULE DIFF_CENTER__genmod
          INTERFACE 
            SUBROUTINE DIFF_CENTER(N,DT,CNPAR,POSCONC,H,BCUP,BCDW,YUP,  &
     &YDW,NUY,LSOUR,QSOUR,TAUR,YOBS,YIN,YOUT)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: DT
              REAL(KIND=8), INTENT(IN) :: CNPAR
              INTEGER(KIND=4), INTENT(IN) :: POSCONC
              REAL(KIND=8), INTENT(IN) :: H(N)
              INTEGER(KIND=4), INTENT(IN) :: BCUP
              INTEGER(KIND=4), INTENT(IN) :: BCDW
              REAL(KIND=8), INTENT(IN) :: YUP
              REAL(KIND=8), INTENT(IN) :: YDW
              REAL(KIND=8), INTENT(IN) :: NUY(0:N)
              REAL(KIND=8), INTENT(IN) :: LSOUR(N)
              REAL(KIND=8), INTENT(IN) :: QSOUR(N)
              REAL(KIND=8), INTENT(IN) :: TAUR(N)
              REAL(KIND=8), INTENT(IN) :: YOBS(N)
              REAL(KIND=8), INTENT(IN) :: YIN(N)
              REAL(KIND=8), INTENT(OUT) :: YOUT(N)
            END SUBROUTINE DIFF_CENTER
          END INTERFACE 
        END MODULE DIFF_CENTER__genmod
