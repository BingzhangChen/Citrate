        !COMPILER-GENERATED INTERFACE MODULE: Tue May  2 15:01:59 2017
        MODULE TRIDIAGONAL__genmod
          INTERFACE 
            SUBROUTINE TRIDIAGONAL(N,AU,BU,CU,DU,FI,LT,VALUE)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: AU(N)
              REAL(KIND=8), INTENT(IN) :: BU(N)
              REAL(KIND=8), INTENT(IN) :: CU(N)
              REAL(KIND=8), INTENT(IN) :: DU(N)
              INTEGER(KIND=4), INTENT(IN) :: FI
              INTEGER(KIND=4), INTENT(IN) :: LT
              REAL(KIND=8), INTENT(OUT) :: VALUE(N)
            END SUBROUTINE TRIDIAGONAL
          END INTERFACE 
        END MODULE TRIDIAGONAL__genmod
