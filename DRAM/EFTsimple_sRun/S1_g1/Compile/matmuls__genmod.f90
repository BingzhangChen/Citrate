        !COMPILER-GENERATED INTERFACE MODULE: Tue May  2 15:02:05 2017
        MODULE MATMULS__genmod
          INTERFACE 
            SUBROUTINE MATMULS(X,V,P)
              USE BIO_MOD, ONLY :                                       &
     &          NPAR
              REAL(KIND=8), INTENT(IN) :: X(NPAR,NPAR)
              REAL(KIND=8), INTENT(IN) :: V(NPAR)
              REAL(KIND=8), INTENT(OUT) :: P(NPAR)
            END SUBROUTINE MATMULS
          END INTERFACE 
        END MODULE MATMULS__genmod
