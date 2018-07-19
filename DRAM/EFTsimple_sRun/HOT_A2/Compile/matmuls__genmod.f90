        !COMPILER-GENERATED INTERFACE MODULE: Thu Apr 13 17:24:35 2017
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
