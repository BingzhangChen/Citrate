        !COMPILER-GENERATED INTERFACE MODULE: Tue Oct 10 17:23:07 2017
        MODULE TRANSFORM__genmod
          INTERFACE 
            SUBROUTINE TRANSFORM(N,X,XMIN,XMAX,CONVERT,Y)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: X(N)
              REAL(KIND=8), INTENT(IN) :: XMIN
              REAL(KIND=8), INTENT(IN) :: XMAX
              INTEGER(KIND=4), INTENT(IN) :: CONVERT
              REAL(KIND=8), INTENT(OUT) :: Y(N)
            END SUBROUTINE TRANSFORM
          END INTERFACE 
        END MODULE TRANSFORM__genmod
