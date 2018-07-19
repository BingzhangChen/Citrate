        !COMPILER-GENERATED INTERFACE MODULE: Thu Apr 13 17:24:32 2017
        MODULE MM__genmod
          INTERFACE 
            SUBROUTINE MM(N,K0,ALPHAK,L,FN,DFNDL,D2FNDL2,D3FNDL3,D4FNDL4&
     &)
              REAL(KIND=8), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: K0
              REAL(KIND=8), INTENT(IN) :: ALPHAK
              REAL(KIND=8), INTENT(IN) :: L
              REAL(KIND=8), INTENT(OUT) :: FN
              REAL(KIND=8), INTENT(OUT) :: DFNDL
              REAL(KIND=8), INTENT(OUT) :: D2FNDL2
              REAL(KIND=8), INTENT(OUT) :: D3FNDL3
              REAL(KIND=8), INTENT(OUT) :: D4FNDL4
            END SUBROUTINE MM
          END INTERFACE 
        END MODULE MM__genmod
