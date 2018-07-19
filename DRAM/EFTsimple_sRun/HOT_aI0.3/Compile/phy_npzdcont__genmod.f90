        !COMPILER-GENERATED INTERFACE MODULE: Thu Apr 13 17:24:32 2017
        MODULE PHY_NPZDCONT__genmod
          INTERFACE 
            SUBROUTINE PHY_NPZDCONT(NO3,PAR_,TEMP_,FE,PMU,MUNET,DMUDL,  &
     &D2MUDL2,D3MUDL3,D4MUDL4,SI,FN,THETA,QN)
              REAL(KIND=8), INTENT(IN) :: NO3
              REAL(KIND=8), INTENT(IN) :: PAR_
              REAL(KIND=8), INTENT(IN) :: TEMP_
              REAL(KIND=8), INTENT(IN) :: FE
              REAL(KIND=8), INTENT(IN) :: PMU
              REAL(KIND=8), INTENT(OUT) :: MUNET
              REAL(KIND=8), INTENT(OUT) :: DMUDL
              REAL(KIND=8), INTENT(OUT) :: D2MUDL2
              REAL(KIND=8), INTENT(OUT) :: D3MUDL3
              REAL(KIND=8), INTENT(OUT) :: D4MUDL4
              REAL(KIND=8), INTENT(OUT) :: SI
              REAL(KIND=8), INTENT(OUT) :: FN
              REAL(KIND=8), INTENT(OUT) :: THETA
              REAL(KIND=8), INTENT(OUT) :: QN
            END SUBROUTINE PHY_NPZDCONT
          END INTERFACE 
        END MODULE PHY_NPZDCONT__genmod
