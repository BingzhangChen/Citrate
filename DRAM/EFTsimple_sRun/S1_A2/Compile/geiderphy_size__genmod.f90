        !COMPILER-GENERATED INTERFACE MODULE: Tue May  2 15:02:03 2017
        MODULE GEIDERPHY_SIZE__genmod
          INTERFACE 
            SUBROUTINE GEIDERPHY_SIZE(PMU,NO3,TC,PAR_,THETA,MUNET,SI,   &
     &LNO3)
              REAL(KIND=8), INTENT(IN) :: PMU
              REAL(KIND=8), INTENT(IN) :: NO3
              REAL(KIND=8), INTENT(IN) :: TC
              REAL(KIND=8), INTENT(IN) :: PAR_
              REAL(KIND=8), INTENT(IN) :: THETA
              REAL(KIND=8), INTENT(OUT) :: MUNET
              REAL(KIND=8), INTENT(OUT) :: SI
              REAL(KIND=8), INTENT(OUT) :: LNO3
            END SUBROUTINE GEIDERPHY_SIZE
          END INTERFACE 
        END MODULE GEIDERPHY_SIZE__genmod
