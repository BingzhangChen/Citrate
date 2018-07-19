        !COMPILER-GENERATED INTERFACE MODULE: Tue May  2 15:02:01 2017
        MODULE MONOD__genmod
          INTERFACE 
            SUBROUTINE MONOD(TEMP,PAR,NO3,MU0,QMIN,AI0,BI0,KN,DFE,KFE,  &
     &MUNET,QN,THETA,SI,LNO3)
              REAL(KIND=8), INTENT(IN) :: TEMP
              REAL(KIND=8), INTENT(IN) :: PAR
              REAL(KIND=8), INTENT(IN) :: NO3
              REAL(KIND=8), INTENT(IN) :: MU0
              REAL(KIND=8), INTENT(IN) :: QMIN
              REAL(KIND=8), INTENT(IN) :: AI0
              REAL(KIND=8), INTENT(IN) :: BI0
              REAL(KIND=8), INTENT(IN) :: KN
              REAL(KIND=8), INTENT(IN) :: DFE
              REAL(KIND=8), INTENT(IN) :: KFE
              REAL(KIND=8), INTENT(OUT) :: MUNET
              REAL(KIND=8), INTENT(OUT) :: QN
              REAL(KIND=8), INTENT(OUT) :: THETA
              REAL(KIND=8), INTENT(OUT) :: SI
              REAL(KIND=8), INTENT(OUT) :: LNO3
            END SUBROUTINE MONOD
          END INTERFACE 
        END MODULE MONOD__genmod
