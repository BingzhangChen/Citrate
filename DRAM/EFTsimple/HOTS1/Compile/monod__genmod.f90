        !COMPILER-GENERATED INTERFACE MODULE: Tue Oct 10 17:22:59 2017
        MODULE MONOD__genmod
          INTERFACE 
            SUBROUTINE MONOD(TEMP,PAR,NO3,PO4,MU0,QNMIN,QPMIN,AI0,BI0,KN&
     &,KP,DFE,KFE,MUNET,QN,QP,THETA,SI,LNO3)
              REAL(KIND=8), INTENT(IN) :: TEMP
              REAL(KIND=8), INTENT(IN) :: PAR
              REAL(KIND=8), INTENT(IN) :: NO3
              REAL(KIND=8), INTENT(IN) :: PO4
              REAL(KIND=8), INTENT(IN) :: MU0
              REAL(KIND=8), INTENT(IN) :: QNMIN
              REAL(KIND=8), INTENT(IN) :: QPMIN
              REAL(KIND=8), INTENT(IN) :: AI0
              REAL(KIND=8), INTENT(IN) :: BI0
              REAL(KIND=8), INTENT(IN) :: KN
              REAL(KIND=8), INTENT(IN) :: KP
              REAL(KIND=8), INTENT(IN) :: DFE
              REAL(KIND=8), INTENT(IN) :: KFE
              REAL(KIND=8), INTENT(OUT) :: MUNET
              REAL(KIND=8), INTENT(OUT) :: QN
              REAL(KIND=8), INTENT(OUT) :: QP
              REAL(KIND=8), INTENT(OUT) :: THETA
              REAL(KIND=8), INTENT(OUT) :: SI
              REAL(KIND=8), INTENT(OUT) :: LNO3
            END SUBROUTINE MONOD
          END INTERFACE 
        END MODULE MONOD__genmod
