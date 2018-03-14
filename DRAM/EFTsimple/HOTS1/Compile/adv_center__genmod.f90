        !COMPILER-GENERATED INTERFACE MODULE: Tue Oct 10 17:23:00 2017
        MODULE ADV_CENTER__genmod
          INTERFACE 
            SUBROUTINE ADV_CENTER(N,DT,H,HO,WW,BCUP,BCDW,YUP,YDW,METHOD,&
     &MODE,Y)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: DT
              REAL(KIND=8), INTENT(IN) :: H(N)
              REAL(KIND=8), INTENT(IN) :: HO(N)
              REAL(KIND=8), INTENT(IN) :: WW(0:N)
              INTEGER(KIND=4), INTENT(IN) :: BCUP
              INTEGER(KIND=4), INTENT(IN) :: BCDW
              REAL(KIND=8), INTENT(IN) :: YUP
              REAL(KIND=8), INTENT(IN) :: YDW
              INTEGER(KIND=4), INTENT(IN) :: METHOD
              INTEGER(KIND=4), INTENT(IN) :: MODE
              REAL(KIND=8), INTENT(INOUT) :: Y(N)
            END SUBROUTINE ADV_CENTER
          END INTERFACE 
        END MODULE ADV_CENTER__genmod
