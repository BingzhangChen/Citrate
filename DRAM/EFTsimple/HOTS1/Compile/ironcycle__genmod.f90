        !COMPILER-GENERATED INTERFACE MODULE: Tue Oct 10 17:23:05 2017
        MODULE IRONCYCLE__genmod
          INTERFACE 
            SUBROUTINE IRONCYCLE(TEMP,DET,PP_NZ,PP_ND,PP_PN,PP_DZ,      &
     &FE_SCAV,DETFE,DFE)
              REAL(KIND=8), INTENT(IN) :: TEMP
              REAL(KIND=8), INTENT(IN) :: DET
              REAL(KIND=8), INTENT(IN) :: PP_NZ
              REAL(KIND=8), INTENT(IN) :: PP_ND
              REAL(KIND=8), INTENT(IN) :: PP_PN
              REAL(KIND=8), INTENT(IN) :: PP_DZ
              REAL(KIND=8), INTENT(OUT) :: FE_SCAV
              REAL(KIND=8), INTENT(INOUT) :: DETFE
              REAL(KIND=8), INTENT(INOUT) :: DFE
            END SUBROUTINE IRONCYCLE
          END INTERFACE 
        END MODULE IRONCYCLE__genmod
