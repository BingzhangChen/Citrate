!========================================================
! Simplify this function to speed up the calculation 2016/9/5
!FUNCTION WAPR (X, NB, L) RESULT (WAP)
! NB = 0, L = 0
FUNCTION WAPR (X) RESULT (WAP)
!
!     WAPR - output
!     X - argument of W(X)
!     NB is the branch of the W function needed:
!        NB = 0 - upper branch
!        NB <> 0 - lower branch
!
!     NERROR is the output error flag:
!        NERROR = 0 -> routine completed successfully
!        NERROR = 1 -> X is out of range
!
!     Range: -exp(-1) <= X for the upper branch of the W function
!            -exp(-1) < X < 0 for the lower branch of the W function
!
!     L - determines how WAPR is to treat the argument X
!        L = 1 -> X is the offset from -exp(-1), so compute
!                 W(X-exp(-1))
!        L <> 1 -> X is the desired X, so compute W(X)
!
!     M - print messages from WAPR?
!         M = 1 -> Yes
!         M <> 1 -> No
!
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     NN is the output device number
!
!     NBITS is the number of bits (less 1) in the mantissa of the
!        floating point number number representation of your machine.
!        It is used to determine the level of accuracy to which the W
!        function should be calculated.
!
!        Most machines use a 24-bit matissa for single precision and
!        53-56 bits for double precision. The IEEE standard is 53
!        bits. The Fujitsu VP2200 uses 56 bits. Long word length
!        machines vary, e.g., the Cray X/MP has a 48-bit mantissa for
!        single precision.
!
IMPLICIT NONE
real, INTENT(in)   :: X 
INTEGER, PARAMETER :: NN=6, NBITS=23, NITER=1
real,    PARAMETER ::EM=-0.367879441171442,&           ! -EXP(-1)
                     EM9=-1.234098040866796E-4,&       ! -EXP(-9)
                     C13=1.E0/3.E0,&
                     C23=2.E0*C13,&
                     EM2=2.E0/EM,&
                     E12=-EM2,&
                     TB=.5E0**NBITS,&
                     TB2=.5E0**(NBITS/2),&       ! SQRT(TB)
                     X0=0.0350769390096679055,&  ! TB**(1/6)*0.5E0
                     X1=0.302120119432784731,&   !(1 - 17*TB**(2/7))*EM
                     AN22=3.6E2/83.E0,&
                     AN11=135./83.E0,&
                     AN3=8.E0/3.E0,&
                     AN4=135.E0/83.E0,&
                     AN5=166.E0/39.E0,&
                     AN6=3167.E0/3549.E0,&
                     S2=1.41421356237310,& ! SQRT(2.E0)
                     S21=2.E0*S2-3.E0,&
                     S22=4.E0-3.E0*S2,&
                     S23=S2-2.E0
real ::  WAP, AN2, DELX,  RETA, ZL, T, TS, ETA, TEMP, TEMP2, ZN
real ::  XX
!INTEGER, INTENT(in) :: NB, L

!        Various mathematical constants
    
!
!     The following COMMON statement is needed only when testing this
!     function using BISECT, otherwise it can be removed.
!
!    COMMON/WAPCOM/NBITS
!    DATA INIT,NITER/0,1/
!     DATA NBITS/23/
!
!     IF(INIT.EQ.0) THEN
!        INIT=1
!
!        Code to calculate NBITS for the host machine. NBITS is the
!        mantissa length less one. This value is chosen to allow for
!        rounding error in the final bit. This need only be run once on
!        any particular machine. It can then be included in the above
!        DATA statement.
!
!        DO I=1,2000
!           B=2.E0**(-I)
!           V=1.E0+B
!           IF(V.EQ.1.E0)THEN
!              NBITS=I-1
!              J=-ALOG10(B)
!              IF(M.EQ.1) WRITE(NN,40)NBITS,J
!              EXIT
!           ENDIF
!        END DO
!
!        Remove to here after NBITS has been calculated once
!
!        The case of large NBITS
!
!        IF(NBITS.GE.56) NITER=2
!
!        Various mathematical constants
!
!        EM=-EXP(-1.E0)
!        EM9=-EXP(-9.E0)
!        C13=1.E0/3.E0
!        C23=2.E0*C13
!        EM2=2.E0/EM
!        E12=-EM2
!        TB=.5E0**NBITS
!        TB2=SQRT(TB)
!        X0=TB**(1.E0/6.E0)*.5E0
!        X1=(1.E0-17.E0*TB**(2.E0/7.E0))*EM
!        AN22=3.6E2/83.E0
!        AN11=135./83.E0
!        AN3=8.E0/3.E0
!        AN4=135.E0/83.E0
!        AN5=166.E0/39.E0
!        AN6=3167.E0/3549.E0
!        S2=SQRT(2.E0)
!        S21=2.E0*S2-3.E0
!        S22=4.E0-3.E0*S2
!        S23=S2-2.E0
!     ENDIF

!    IF(L.EQ.1) THEN
!       DELX=X
!       IF(DELX.LT.0.E0) THEN
!          WAP = 1./0.
!          RETURN
!       END IF
!       XX=X+EM
!!        IF(E12*DELX.LT.TB**2.AND.M.EQ.1) WRITE(NN,60)DELX
!    ELSE
       IF(X.LT.EM) THEN
          WAP = -1.
          RETURN
       END IF
       IF(X.EQ.EM) THEN
          WAP=-1.E0
          RETURN
       ENDIF
       XX  =X
       DELX=XX-EM
!        IF(DELX.LT.TB2.AND.M.EQ.1) WRITE(NN,70)XX
!    ENDIF
!    IF(NB.EQ.0) THEN
!
!        Calculations for Wp
!
       IF(ABS(XX).LE.X0) THEN
          WAP=XX/(1.E0+XX/(1.E0+XX/(2.E0+XX/(.6E0+.34E0*XX))))
          RETURN
       ELSE IF(XX.LE.X1) THEN
          RETA=SQRT(E12*DELX)
          WAP=RETA/(1.E0+RETA/(3.E0+RETA/(RETA/(AN4+RETA/(RETA*&
               AN6+AN5))+AN3)))-1.E0
          RETURN
       ELSE IF(XX.LE.2.E1) THEN
          RETA=S2*SQRT(1.E0-XX/EM)
          AN2=4.612634277343749E0*SQRT(SQRT(RETA+&
               1.09556884765625E0))
          WAP=RETA/(1.E0+RETA/(3.E0+(S21*AN2+S22)*RETA/&
               (S23*(AN2+RETA))))-1.E0
       ELSE
          ZL =ALOG(XX)
          WAP=ALOG(XX/ALOG(XX/ZL**EXP(-1.124491989777808E0/&
               (.4225028202459761E0+ZL))))
       ENDIF
!    ELSE
!!
!!        Calculations for Wm
!!
!       IF(XX.GE.0.E0) THEN
!          WAP = -1./0.
!          RETURN
!       END IF
!       IF(XX.LE.X1) THEN
!          RETA=SQRT(E12*DELX)
!          WAP=RETA/(RETA/(3.E0+RETA/(RETA/(AN4+RETA/(RETA*&
!               AN6-AN5))-AN3))-1.E0)-1.E0
!          RETURN
!       ELSE IF(XX.LE.EM9) THEN
!          ZL=ALOG(-XX)
!          T=-1.E0-ZL
!          TS=SQRT(T)
!          WAP=ZL-(2.E0*TS)/(S2+(C13-T/(2.7E2+&
!               TS*127.0471381349219E0))*TS)
!       ELSE
!          ZL=ALOG(-XX)
!          ETA=2.E0-EM2*XX
!          WAP=ALOG(XX/ALOG(-XX/((1.E0-.5043921323068457E0*&
!               (ZL+1.E0))*(SQRT(ETA)+ETA/3.E0)+1.E0)))
!       ENDIF
!    ENDIF
!     DO I=1,NITER
       ZN   =ALOG(XX/WAP)-WAP
       TEMP =1D0+WAP
       TEMP2=TEMP+C23*ZN
       TEMP2=2D0*TEMP*TEMP2
       WAP  =WAP*(1D0+(ZN/TEMP)*(TEMP2-ZN)/(TEMP2-2.E0*ZN))
!     END DO
    RETURN
END FUNCTION WAPR

