C**********************************************************************
C     WAVE EQUATION    EXPLICIT METHOD                                *
C**********************************************************************
      PARAMETER(NX=51)
      DIMENSION U(NX),UU(NX),V(NX),VV(NX)
C
C**** INPUT & CALCULATE PARAMETERS
C
      WRITE(*,*) 'Input number of mesh (<51)  (20)'
       READ(*,*) KX
         MX = KX + 1
      WRITE(*,*) 'Input number of time step  (1000)'
       READ(*,*) KM
   99 WRITE(*,*) 'Input time increment dT  (.001)'
       READ(*,*) DT
        DX = 1./FLOAT(MX-1)
         R = DT/DX
        IF(R.GT.1) GO TO 99
        I5 = .125/DT
        IF(I5.EQ.0) I5=1
        PAI= ATAN(1.)*4.
        CALL OUTPUT(U,NX,MX,1)
C
C***  Initial Condition
      DO 10 I = 1,MX
            X = FLOAT(I-1)/FLOAT(MX-1)
         U(I) = 0.25*SIN(PAI*X)
         V(I) = 0.25*SIN(PAI*X)
   10 CONTINUE
C
C***  Time Marching
      DO 20 K = 1,KM
         U(1) = 0.
        U(MX) = 0.
C
        DO 60 I = 1,MX
          VV(I) = U(I)+.3
   60   CONTINUE
C
        IF(MOD(K,I5).EQ.1) CALL OUTPUT(VV,NX,MX,2)
C
        DO 30 I = 2,MX-1
          UU(I) = R*R*(U(I-1)-2.*U(I)+U(I+1))+2.*U(I)-V(I)
   30   CONTINUE
C
        DO 50 I = 2,MX-1
           V(I) = U(I)
           U(I) = UU(I)
   50   CONTINUE
C
   20 CONTINUE
C
        CALL OUTPUT(VV,NX,MX,3)
      STOP
      END
C
C**********************************************************************
C     SUBROUTINE FOR OUTPUT                                           *
C**********************************************************************
      SUBROUTINE OUTPUT(U,NX,MX,MM)
      CHARACTER Z(60,23)
      DIMENSION U(NX)
C
      IF(MM.EQ.1) THEN
C
        DO 10 J  = 1,23
        DO 10 I  = 2,MX-1
          Z(I,J) =' '
   10   CONTINUE
        DO 20 J  = 1,23
          Z(1,J) = ':'
          Z(MX,J)= ':'
   20   CONTINUE
        DO 30 I  = 1,MX
          Z(I,1) = '-'
          Z(I,23)= '-'
   30   CONTINUE
C
      ELSE IF(MM.EQ.2) THEN
C
        DO 40 I = 1,MX
              J = U(I)*40+.001
              IF(J.GT.23.OR.J.LE.0) GO TO 40
              Z(I,J) = '*' 
   40   CONTINUE
C
      ELSE
C
        DO 50 J = 23,1,-1
          WRITE(*,600) (Z(I,J),I=1,60)
   50   CONTINUE
  600   FORMAT(1H ,60A1)
C
      END IF
C
      RETURN
      END
