C**********************************************************************
C     DIFFUSION EQUATION   EULER EXPLICIT METHOD                      *
C**********************************************************************
      PARAMETER(NX=51)
      DIMENSION U(NX),UU(NX)
C
C**** INPUT & CALCULATE PARAMETERS
      WRITE(*,*) 'Input number of mesh. (20)'
       READ(*,*) KX
       MX = KX + 1
      WRITE(*,*) 'Input number of time steps. (250)'
       READ(*,*) KM
      WRITE(*,*) 'Input time increment DELTA T. (0.001)'
       READ(*,*) DT
C
         DX = 1./FLOAT(MX-1)
          R = DT/DX**2
         IH = (MX+1)/2
         I5 = .05/DT
            IF(I5.EQ.0) I5=1
         CALL OUTPUT(U,NX,MX,1)
C
C**** INITIAL CONDITION
        DO 10 I = 1,MX
              X = FLOAT(I-1)/FLOAT(MX-1)
         IF(I.LE.IH) THEN
           U(I) = X
         ELSE
           U(I) = 1.-X
         END IF
   10   CONTINUE
C
C**** MAIN LOOP
      DO 20 K = 1,KM
         U(1) = 0.
         U(MX)= 0.
         IF(MOD(K,I5).EQ.1) CALL OUTPUT(U,NX,MX,2)
C
         DO 30 I = 2,MX-1
           UU(I) = R*U(I-1)+(1.-2*R)*U(I)+R*U(I+1)
   30    CONTINUE
         DO 40 I = 2,MX-1
            U(I) = UU(I)
   40    CONTINUE
C
         IF(ABS(U(IH)).GE.10000.) THEN
             WRITE(*,*) 'DIVERGE!'
             STOP
         END IF
C
   20 CONTINUE
C
      CALL OUTPUT(U,NX,MX,3)
      STOP
      END
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