C**********************************************************************
C     DIFFUSION EQUATION   EULER IMPLICIT METHOD                      *
C**********************************************************************
      PARAMETER(NX=51)
      DIMENSION U(NX),UU(NX),A(NX),B(NX),C(NX),D(NX)
C
C**** INPUT & CALCULATE PARAMETERS
      WRITE(*,*) 'Input number of mesh (<51)... (40) '
       READ(*,*) KX
       MX = KX + 1
      WRITE(*,*) 'Input number of time steps (50)'
       READ(*,*) KM
      WRITE(*,*) 'Input time increment DELTA T. (0.01)'
       READ(*,*) DT
C
          DX = 1./FLOAT(MX-1)
           R = DT/DX**2
          IH = (MX+1)/2
          I5 =.05/DT
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
   10 CONTINUE
C
C****  IMPLICIT METHOD
      DO 20 K = 1,KM
        U(1)  = 0.
        U(MX) = 0.
        IF(MOD(K,I5).EQ.1)  CALL OUTPUT(U,NX,MX,2)
C
        DO 90 I = 1,MX
           A(I) = R
           B(I) = -2.*R-1.
           C(I) = R
           D(I) = -U(I)
   90   CONTINUE
C
        CALL THOMAS(1,MX,A,B,C,D)
C
        DO 50 I = 2,MX-1
           U(I) = D(I)
   50   CONTINUE
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
C
C**********************************************************************
C     SUBROUTINE FOR SOLVING TRI-DIAGONAL MATRIX                      *
C**********************************************************************
      SUBROUTINE THOMAS(IL,IU,A,B,C,D)
      DIMENSION A(1),B(1),C(1),D(1)
C
      IP=IL+1
      DO 10 I = IP,IU
            R = C(I)/B(I-1)
         B(I) = B(I)-R*A(I-1)
         D(I) = D(I)-R*D(I-1)
   10 CONTINUE
C
      D(IU) = D(IU)/B(IU)
      DO 20 I = IP,IU
            J = IU-I+IL
         D(J) = (D(J)-A(J)*D(J+1))/B(J)
   20 CONTINUE
C
      RETURN
      END
