C**********************************************************************
C      FLOW IN THE CUBIC CAVITY
C**********************************************************************
      PARAMETER (MX=21,MY=21,MZ=21)
      COMMON U(MX,MY,MZ),V(MX,MY,MZ),W(MX,MY,MZ),P(MX,MY,MZ)
     1      ,Q(MX,MY,MZ),UL(MX),R(MX,MY)
C***  READ DATA
        WRITE(*,*) 'NUMBERS OF MESH (X,Y,Z DIRECTIONS)? (20,20,20)'
         READ(*,*) I20,J20,K20
        WRITE(*,*) 'TIME INCREMENT DT & REYNOLDS NUMBER RE? (.01,40)'
         READ(*,*) DT,RE
        WRITE(*,*) 'ATTACK ANGLE TETA (DEGREE)? (0.)'
         READ(*,*) TETA
        WRITE(*,*) 'NUMBER OF TIME STEP LM (100) '
         READ(*,*) LM
C***  CALCULATIN OF PARAMETERS
          I21 = I20 + 1
          J21 = J20 + 1
          K21 = K20 + 1
          I19 = I20 - 1
          J19 = J20 - 1
          K19 = K20 - 1
          TD  = 1./DT
          R1  = 1./RE
          KM  = 10
          TET = ATAN(1.)*4.*TETA/180.
          DX  = 1./FLOAT(I20)
          DY  = 1./FLOAT(J20)
          DZ  = 1./FLOAT(K20)
          XD  = 1./DX
          YD  = 1./DY
          ZD  = 1./DZ
          DIJK=1./(2./(DX*DX)+2./(DY*DY)+2./(DZ*DZ))
C
C***  INITIAL CONDITION
         DO 20  K = 1,K21
         DO 20  J = 1,J21
         DO 20  I = 1,I21
           U(I,J,K) = 0.
           V(I,J,K) = 0.
           W(I,J,K) = 0.
           P(I,J,K) = 0.
   20    CONTINUE
C
C***  MAIN LOOP
      DO 30  L = 1,LM
C
C***  BOUNDARY CONDITION FOR VELOCITY
C***  X-Z
         DO 40 K = 1,K21
         DO 40 I = 1,I21
           U(I,J20,K) = -U(I,J19,K)
C          U(I,J20,K) =  0.
           V(I,J21,K) =  V(I,J19,K)
           W(I,J20,K) = -W(I,J19,K)
           V(I,J20,K) =  0.
           U(I,1,K)   = -U(I,2,K)
           V(I,1,K)   =  V(I,3,K)
           W(I,1,K)   = -W(I,2,K)
           V(I,2,K)   =  0.
   40    CONTINUE
C*** Y-Z
         DO 41 K = 1,K21
         DO 41 J = 1,J21
           U(I21,J,K) =  U(I19,J,K)
           V(I20,J,K) = -V(I19,J,K)
           W(I20,J,K) = -W(I19,J,K)
           U(I20,J,K) =  0.
           U(1,J,K)   =  U(3,J,K)
           V(1,J,K)   = -V(2,J,K)
           W(1,J,K)   = -W(2,J,K)
           U(2,J,K)   =  0.
   41    CONTINUE
C*** X-Y
         DO 43 J = 1,J21
         DO 43 I = 1,I21
           U(I,J,1)   = -U(I,J,2)
           V(I,J,1)   = -V(I,J,2)
           W(I,J,1)   =  W(I,J,3)
           W(I,J,2)   =  0.
C          U(I,J,K20) =  2.0 -U(I,J,K19)
           U(I,J,K20) =  COS(TET)
           U(I,J,K21) =  COS(TET)
           V(I,J,K20) =  SIN(TET)
           V(I,J,K21) =  SIN(TET)
           W(I,J,K21) =  W(I,J,K19)
           W(I,J,K20) =  0.
   43    CONTINUE
C
C***  CALCULATION FOR RHS OF POISSON EQ.
         DO 22 K = 2,K20
         DO 22 J = 2,J20
         DO 22 I = 2,I20
             U1  =  (U(I+1,J,K)-U(I,J,K))*XD
             U2  =  (U(I,J+1,K)-U(I,J,K))*YD
             U3  =  (U(I,J,K+1)-U(I,J,K))*ZD
             V1  =  (V(I+1,J,K)-V(I,J,K))*XD
             V2  =  (V(I,J+1,K)-V(I,J,K))*YD
             V3  =  (V(I,J,K+1)-V(I,J,K))*ZD
             W1  =  (W(I+1,J,K)-W(I,J,K))*XD
             W2  =  (W(I,J+1,K)-W(I,J,K))*YD
             W3  =  (W(I,J,K+1)-W(I,J,K))*ZD
         Q(I,J,K)= -U1*U1-V2*V2-W3*W3-2.*(U2*V1+V3*W2+W1*U3)
     1             +TD*(U1+V2+W3)
   22    CONTINUE
C
C***  SOLUTION OF POISSON EQ. FOR PRESSURE BY G-S METHOD
        DO 23  KK = 1,KM
          G1  = 0.
C
C***  BOUNDARY CONDITION FOR PRESSURE
C***  Y-Z
         DO 24 K = 1,K21
         DO 24 J = 1,J21
           P(1,J,K)   = P(2,J,K)
           P(I20,J,K) = P(I19,J,K)
   24    CONTINUE
C***  X-Z
         DO 25 K = 1,K21
         DO 25 I = 1,I21
           P(I,1,K)   = P(I,2,K)
           P(I,J20,K) = P(I,J19,K)
   25    CONTINUE
C***  X-Y
         DO 29 J = 1,J21
         DO 29 I = 1,I21
           P(I,J,1)   = P(I,J,2)
           P(I,J,K20) = P(I,J,K19)
   29    CONTINUE
C***  GAUSS-SEIDEL METHOD 
         DO 26 K = 2,K19
         DO 26 J = 2,J19
         DO 26 I = 2,I19
             PCOR   = DIJK*((P(I+1,J,K)+P(I-1,J,K))/(DX*DX)
     1                +(P(I,J+1,K)+P(I,J-1,K))/(DY*DY)
     2                +(P(I,J,K+1)+P(I,J,K-1))/(DZ*DZ)-Q(I,J,K)) 
     3                -P(I,J,K)
             G1     = G1 + PCOR*PCOR
          P(I,J,K) = P(I,J,K) + PCOR
   26    CONTINUE
         IF( G1.LE..001) GO TO 31
   23 CONTINUE
   31   IF(L-L/20*20.eq.0) WRITE(*,*) L,K,G1
C
C***  CALCULATION OF NEW VELOCITY (N-S EQ.)
         DO 27 K = 2,K19
         DO 27 J = 2,J19
         DO 27 I = 2,I19
C
            UV = (U(I+1,J,K)-2.*U(I,J,K)+U(I-1,J,K))*XD*XD
     1          +(U(I,J+1,K)-2.*U(I,J,K)+U(I,J-1,K))*YD*YD
     2          +(U(I,J,K+1)-2.*U(I,J,K)+U(I,J,K-1))*ZD*ZD
            VV = (V(I+1,J,K)-2.*V(I,J,K)+V(I-1,J,K))*XD*XD
     1          +(V(I,J+1,K)-2.*V(I,J,K)+V(I,J-1,K))*YD*YD
     2          +(V(I,J,K+1)-2.*V(I,J,K)+V(I,J,K-1))*ZD*ZD
            WV = (W(I+1,J,K)-2.*W(I,J,K)+W(I-1,J,K))*XD*XD
     1          +(W(I,J+1,K)-2.*W(I,J,K)+W(I,J-1,K))*YD*YD
     2          +(W(I,J,K+1)-2.*W(I,J,K)+W(I,J,K-1))*ZD*ZD
C
            UN = U(I,J,K)*(U(I+1,J,K)-U(I-1,J,K))*.5*XD
     1          +(V(I-1,J,K)+V(I,J,K)+V(I-1,J+1,K)+V(I,J+1,K))*.25
     2          *(U(I,J+1,K)-U(I,J-1,K))*.5*YD
     3          +(W(I-1,J,K)+W(I,J,K)+W(I-1,J,K+1)+W(I,J,K+1))*.25
     4          *(U(I,J,K+1)-U(I,J,K-1))*.5*ZD
            VN = V(I,J,K)*(V(I,J+1,K)-V(I,J-1,K))*.5*YD
     1          +(U(I,J-1,K)+U(I,J,K)+U(I+1,J-1,K)+U(I+1,J,K))*.25
     2          *(V(I+1,J,K)-V(I-1,J,K))*.5*XD
     3          +(W(I,J-1,K)+W(I,J,K)+W(I,J-1,K+1)+W(I,J,K+1))*.25
     4          *(V(I,J,K+1)-V(I,J,K-1))*.5*ZD
            WN = W(I,J,K)*(W(I,J,K+1)-W(I,J,K-1))*.5*ZD
     1          +(V(I,J,K-1)+V(I,J,K)+V(I,J+1,K-1)+V(I,J+1,K))*.25
     2          *(W(I,J+1,K)-W(I,J-1,K))*.5*YD
     3          +(U(I,J,K-1)+U(I,J,K)+U(I+1,J,K-1)+U(I+1,J,K))*.25
     4          *(W(I+1,J,K)-W(I-1,J,K))*.5*XD
C
        U(I,J,K) = U(I,J,K)+DT*(-UN-(P(I,J,K)-P(I-1,J,K))*XD+R1*UV)
        V(I,J,K) = V(I,J,K)+DT*(-VN-(P(I,J,K)-P(I,J-1,K))*YD+R1*VV)
        W(I,J,K) = W(I,J,K)+DT*(-WN-(P(I,J,K)-P(I,J,K-1))*ZD+R1*WV)
   27   CONTINUE
   30 CONTINUE
C
C***  CALCULATION OF STREAM FUNCTION USING
        DO 28 K = 1,K21
        DO 28 I = 1,I21
          R(I,K) = U(I,6,K)
   28   CONTINUE
          CALL PSI(R,MX,MZ,I21,K21)
        WRITE(*,*) 'INPUT ANY NUMBER'
         READ(*,*) AA
        DO 281 K = 1,K21
        DO 281 I = 1,I21
          R(I,K) = U(I,3,K)
  281   CONTINUE
          CALL PSI(R,MX,MZ,I21,K21)
        WRITE(*,*) 'INPUT ANY NUMBER'
         READ(*,*) AA
C
C***  WRITE TYPICAL DATA
          WRITE(*,*) '---U---'
          WRITE(*,*) ((U(I,J,5),I=2,I20,4),J=2,J20,4)
          WRITE(*,*) '---V---'
          WRITE(*,*) ((V(I,J,5),I=2,I20,4),J=2,J20,4)
          WRITE(*,*) '---P---'
          WRITE(*,*) ((P(I,J,5),I=2,I20,4),J=2,J20,4)
C
C***   SAVE DATA
        WRITE(*,*) 'DO YOU WANT TO SAVE DATA? (YES=1 NO=0)  '
         READ(*,*) ISAVE
         IF (ISAVE.EQ.1) THEN
           WRITE(*,*) 'File Number ? (3)'
           READ(*,*) IS
           WRITE(IS,*) I21,J21,K21
           DO 777 K=1,K21
           DO 777 J=1,J21
           DO 777 I=1,I21
              WRITE(IS,*) U(I,J,K),V(I,J,K),W(I,J,K),P(I,J,K)
  777      CONTINUE
C
        END IF
      STOP
      END
C**********************************************************************
C     SUBROUTINE FOR CALCULATING STREAM FUNCTION
C**********************************************************************
      SUBROUTINE PSI(U,NX,NY,I21,J21)
      DIMENSION U(NX,NY)
      DIMENSION R(51,51)
C
C***  CALCULATING PSI 
        DO 70 I = 1,I21
         R(I,1) = 0.
   70   CONTINUE
        DY=0.1
        DO 71 J = 2,J21
        DO 71 I = 1,I21
         R(I,J) = R(I,J-1)+U(I,J)*DY
   71   CONTINUE
        DO 72 J = 1,J21
        DO 72 I = 1,I21
         U(I,J) = R(I,J)
   72   CONTINUE
           CALL OUT(U,NX,NY,I21,J21)
C
      RETURN
      END
C
C**********************************************************************
C     Subroutine for output                                           *
C**********************************************************************
      SUBROUTINE OUT(U,NX,NY,MX,MY)
      DIMENSION U(NX,NY)
      DIMENSION INDEX(80)
      CHARACTER*2 Z(80)
C
        UMAX = U(1,1)
        UMIN = U(1,1)
        DO 10 J = 1,MY
          DO 10 I = 1,MX
            IF(UMAX.LT.U(I,J)) UMAX = U(I,J)
            IF(UMIN.GT.U(I,J)) UMIN = U(I,J)
   10   CONTINUE
C
        DO 40 J = MY,1,-1
          DO 50 I = 1,MX
           IND = INT((U(I,J)-UMIN)/(UMAX-UMIN)*7.9999)*11
           IF(MOD(IND,2).EQ.0) Z(I)='**'
           IF(MOD(IND,2).EQ.1) Z(I)='  '
   50     CONTINUE
          Z(1)  = '| '
          Z(MX) = ' |'
           IF(J.EQ.MY.OR.J.EQ.1) THEN
           DO 51 I = 1,MX
            Z(I)='--'
   51      CONTINUE
           END IF
           WRITE(*,601) (Z(I),I=1,MX)
   40   CONTINUE
  601 FORMAT(1H ,80A2)
C
       WRITE(*,*) 'INPUT ANY NUMBER (1)'
        READ(*,*) AA
C
        DO 20 J = MY,1,-1
          DO 30 I = 1,MX
           INDEX(I) = INT((U(I,J)-UMIN)/(UMAX-UMIN)*9.9999)*11
   30     CONTINUE
           WRITE(*,600) (INDEX(I),I=1,MX)
   20   CONTINUE
  600 FORMAT(1H ,80I2)
C
      RETURN
      END
