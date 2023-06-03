C**********************************************************************
C     FLOW PAST MANY OBSTACLES                                        *
C**********************************************************************
      PARAMETER (ID=61,JD=31)
      COMMON U(ID,JD),V(ID,JD),P(ID,JD)
     1      ,Q(ID,JD),D(ID,JD),S(ID,JD)
     1      ,IFL(ID,JD),IBD(ID,JD),IWT(ID,JD)
     1      ,UL(ID),X(ID),XA(ID),XB(ID),XC(ID),XD(ID),XG(ID)
     1      ,Y(JD),YA(JD),YB(JD),YC(JD),YD(JD),YG(JD)
      DIMENSION L1(4),L2(4),L3(4)
      DATA L1/11,15,7,16/,L2/25,29,12,24/,L3/38,43,9,17/
      WRITE(*,*) 'INPUT TIME INCREMENT DT  (0.005) '
       READ(*,*) DT
      WRITE(*,*) 'INPUT REYNOLDS NUMBER (40) '
       READ(*,*) RE
      WRITE(*,*) 'INPUT NUMBER OF TIME STEP (400) '
       READ(*,*) LM
      WRITE(*,*) 'Do you need a test calculation? (YES=1,2 NO=0) (1)'
      WRITE(*,*) ' 1... Many obstacles  2... Circular cylinder '
       READ(*,*) ITEST
C
       IF(ITEST.EQ.1) THEN
C***   X AND Y
         I21=ID
         J21=JD
         DO 210 I = 1,I21-10
             X(I) = 4.*FLOAT(I-1)/FLOAT(I21-1)
  210    CONTINUE
         DO 215 I = I21-9,I21
             X(I) = X(I-1) + (X(I21-10)-X(I21-11))*1.2**(I-I21+10)
  215    CONTINUE
         DO 220 J = 1,J21
             Y(J) = FLOAT(J-1)/FLOAT(J21-1)
  220    CONTINUE
C***   IFL
         DO 235 J = 1,J21
           DO 235 I = 1,I21
             IFL(I,J) = 1
  235    CONTINUE
         DO 240 J = L1(3),L1(4)
           DO 240 I = L1(1),L1(2)
             IFL(I,J) = 0
  240    CONTINUE
         DO 255 J = L2(3),L2(4)
           DO 255 I = L2(1),L2(2)
             IFL(I,J) = 0
  255    CONTINUE
         DO 260 J = L3(3),L3(4)
           DO 260 I = L3(1),L3(2)
             IFL(I,J) = 0
  260    CONTINUE
C
       ELSE IF (ITEST.EQ.2) THEN
C
       I21=ID
       J21=JD
         I11  =  (I21+1)/3
         J11  =  (J21+1)/2
       DO 301 I=1,ID
       X(I)=FLOAT(I-1)/FLOAT(I21-1)*4.
  301  CONTINUE
       DO 302 J=1,JD
       Y(J)=FLOAT(J-1)/FLOAT(J21-1)
  302  CONTINUE
       DO 303 J=1,J21
       DO 303 I=1,I21
          STL=SQRT((X(I)-1.3)**2+(Y(J)-0.5)**2)
          IFL(I,J)=1
          IF(STL.LT..2) IFL(I,J)=0
  303  CONTINUE
C
       ELSE
C
         READ(35,*) NX,NY
         READ(35,*) ((IFL(I,J),I=1,NX),J=1,NY)
         READ(35,*) (X(I),I=1,NX)
         READ(35,*) (Y(J),J=1,NY)
         I21  =  NX
         J21  =  NY
C
       END IF
C
         TD   =  1./DT
         R1   =  1./RE
         I20  =  I21-1
         I19  =  I21-2
         J20  =  J21-1
         J19  =  J21-2
         KM   =  10
C
           DO 222 J = 1,J21
           DO 222 I = 1,I21
               IBD(I,J) = 0
               IWT(I,J) = 0
  222      CONTINUE
C
C***  CALCULATION OF METRICS (ONE-DIMENSION)
           DO 11 I = 2,I20
               XG(I) = 2./(X(I+1)-X(I-1))
               XA(I) = XG(I)*XG(I)
               XB(I) = (X(I+1)-2.*X(I)+X(I-1))*XG(I)**3
               XC(I) = XA(I)-XB(I)*.5
               XD(I) = XA(I)+XB(I)*.5
   11      CONTINUE
           DO 91 J = 2,J20
               YG(J) = 2./(Y(J+1)-Y(J-1))
               YA(J) = YG(J)*YG(J)
               YB(J) = (Y(J+1)-2.*Y(J)+Y(J-1))*YG(J)**3
               YC(J) = YA(J)-YB(J)*.5
               YD(J) = YA(J)+YB(J)*.5
   91      CONTINUE
C
C***  CALCULATION OF THE POSITION OF THE OBSTACLE AND WEIGHT
           DO 230 J = 1,J21
              DO 231 I = 1,I21-1
                IF(IFL(I+1,J)-IFL(I,J).EQ.-1) IBD(I,J) = 1
  231         CONTINUE
              DO 232 I = I21,2,-1
                IF(IFL(I,J)-IFL(I-1,J).EQ.1) IBD(I,J) = 1
  232         CONTINUE
  230      CONTINUE
           DO 250 I = 1,I21
              DO 251 J = 1,J21-1
                IF(IFL(I,J+1)-IFL(I,J).EQ.-1) IBD(I,J) = 1
  251         CONTINUE
              DO 252 J = J21,2,-1
                IF(IFL(I,J)-IFL(I,J-1).EQ.1) IBD(I,J)=1
  252         CONTINUE
  250      CONTINUE
C
           DO 291 J = 1,J21
           DO 291 I = 1,I21
              IBD(I,J) = IBD(I,J)*IFL(I,J)
  291      CONTINUE
c           DO 293 J = 1,J21
c             WRITE(*,601) (IFL(I,J),I=1,I21),(IBD(I,J),I=1,I21)
c  293      CONTINUE
c  600      FORMAT(1H ,'  ',41I2)
c  601      FORMAT(1H ,'  ',41I1,'  ',41I1)
C
           DO 270 J = 2,J21-1
           DO 270 I = 2,I21-1
              IQQ = IBD(I+1,J)+IBD(I-1,J)+IBD(I,J+1)+IBD(I,J-1)
              IQQ = IQQ*(1-IFL(I,J))
                  IF(IQQ.EQ.1) IWT(I,J) = 12
                  IF(IQQ.EQ.2) IWT(I,J) = 6
                  IF(IQQ.EQ.3) IWT(I,J) = 4
                  IF(IQQ.EQ.4) IWT(I,J) = 3
  270      CONTINUE
C
C***  INITIAL CONDITION
           DO 20 J = 1,J21
           DO 20 I = 1,I21
               U(I,J) = 1.
               V(I,J) = 0.
               P(I,J) = 0.
   20      CONTINUE
C
C***  MAIN LOOP
      DO 30 L = 1,LM
C
C***  BOUNDARY CONDITION FOR VELOCITY
         DO 40 I = 1,I21
            U(I,J21) = U(I,J19)
            V(I,J21) =-V(I,J19)
            U(I,1)   = U(I,3)
            V(I,1)   =-V(I,3)
            U(I,J20) = U(I,J19)
            V(I,J20) = 0.
            U(I,2)   = U(I,3)
            V(I,2)   = 0.
   40    CONTINUE
         DO 41 J = 1,J21
            U(I21,J) = U(I19,J)
            V(I21,J) = V(I19,J)
            U(I20,J) = U(I19,J)
            V(I20,J) = V(I19,J)
C           U(1,J)  = U(3,J)*0.+1.
C           IF(J.LT.J11) U(1,J)=0.
C           V(1,J)  =-V(3,J)*0.
   41    CONTINUE
C
C***  CALCULATION OF RIGHT HAND SIDE OF POISSON EQUATION
         DO 22 J = 2,J20
         DO 22 I = 2,I20
            U1     = (U(I+1,J)-U(I-1,J))*.5*XG(I)
            U2     = (U(I,J+1)-U(I,J-1))*.5*YG(J)
            V1     = (V(I+1,J)-V(I-1,J))*.5*XG(I)
            V2     = (V(I,J+1)-V(I,J-1))*.5*YG(J)
            Q(I,J) =-U1*U1-V2*V2-2.*U2*V1+TD*(U1+V2)
   22 CONTINUE
         DO 15 J = 2,J20
         DO 15 I = 2,I20
            D(I,J) = .5/(XA(I)+YA(J))
   15    CONTINUE
C
C***  SOLVE POISSON EQ. FOR PRESSURE BY SOR METHOD
        DO 23 KK = 1,KM
           G1 = 0.
C
           DO 285 J = 2,J21-1
           DO 285 I = 2,I21-1
             S(I,J) = IBD(I+1,J)*P(I+1,J)+IBD(I-1,J)*P(I-1,J)
     1               +IBD(I,J+1)*P(I,J+1)+IBD(I,J-1)*P(I,J-1)
  285 CONTINUE
           DO 286 J = 2,J21-1
           DO 286 I = 2,I21-1
             P(I,J) = S(I,J)*(1-IFL(I,J))*IWT(I,J)/12.+P(I,J)*IFL(I,J)
  286      CONTINUE
C
C***  BOUNDARY CONDITION FOR PRESSURE
           DO 24 J = 1,J21
              P(2,J)   = P(3,J)
              P(I20,J) = P(I19,J)
   24      CONTINUE
           DO 25 I = 1,I21
              P(I,2)   = P(I,3)
              P(I,J20) = P(I,J19)
   25      CONTINUE
C
           DO 26  J = 3,J19
             DO 216 I = 3,I19
                UL(I) = D(I,J)*(XC(I)*P(I+1,J)+YC(J)*P(I,J+1)
     1                 +XD(I)*P(I-1,J)+YD(J)*P(I,J-1)-Q(I,J)) -P(I,J)
  216        CONTINUE
             DO 217 I = 3,I19
  217              G1 = G1+UL(I)*UL(I)*IFL(I,J)
             DO 218 I = 3,I19
  218          P(I,J) = UL(I)+P(I,J)*IFL(I,J)
   26      CONTINUE
           IF( G1.LE..0001) GO TO 31
   23   CONTINUE
   31   WRITE(*,*) L,KK,G1
C
C***  TIME INTEGRATION OF NS EQ.
        DO 27 J = 3,J19
        DO 27 I = 3,I19
            UNX = U(I,J)*(U(I+1,J)-U(I-1,J))/2.
     1           -ABS(U(I,J))*(U(I+1,J)-2.*U(I,J)+U(I-1,J))/2.
            UNY = V(I,J)*(U(I,J+1)-U(I,J-1))/2.
     1           -ABS(V(I,J))*(U(I,J+1)-2.*U(I,J)+U(I,J-1))/2.
C
            VNX = U(I,J)*(V(I+1,J)-V(I-1,J))/2.
     1           -ABS(U(I,J))*(V(I+1,J)-2.*V(I,J)+V(I-1,J))/2.
            VNY = V(I,J)*(V(I,J+1)-V(I,J-1))/2.
     1           -ABS(V(I,J))*(V(I,J+1)-2.*V(I,J)+V(I,J-1))/2.
C
            UV  = (U(I+1,J)-2.*U(I,J)+U(I-1,J))*XA(I)
     1           +(U(I,J+1)-2.*U(I,J)+U(I,J-1))*YA(J)
     2           -(U(I+1,J)-U(I-1,J))*.5*XB(I)
     3           -(U(I,J+1)-U(I,J-1))*.5*YB(J)
            VV  = (V(I+1,J)-2.*V(I,J)+V(I-1,J))*XA(I)
     1           +(V(I,J+1)-2.*V(I,J)+V(I,J-1))*YA(J)
     2           -(V(I+1,J)-V(I-1,J))*.5*XB(I)
     3           -(V(I,J+1)-V(I,J-1))*.5*YB(J)
C
         D(I,J) = U(I,J)+DT*(-(UNX*XG(I)+UNY*YG(J))
     1           -(P(I+1,J)-P(I-1,J))*.5*XG(I)+R1*UV)
         Q(I,J) = V(I,J)+DT*(-(VNX*XG(I)+VNY*YG(J))
     1           -(P(I,J+1)-P(I,J-1))*.5*YG(J)+R1*VV)
   27   CONTINUE
        DO 245 J = 3,J19
        DO 245 I = 3,I19
          U(I,J) = D(I,J)*IFL(I,J)
          V(I,J) = Q(I,J)*IFL(I,J)
  245   CONTINUE
   30 CONTINUE
C
C***  WRITE DATA
       WRITE(*,*) 'WRITE TYPICAL DATA (YES=1)  '
        READ(*,*) IWRITE
          IF(IWRITE.EQ.1) THEN
              WRITE(*,*) '---U---'
              WRITE(*,*) ((U(I,J),I=2,I21,4),J=2,J21,4)
              WRITE(*,*) '---V---'
              WRITE(*,*) ((V(I,J),I=2,I21,4),J=2,J21,4)
              WRITE(*,*) '---P---'
              WRITE(*,*) ((P(I,J),I=2,I21,4),J=2,J21,4)
          END IF
       WRITE(*,*) 'Save Data? (Yes=1,No=0) '
        READ(*,*) ISAVE
        IF(ISAVE.EQ.1) THEN
          WRITE(7,*) I21,J21
           DO 120 J = 1,J21
           DO 120 I = 1,I21
              WRITE(7,*) U(I,J),V(I,J),P(I,J)
  120      CONTINUE
        END IF
C
        CALL OUT3(I21,J21)
      STOP
      END
C**********************************************************************
C     SUBPOUTINE FOR OUTPUT                                           *
C**********************************************************************
      SUBROUTINE OUT3(I21,J21)
      PARAMETER (ID=61,JD=31)
      COMMON U(ID,JD),V(ID,JD),P(ID,JD)
     1      ,Q(ID,JD),D(ID,JD),S(ID,JD)
     1      ,IFL(ID,JD),IBD(ID,JD),IWT(ID,JD)
     1      ,UL(ID),X(ID),XA(ID),XB(ID),XC(ID),XD(ID),XG(ID)
     1      ,Y(JD),YA(JD),YB(JD),YC(JD),YD(JD),YG(JD)
      CHARACTER*1 Z(ID,JD)
      DIMENSION INDEX(ID,JD)
C
        DO 310 I = 1,I21
          P(I,1) = 0.
            DO 310 J = 2,J21
              P(I,J) = P(I,J-1) + .5*(U(I,J)+U(I,J-1))*(Y(J)-Y(J-1))
  310   CONTINUE
      PMIN=P(1,1)
      DO 10 J=1,J21
      DO 10 I=1,I21
        IF(P(I,J).LT.PMIN) PMIN=P(I,J)
   10 CONTINUE
      PMAX=P(1,1)
      DO 30 J=1,J21
      DO 30 I=1,I21
        IF(P(I,J).GT.PMAX) PMAX=P(I,J)
   30 CONTINUE
      WRITE(*,*) 'PMAX=',PMAX,'   PMIN=',PMIN
      DO 40 J2=1,J21
      DO 40 I2=1,I21
        INDEX(I2,J2)=0
        XQ=FLOAT(I2-1)/FLOAT(I21-1)*4.
        YQ=FLOAT(J2-1)/FLOAT(J21-1)
        DO 60 I=1,I21
        IF(X(I).GE.XQ) THEN
        IQ = I
        GO TO 65
        END IF
   60   CONTINUE
   65   CONTINUE
        DO 70 J=1,J21
        IF(Y(J).GE.YQ) THEN
        JQ = J
        GO TO 75
        END IF
   70   CONTINUE
   75   CONTINUE
            F0=(X(IQ+1)-X(IQ))*(Y(JQ+1)-Y(JQ))
            F1=(XQ-X(IQ))*(YQ-Y(JQ))/F0
            F2=(X(IQ+1)-XQ)*(YQ-Y(JQ))/F0
            F3=(X(IQ+1)-XQ)*(Y(JQ+1)-YQ)/F0
            F4=(XQ-X(IQ))*(Y(JQ+1)-YQ)/F0
        PP=P(IQ,JQ)*F3+P(IQ+1,JQ)*F4+P(IQ+1,JQ+1)*F1+P(IQ,JQ+1)*F2
        IND=INT((PP-PMIN)/(PMAX-PMIN)*15.)
        INDEX(I2,J2)=MOD(IND,9)+1
   40   CONTINUE
        DO 80 J=1,J21
        DO 80 I=1,I21
        INDEX(I,J)=INDEX(I,J)*IFL(I,J)
   80   CONTINUE
      DO 85 J=1,J21
      DO 85 I=1,I21
        IND=INDEX(I,J)
        IF(IND.EQ.0) THEN
          Z(I,J)='//'
        ELSE
          IF(MOD(IND,2).EQ.0) Z(I,J)='..'
          IF(MOD(IND,2).EQ.1) Z(I,J)='**'
        END IF
   85  CONTINUE
       DO 90 J=J21,1,-1
        WRITE(*,600) (INDEX(I,J),I=1,I21)
   90  CONTINUE
  600   FORMAT(1H ,80I1)
      WRITE(*,*) 'INPUT ANY NUMBER '
       READ(*,*) AA
      WRITE(*,*) '  '
      DO 95 J=J21,1,-1
        WRITE(*,601) (Z(I,J),I=1,I21)
   95 CONTINUE
  601   FORMAT(1H ,80A1)
      RETURN
      END
