C**********************************************************************
C     UNSTEADY FLOW & HEAT TRANFER IN THE ROOM                        *
C           MAC METHOD     NON-UNIFORM MESH SIZE                      *
C**********************************************************************
C
      PARAMETER(NX=51,NY=51,NXP=NX+1,NYP=NY+1)
      COMMON U(NX,NY),V(NX,NY),P(NX,NY),R(NX,NY),D(NX,NY),T(NX,NY)
     1,      X(NXP),Y(NYP),A1(NX),B1(NX),C1(NX),A2(NX),B2(NX),C2(NX)
     2,                    A3(NY),B3(NY),C3(NY),A4(NY),B4(NY),C4(NY)
C
C***  READ PARAMETERS
          PRN = 0.71
  99  WRITE(*,*)'INPUT NUMBER OF MESH FOR X & Y DIRECTION(<51) (20,20)'
       READ(*,*) NA,NB
          I21 = NA + 1
          J21 = NB + 1
          I20 = I21-1
          I19 = I21-2
          J20 = J21-1
          J19 = J21-2
      WRITE(*,*)'INPUT REYNOLDS NUMBER RE & TIME INCREMENT DT(100,.01)'
       READ(*,*) RE,DT
          TD  = 1./DT
          R1  = 1./RE
          R2  = 1./(RE*PRN)
      WRITE(*,*) 'INPUT HEAT PARAMETER --AR--  (5.0)'
       READ(*,*) AR
      WRITE(*,*) 'INPUT NUMBER OF TIME STEP (400)'
       READ(*,*) LM
      WRITE(*,*)'INPUT MAXIMUM NUMBER OF ITERATION FOR POISSON EQ.(20)'
       READ(*,*) KM
      WRITE(*,*) 'INPUT MAXIMUM ERROR EPS   (0.00001)'
       READ(*,*) EPS
      WRITE(*,*) 
     1 'INPUT CONCENTRATION PARAMETER FOR X & Y DIRECTION? (.98,.98)'
       READ(*,*) BX,BY
      WRITE(*,*) 'INPUT JA,JB -- ENTRANCE(JA-NX) EXIT(1-JB) (17,5)'
       READ(*,*) JA,JB
C
C***  GRID
               FA   = (EXP(BX)+1.)/(EXP(BX)-1.)
               FB   = (EXP(BY)+1.)/(EXP(BY)-1.)
           DO 110 I = 1,I21 
               BXA  = BX*FLOAT(I-1)/FLOAT(I21-1)
               X(I) = FA*(EXP(BXA)-1.)/(EXP(BXA)+1.)
  110      CONTINUE
           X00      = 2.*X(1)-X(2)
           X(I21+1) = 2.*X(I21)-X(I21-1)
           DO 120 J = 1,J21 
               BYA  = BY*FLOAT(J-1)/FLOAT(J21-1)
               Y(J) = FB*(EXP(BYA)-1.)/(EXP(BYA)+1.)
  120      CONTINUE
           Y00      = 2.*Y(1)-Y(2)
           Y(J21+1) = 2.*Y(J21)-Y(J21-1)
C
C***  CALCULATE METRICS (ONE-DIMENSION)
           DO 130 I = 1,I21
              X1    = X(I+1)-X(I)
              IF(I.GT.1) THEN
                 X2    = X(I)-X(I-1)
                 X3    = X(I+1)-X(I-1)
                 X4    = X(I+1)-2.*X(I)+X(I-1)
              ELSE 
                 X2    = X(I)-X00
                 X3    = X(I+1)-X00
                 X4    = X(I+1)-2.*X(I)+X00
              END IF
              A1(I) =-X1/(X2*X3)
              B1(I) = X4/(X1*X2)
              C1(I) = X2/(X1*X3)
              A2(I) = 2./(X2*X3)
              B2(I) =-2./(X1*X2)
              C2(I) = 2./(X1*X3)
  130      CONTINUE
           DO 140 J = 1,J21
              Y1    = Y(J+1)-Y(J)
              IF(J.GT.1) THEN
                 Y2    = Y(J)-Y(J-1)
                 Y3    = Y(J+1)-Y(J-1)
                 Y4    = Y(J+1)-2.*Y(J)+Y(J-1)
              ELSE
                 Y2    = Y(J)-Y00
                 Y3    = Y(J+1)-Y00
                 Y4    = Y(J+1)-2.*Y(J)+Y00
              END IF
              A3(J) =-Y1/(Y2*Y3)
              B3(J) = Y4/(Y1*Y2)
              C3(J) = Y2/(Y1*Y3)
              A4(J) = 2./(Y2*Y3)
              B4(J) =-2./(Y1*Y2)
              C4(J) = 2./(Y1*Y3)
  140      CONTINUE
C
C***  INITIAL CONDITION
           DO 20 J  = 1,J21
           DO 20 I  = 1,I21
             U(I,J) = 0.
             V(I,J) = 0.
             P(I,J) = 0.
             T(I,J) = 0.
   20      CONTINUE
C
C***  MAIN LOOP
      DO 30 L = 1,LM
C
C***  BOUNDARY CONDITION FOR VELOCITY & TEMPERATURE
           DO 40 I    = 1,I21
             U(I,J21) = 0.
             V(I,J21) = 0.
             V(I,1)   = 0.
             U(I,1)   = 0.
             T(I,J21) = 0.
             T(I,1)   = 1.
   40      CONTINUE
           DO 41 J    = JB+1,J21
             U(I21,J) = 0.
             V(I21,J) = 0.
             T(I21,J) = T(I20,J)
   41      CONTINUE
           DO 43 J    = 1,JA-1
             U(1,J)   = 0.
             V(1,J)   = 0.
             T(1,J)   = T(2,J)
   43      CONTINUE
           DO 44 J    = JA,J21-1
             U(1,J)   = 1.
             V(1,J)   = 0.
             T(1,J)   = 0.5
   44      CONTINUE
           DO 45 J    = 1,JB
             U(I21,J) = U(I20,J)*0.+1.
             V(I21,J) = V(I20,J)*0.
             T(I21,J) = T(I20,J)
   45      CONTINUE
C
C***  CALCULATIN OF RIGHT HAND SIDE OF POISSON EQ.
           DIVV = 0.
           DO 22 J    = 2,J20
           DO 22 I    = 2,I20
                 U1   = A1(I)*U(I-1,J)+B1(I)*U(I,J)+C1(I)*U(I+1,J)
                 U2   = A3(J)*U(I,J-1)+B3(J)*U(I,J)+C3(J)*U(I,J+1)
                 V1   = A1(I)*V(I-1,J)+B1(I)*V(I,J)+C1(I)*V(I+1,J)
                 V2   = A3(J)*V(I,J-1)+B3(J)*V(I,J)+C3(J)*V(I,J+1)
                 TY   = A3(J)*T(I,J-1)+B3(J)*T(I,J)+C3(J)*T(I,J+1)
               R(I,J) =-U1*U1-2.*U2*V1-V2*V2+TD*(U1+V2)+AR*TY
   22      CONTINUE
C
C***  SOLVING POIISON EQUATION FOR PRESSURE BY SOR METHOD
       DO 50 K = 1,KM
           G2  = 0.
C***  BOUNDARY CONDITION FOR PRESSURE
           DO 51 J    = 1,J21
             P(1,J)   = P(2,J)
             P(I21,J) = P(I20,J)
   51      CONTINUE
           DO 52 I    = 1,I21
             P(I,1)   = P(I,2)
             P(I,J21) = P(I,J20)
   52      CONTINUE
           DO 53 J    = 2,J20
           DO 53 I    = 2,I20
              ULI     = -(A2(I)*P(I-1,J)+C2(I)*P(I+1,J)+A4(J)*P(I,J-1)
     1                  +C4(J)*P(I,J+1)-R(I,J))/(B2(I)+B4(J)) -P(I,J)
              G2      = G2+ULI*ULI
              P(I,J)  = ULI+P(I,J)
   53      CONTINUE
             IF( G2.LE.EPS) GO TO 54
   50  CONTINUE
   54  CONTINUE
             IF(MOD(L,2).EQ.0) WRITE(*,*) L,K,G2
C
C***  TIME INTEGRATION OF NS EQUATION
       DO 55 J = 2,J20
       DO 55 I = 2,I20
           UN  = U(I,J)*(A1(I)*U(I-1,J)+B1(I)*U(I,J)+C1(I)*U(I+1,J))
     1          +V(I,J)*(A3(J)*U(I,J-1)+B3(J)*U(I,J)+C3(J)*U(I,J+1))
           UV  = A2(I)*U(I-1,J)+B2(I)*U(I,J)+C2(I)*U(I+1,J)
     1          +A4(J)*U(I,J-1)+B4(J)*U(I,J)+C4(J)*U(I,J+1)
           PX  = A1(I)*P(I-1,J)+B1(I)*P(I,J)+C1(I)*P(I+1,J) 
C
        D(I,J) = U(I,J)+DT*(-UN-PX+R1*UV)
   55  CONTINUE
       DO 56 J = 2,J20
       DO 56 I = 2,I20
           VN = U(I,J)*(A1(I)*V(I-1,J)+B1(I)*V(I,J)+C1(I)*V(I+1,J))
     1         +V(I,J)*(A3(J)*V(I,J-1)+B3(J)*V(I,J)+C3(J)*V(I,J+1))
           VV = A2(I)*V(I-1,J)+B2(I)*V(I,J)+C2(I)*V(I+1,J)
     1         +A4(J)*V(I,J-1)+B4(J)*V(I,J)+C4(J)*V(I,J+1)
           PY = A3(J)*P(I,J-1)+B3(J)*P(I,J)+C3(J)*P(I,J+1)
C
       R(I,J) = V(I,J)+DT*(-VN-PY+R1*VV+AR*T(I,J))
   56 CONTINUE
C
        DO 57 J = 2,J20
        DO 57 I = 2,I20
         U(I,J) = D(I,J)
   57   CONTINUE
        DO 58 J = 2,J20
        DO 58 I = 2,I20
         V(I,J) = R(I,J)
   58   CONTINUE
C
       DO 85 J = 2,J20
       DO 85 I = 2,I20
           TN = U(I,J)*(A1(I)*T(I-1,J)+B1(I)*T(I,J)+C1(I)*T(I+1,J))
     1         +V(I,J)*(A3(J)*T(I,J-1)+B3(J)*T(I,J)+C3(J)*T(I,J+1))
           TV = A2(I)*T(I-1,J)+B2(I)*T(I,J)+C2(I)*T(I+1,J)
     1         +A4(J)*T(I,J-1)+B4(J)*T(I,J)+C4(J)*T(I,J+1)
C
       D(I,J) = T(I,J)+DT*(-TN+R2*TV)
   85 CONTINUE
C
        DO 59 J = 2,J20
        DO 59 I = 2,I20
         T(I,J) = D(I,J)
   59   CONTINUE
   30 CONTINUE
C
      WRITE(*,*) '---U---'
      WRITE(*,600) ((U(I,J),I=2,20,2),J=2,20,2)
      WRITE(*,*) '---V---'
      WRITE(*,600) ((V(I,J),I=2,20,2),J=2,20,2)
      WRITE(*,*) '---P---'
      WRITE(*,600) ((P(I,J),I=2,20,2),J=2,20,2)
  600 FORMAT(1H ,10F8.2)
C
C***  CALCULATION OF STREAM FUNCTION FOR OUTPUT
        DO 66 I = 1,I21
         D(I,1) = 0.
   66   CONTINUE
        DO 65 J = 2,J21
        DO 65 I = 1,I21
         D(I,J) = D(I,J-1)+(Y(J)-Y(J-1))*U(I,J)
   65   CONTINUE
         CALL OUT(T,NX,NY,I21,J21)
         CALL OUT(D,NX,NY,I21,J21)
C
      WRITE(*,*) 'Save data? Yes=1,No=0'
       READ(*,*) ISAVE
        IF(ISAVE.EQ.1) THEN
          DO 67 J = 1,J21
          DO 67 I = 1,I21
             WRITE(17,*) U(I,J),V(I,J),P(I,J),T(I,J)
   67     CONTINUE
        END IF
      STOP
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
  601   FORMAT(1H ,80A2)
        WRITE(*,*) 'INPUT ANY NUMBER (1) '
         READ(*,*) AA
C
        DO 20 J = MY,1,-1
          DO 30 I = 1,MX
           INDEX(I) = INT((U(I,J)-UMIN)/(UMAX-UMIN)*9.9999)*11
   30     CONTINUE
           WRITE(*,600) (INDEX(I),I=1,MX)
   20   CONTINUE
  600   FORMAT(1H ,80I2)
C
      RETURN
      END
