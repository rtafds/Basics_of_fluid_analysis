C**********************************************************************
C     FLOW IN THE ROOM
C**********************************************************************
      PARAMETER (MX=21,MY=21,MZ=21)
      COMMON U(MX,MY,MZ),V(MX,MY,MZ),W(MX,MY,MZ),P(MX,MY,MZ)
     1      ,Q(MX,MY,MZ),D(MX,MY,MZ),T(MX,MY,MZ),R(MX,MY),UL(MX)
     1      ,X(MX),XA(MX),XB(MX),XC(MX),XD(MX),XG(MX)
     1      ,Y(MY),YA(MY),YB(MY),YC(MY),YD(MY),YG(MY)
     1      ,Z(MZ),ZA(MZ),ZB(MY),ZC(MZ),ZD(MZ),ZG(MZ)
C
C***  READ DATA
      WRITE(*,*) 'NUMBERS OF GRID (X,Y,Z DIRECTION)? (20,20,20)'
       READ(*,*) I20,J20,K20
      WRITE(*,*) 'TIME INCREMENT DT & REYNOLDS NUMBER RE? (.005,40)'
       READ(*,*) DT,RE
      WRITE(*,*) 'CONCENTRATION PARAM. FOR X,Y,Z DIRECTION? (.9,.9,.9)'
       READ(*,*) BX,BY,BZ
      WRITE(*,*) 'NUMBER OF TIME STEP'
       READ(*,*) LM
C
C*** CALCULATION OF PARAMETER
         I21  =  I20 + 1
         J21  =  J20 + 1
         K21  =  K20 + 1
         TD   =  1./DT
         R1   =  1./RE
         I19  =  I21-2
         J19  =  J21-2
         K19  =  K21-2
         KM   =  10
         BX   = ALOG((1.+BX)/(1.-BX))
         BY   = ALOG((1.+BY)/(1.-BY))
         BZ   = ALOG((1.+BZ)/(1.-BZ))
         FA   = (EXP(BX)+1.)/(EXP(BX)-1.)
         FB   = (EXP(BY)+1.)/(EXP(BY)-1.)
         FC   = (EXP(BZ)+1.)/(EXP(BZ)-1.)
         PR   = 0.71
         R2   = 1./(RE*PR)
         AR   = 0.2
C
C***  GRID
           DO 110 I = 2,I20 
               BXA  = BX*(2.*FLOAT(I-2)/FLOAT(I20-2)-1.)
               X(I) = .5*(1.+FA*(EXP(BXA)-1.)/(EXP(BXA)+1.))
  110      CONTINUE
            X(1)    = 2.*X(2)-X(3)
            X(I21)  = 2.*X(I20)-X(I19)
           DO 115 J = 2,J20 
               BYA  = BY*(2.*FLOAT(J-2)/FLOAT(J20-2)-1.)
               Y(J) = .5*(1.+FB*(EXP(BYA)-1.)/(EXP(BYA)+1.))
  115      CONTINUE
            Y(1)    = 2.*Y(2)-Y(3)
            Y(J21)  = 2.*Y(J20)-Y(J19)
           DO 120 K = 2,K20 
               BZA  = BZ*(2.*FLOAT(K-2)/FLOAT(K20-2)-1.)
               Z(K) = .5*(1.+FC*(EXP(BZA)-1.)/(EXP(BZA)+1.))
  120      CONTINUE
            Z(1)    = 2.*Z(2)-Z(3)
            Z(K21)  = 2.*Z(K20)-Z(K19)
C
C***  CALCULATION OF METRICS
           DO 15 I = 2,I20
             XG(I)  = 2./(X(I+1)-X(I-1))
             XA(I)  = XG(I)*XG(I)
             XB(I)  = (X(I+1)-2.*X(I)+X(I-1))*XG(I)**3
             XC(I)  = XA(I)-XB(I)*.5
             XD(I)  = XA(I)+XB(I)*.5
   15      CONTINUE
           DO 16 J = 2,J20
             YG(J)  = 2./(Y(J+1)-Y(J-1))
             YA(J)  = YG(J)*YG(J)
             YB(J)  = (Y(J+1)-2.*Y(J)+Y(J-1))*YG(J)**3
             YC(J)  = YA(J)-YB(J)*.5
             YD(J)  = YA(J)+YB(J)*.5
   16      CONTINUE
           DO 17 K = 2,K20
             ZG(K)  = 2./(Z(K+1)-Z(K-1))
             ZA(K)  = ZG(K)*ZG(K)
             ZB(K)  = (Z(K+1)-2.*Z(K)+Z(K-1))*ZG(K)**3
             ZC(K)  = ZA(K)-ZB(K)*.5
             ZD(K)  = ZA(K)+ZB(K)*.5
   17      CONTINUE
           DO 39 K = 2,K20
           DO 39 J = 2,J20
           DO 39 I = 2,I20
             D(I,J,K) = .5/(XA(I)+YA(J)+ZA(K))
   39      CONTINUE
C
C***  INITIAL CONDITION
           DO 20 K = 1,K21
           DO 20 J = 1,J21
           DO 20 I = 1,I21
             U(I,J,K) = 0.
             V(I,J,K) = 0.
             W(I,J,K) = 0.
             P(I,J,K) = 0.
   20      CONTINUE
C
C***  MAIN LOOP
       DO 30 L = 1,LM
C
C***  BOUNDARY CONDITION FOR VELOCITY
C***  X-Z
          DO 40 K = 1,K21
          DO 40 I = 1,I21
            U(I,J21,K) = -U(I,J19,K)
            V(I,J21,K) =  V(I,J19,K)
            W(I,J21,K) = -W(I,J19,K)
            U(I,1,K)   = -U(I,3,K)
            V(I,1,K)   =  V(I,3,K)
            W(I,1,K)   = -W(I,3,K)
            U(I,J20,K) =  0.
            V(I,J20,K) =  0.
            W(I,J20,K) =  0.
            T(I,J20,K) =  T(I,J19,K)
            U(I,2,K)   =  0.
            V(I,2,K)   =  0.
            W(I,2,K)   =  0.
            T(I,2,K)   =  T(I,3,K)
   40     CONTINUE
C***  Y-Z
          DO 41 K = 1,K21
          DO 41 J = 1,J21
            U(I21,J,K) = U(I19,J,K)
            V(I21,J,K) = -V(I19,J,K)
            W(I21,J,K) = -W(I19,J,K)
            U(I20,J,K) =  0.
            V(I20,J,K) =  0.
            W(I20,J,K) =  0.
            T(I20,J,K) =  T(I19,J,K)
            U(1,J,K)   =  U(3,J,K)
            V(1,J,K)   = -V(3,J,K)
            W(1,J,K)   = -W(3,J,K)
            U(2,J,K)   =  0.
            V(2,J,K)   =  0.
            W(2,J,K)   =  0.
            T(2,J,K)   =  T(3,J,K)
   41     CONTINUE
          J7   = 7
          J14  = 14
          K4   = 4
          K6   = 6
          K16  = 16
          K18  = 18
          DO 44 K = K16,K18
          DO 44 J = J7,J14
            U(1,J,K)  = 1.
            U(2,J,K)  = 1.
            T(1,J,K)  = 1.
            T(2,J,K)  = 1.
   44     CONTINUE
          DO 45 K = K4,K6
          DO 45 J = J7,J14
            U(2,J,K)  = -1.*(Z(K18)-Z(K16))/(Z(K6)-Z(K4))
            U(1,J,K)  = U(2,J,K)
   45     CONTINUE
C***  X-Y
          DO 43 J = 1,J21
          DO 43 I = 1,I21
            U(I,J,1)   =  U(I,J,3)
            V(I,J,1)   =  V(I,J,3)
            W(I,J,1)   = -W(I,J,3)
            U(I,J,2)   =  0.
            V(I,J,2)   =  0.
            W(I,J,2)   =  0.
            T(I,J,2)   =  1.
C            U(I,J,K21) =  2.- U(I,J,K19)
            U(I,J,K21) = U(I,J,K19)
            V(I,J,K21) =  V(I,J,K19)
            W(I,J,K21) = W(I,J,K19)
            U(I,J,K20) = 1.*0.
            V(I,J,K20) = 0.
            W(I,J,K20) = 0.
            T(I,J,K20) = 0.
   43     CONTINUE
C
C***  CALCULATION FOR RHS OF POISSON EQ.
          DO 22 K = 2,K20
          DO 22 J = 2,J20
          DO 22 I = 2,I20
               U1 = (U(I+1,J,K)-U(I-1,J,K))*.5*XG(I)
               U2 = (U(I,J+1,K)-U(I,J-1,K))*.5*YG(J)
               U3 = (U(I,J,K+1)-U(I,J,K-1))*.5*ZG(K)
               V1 = (V(I+1,J,K)-V(I-1,J,K))*.5*XG(I)
               V2 = (V(I,J+1,K)-V(I,J-1,K))*.5*YG(J)
               V3 = (V(I,J,K+1)-V(I,J,K-1))*.5*ZG(K)
               W1 = (W(I+1,J,K)-W(I-1,J,K))*.5*XG(I)
               W2 = (W(I,J+1,K)-W(I,J-1,K))*.5*YG(J)
               W3 = (W(I,J,K+1)-W(I,J,K-1))*.5*ZG(K)
               T3 = (T(I,J,K+1)-T(I,J,K-1))*.5*ZG(K)
         Q(I,J,K) = -U1*U1-V2*V2-W3*W3-2.*(U2*V1+V3*W2+W1*U3)
     1              +TD*(U1+V2+W3)-AR*T3
   22     CONTINUE
C
C***  SOLUTION OF POISSON EQ. FOR PRESSURE BY G-S METHOD
        DO 23 KK = 1,KM
             G1  = 0.
C
C***  BOUNDARY CONDITION FOR PRESSURE
C***  Y-Z
          DO 24 K = 1,K21
          DO 24 J = 1,J21
            P(2,J,K)   = P(3,J,K)
            P(I20,J,K) = P(I19,J,K)
   24     CONTINUE
C***  X-Z
          DO 25 K = 1,K21
          DO 25 I = 1,I21
            P(I,2,K)   = P(I,3,K)
            P(I,J20,K) = P(I,J19,K)
   25     CONTINUE
C***  X-Y
          DO 29 J = 1,J21
          DO 29 I = 1,I21
            P(I,J,2)   = P(I,J,3)
            P(I,J,K20) = P(I,J,K19)
   29     CONTINUE
C
C***  GAUSS-SEIDEL METHOD
          DO 26 K = 3,K19
          DO 26 J = 3,J19
          DO 26 I = 3,I19
             PCOR = D(I,J,K)*(XC(I)*P(I+1,J,K)+YC(J)*P(I,J+1,K)
     1             +XD(I)*P(I-1,J,K)+YD(J)*P(I,J-1,K)
     2             +ZC(K)*P(I,J,K+1)+ZD(K)*P(I,J,K-1)-Q(I,J,K))
     3             -P(I,J,K)
             G1   = G1 + PCOR*PCOR
             P(I,J,K) = P(I,J,K) + PCOR
   26     CONTINUE
          IF( G1.LE..001) GO TO 31
   23   CONTINUE
   31 WRITE(*,*) L,K,G1
C
C***  CALCURATION OF NEW VELOCITY (N-S EQUATION)
          DO 27 K = 3,K19
          DO 27 J = 3,J19
          DO 27 I = 3,I19
C
            UNX = U(I,J,K)*(U(I+1,J,K)-U(I-1,J,K))*.5
            UNY = V(I,J,K)*(U(I,J+1,K)-U(I,J-1,K))*.5
            UNZ = W(I,J,K)*(U(I,J,K+1)-U(I,J,K-1))*.5
C
            VNX = U(I,J,K)*(V(I+1,J,K)-V(I-1,J,K))*.5
            VNY = V(I,J,K)*(V(I,J+1,K)-V(I,J-1,K))*.5
            VNZ = W(I,J,K)*(V(I,J,K+1)-V(I,J,K-1))*.5
C
            WNX = U(I,J,K)*(W(I+1,J,K)-W(I-1,J,K))*.5
            WNY = V(I,J,K)*(W(I,J+1,K)-W(I,J-1,K))*.5
            WNZ = W(I,J,K)*(W(I,J,K+1)-W(I,J,K-1))*.5
C
            TNX = U(I,J,K)*(T(I+1,J,K)-T(I-1,J,K))*.5
            TNY = V(I,J,K)*(T(I,J+1,K)-T(I,J-1,K))*.5
            TNZ = W(I,J,K)*(T(I,J,K+1)-T(I,J,K-1))*.5
C
             UV = (U(I+1,J,K)-2.*U(I,J,K)+U(I-1,J,K))*XA(I)
     1           +(U(I,J+1,K)-2.*U(I,J,K)+U(I,J-1,K))*YA(J)
     2           +(U(I,J,K+1)-2.*U(I,J,K)+U(I,J,K-1))*ZA(K)
     3           -(U(I+1,J,K)-U(I-1,J,K))*.5*XB(I)
     4           -(U(I,J+1,K)-U(I,J-1,K))*.5*YB(J)
     5           -(U(I,J,K+1)-U(I,J,K-1))*.5*ZB(K)
             VV = (V(I+1,J,K)-2.*V(I,J,K)+V(I-1,J,K))*XA(I)
     1           +(V(I,J+1,K)-2.*V(I,J,K)+V(I,J-1,K))*YA(J)
     2           +(V(I,J,K+1)-2.*V(I,J,K)+V(I,J,K-1))*ZA(K)
     3           -(V(I+1,J,K)-V(I-1,J,K))*.5*XB(I)
     4           -(V(I,J+1,K)-V(I,J-1,K))*.5*YB(J)
     5           -(W(I,J,K+1)-W(I,J,K-1))*.5*ZB(K)
             WV = (W(I+1,J,K)-2.*W(I,J,K)+W(I-1,J,K))*XA(I)
     1           +(W(I,J+1,K)-2.*W(I,J,K)+W(I,J-1,K))*YA(J)
     2           +(W(I,J,K+1)-2.*W(I,J,K)+W(I,J,K-1))*ZA(K)
     3           -(W(I+1,J,K)-W(I-1,J,K))*.5*XB(I)
     4           -(W(I,J+1,K)-W(I,J-1,K))*.5*YB(J)
     5           -(W(I,J,K+1)-W(I,J,K-1))*.5*ZB(K)
             TV = (T(I+1,J,K)-2.*T(I,J,K)+T(I-1,J,K))*XA(I)
     1           +(T(I,J+1,K)-2.*T(I,J,K)+T(I,J-1,K))*YA(J)
     2           +(T(I,J,K+1)-2.*T(I,J,K)+T(I,J,K-1))*ZA(K)
     3           -(T(I+1,J,K)-T(I-1,J,K))*.5*XB(I)
     4           -(T(I,J+1,K)-T(I,J-1,K))*.5*YB(J)
     5           -(T(I,J,K+1)-T(I,J,K-1))*.5*ZB(K)
C
          U(I,J,K) = U(I,J,K)+DT*(-(UNX*XG(I)+UNY*YG(J)+UNZ*ZG(K))
     1              -(P(I+1,J,K)-P(I-1,J,K))*.5*XG(I)+R1*UV)
          V(I,J,K) = V(I,J,K)+DT*(-(VNX*XG(I)+VNY*YG(J)+VNZ*ZG(K))
     1              -(P(I,J+1,K)-P(I,J-1,K))*.5*YG(J)+R1*VV)
          W(I,J,K) = W(I,J,K)+DT*(-(WNX*XG(I)+WNY*YG(J)+WNZ*ZG(K))
     1              -(P(I,J,K+1)-P(I,J,K-1))*.5*ZG(K)+R1*WV-AR*T(I,J,K))
          T(I,J,K) = T(I,J,K)+DT*(-TNX*XG(I)-TNY*YG(J)-TNZ*ZG(K)+R2*TV)
   27     CONTINUE
   30   CONTINUE
C
C***  CALCULATION OF STREAM FUNCTION
        DO 28 K = 1,K21
        DO 28 I = 1,I21
          R(I,K) = U(I,11,K)
   28   CONTINUE
          CALL PSI(R,Z,MX,MZ,I21,K21)
        WRITE(*,*) 'INPUT ANY NUMBER'
         READ(*,*) AA
        DO 281 K = 1,K21
        DO 281 I = 1,I21
          R(I,K) = T(I,11,K)
  281   CONTINUE
        CALL OUT(R,MX,MZ,I21,K21)
        WRITE(*,*) 'INPUT ANY NUMBER'
         READ(*,*) AA
          WRITE(*,*) '---U---'
          WRITE(*,*) ((U(I,J,5),I=2,I20,4),J=2,J20,4)
          WRITE(*,*) '---W---'
          WRITE(*,*) ((W(I,J,5),I=2,I20,4),J=2,J20,4)
          WRITE(*,*) '---P---'
          WRITE(*,*) ((P(I,J,5),I=2,I20,4),J=2,J20,4)
C
C***   SAVE DATA
        WRITE(*,*) 'DO YOU WANT TO SAVE DATA? (YES=1 NO=0)  '
         READ(*,*) ISAVE
         IF (ISAVE.EQ.1) THEN
          WRITE(*,*) 'INPUT FILE NUMBER (3) '
           READ(*,*) IS
            WRITE(IS,*) I21,J21,K21
            DO 301 I=1,I21
              WRITE(IS,*) X(I)
  301       CONTINUE
            DO 302 J=1,J21
              WRITE(IS,*) Y(J)
  302       CONTINUE
            DO 303 K=1,K21
              WRITE(IS,*) Z(K)
  303       CONTINUE
            DO 300 K=1,K21
            DO 300 J=1,J21
            DO 300 I=1,I21
              WRITE(IS,*) U(I,J,K),V(I,J,K),W(I,J,K),P(I,J,K),T(I,J,K)
  300       CONTINUE
         END IF
        STOP
      END
C**********************************************************************
C     SUBROUTINE FOR CALCURATION OF STREAM FUNCTION
C**********************************************************************
      SUBROUTINE PSI(U,Y,NX,NY,I21,J21)
      DIMENSION U(NX,NY),Y(NY)
      DIMENSION R(51,51)
C
C***  CALCULATING PSI 
        DO 70 I = 1,I21
         R(I,1) = 0.
   70   CONTINUE
        DO 71 J = 2,J21
        DO 71 I = 1,I21
         R(I,J) = R(I,J-1)+U(I,J)*(Y(J)-Y(J-1))
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
