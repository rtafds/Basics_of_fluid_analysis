C**********************************************************************
C     FLOW PAST MANY OBSTACLES
C**********************************************************************
      PARAMETER (ID=31,JD=17,KD=17)
      COMMON U(ID,JD,KD),V(ID,JD,KD),W(ID,JD,KD),P(ID,JD,KD)
     1      ,Q(ID,JD,KD),D(ID,JD,KD),S(ID,JD,KD),R(ID,KD)
     1      ,IFL(ID,JD,KD),IBD(ID,JD,KD),IWT(ID,JD,KD)
     1      ,UL(ID),X(ID),XA(ID),XB(ID),XC(ID),XD(ID),XG(ID)
     1      ,Y(JD),YA(JD),YB(JD),YC(JD),YD(JD),YG(JD)
     1      ,Z(KD),ZA(KD),ZB(KD),ZC(KD),ZD(KD),ZG(KD)
C
C***  INPUT DATA
      WRITE(*,*) 'NUMBER OF MESH (X,Y,Z DIRECTION)? (30,16,16)'
       READ(*,*) I20,J20,K20
      WRITE(*,*) 'TIME INCREMENT DT & REYNOLDS NUMBER RE? (.01,100)'
       READ(*,*) DT,RE
      WRITE(*,*) 'NUMBER OF TIME STEPS? (200)'
       READ(*,*) LM
C
C***  CALCULATION OF PARAMETERS
          I21  =  I20 + 1
          J21  =  J20 + 1
          K21  =  K20 + 1
          I17  =  17
          I11  =  (I17+1)/2
          J11  =  (J21+1)/2
          K11  =  (K21+1)/2
          DDH  =  1./FLOAT(K21-3)
          TD   =  1./DT
          R1   =  1./RE
          I19  =  I21-2
          J19  =  J21-2
          K19  =  K21-2
          KM   =  10
          ITEST=  1
C
      IF(ITEST.EQ.0) THEN
C
      WRITE(*,*) 'CONCENTRATION PARAM. FOR X,Y,Z DIRECTION? (.9,.9,.9)'
       READ(*,*) BX,BY,BZ
         BX   = ALOG((1.+BX)/(1.-BX))
         BY   = ALOG((1.+BY)/(1.-BY))
         BZ   = ALOG((1.+BZ)/(1.-BZ))
         FA   = (EXP(BX)+1.)/(EXP(BX)-1.)
         FB   = (EXP(BY)+1.)/(EXP(BY)-1.)
         FC   = (EXP(BZ)+1.)/(EXP(BZ)-1.)
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
      ELSE
C
           DO 201 I = 1,I17
  201         X(I) = FLOAT(I-I11)/FLOAT(I17-I11)
           DO 202 J = 1,J21
  202         Y(J) = FLOAT(J-J11)/FLOAT(J21-J11)
           DO 203 K = 1,K21
  203         Z(K) = FLOAT(K-K11)/FLOAT(K21-K11)
           XL = X(I17) - X(I17-1)
           DO 204 I = I17+1,I21
              XL   = XL*1.1
              X(I) = X(I-1) + XL
  204      CONTINUE
C
      END IF
C
C*** CALCULATION OF METRICS
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
           DO 92 K = 2,K20
             ZG(K) = 2./(Z(K+1)-Z(K-1))
             ZA(K) = ZG(K)*ZG(K)
             ZB(K) = (Z(K+1)-2.*Z(K)+Z(K-1))*ZG(K)**3
             ZC(K) = ZA(K)-ZB(K)*.5
             ZD(K) = ZA(K)+ZB(K)*.5
   92      CONTINUE
C
C***  DEFINITION OF OBSTACLE (SPHERE)
           DO 220 K = 1,K21
           DO 220 J = 1,J21
           DO 220 I = 1,I21
             STL = SQRT(X(I)**2+Y(J)**2+Z(K)**2)
             IFL(I,J,K) = 1
             IF(STL.LT..4) IFL(I,J,K) = 0
               IBD(I,J,K) = 0
               IWT(I,J,K) = 0
  220      CONTINUE
           DO 205 K=4,7
           DO 205 J=4,7
           DO 205 I=I17+2,I17+5
              IFL(I,J,K)=0
  205      CONTINUE
C
C***  CALCULATION OF THE POSITION OF THE BOUNDARY
           DO 230 K = 1,K21
           DO 230 J = 1,J21
             DO 231 I = 1,I21-1
               IF(IFL(I+1,J,K)-IFL(I,J,K).EQ.-1) IBD(I,J,K) = 1
  231        CONTINUE
             DO 232 I = I21,2,-1
               IF(IFL(I,J,K)-IFL(I-1,J,K).EQ.1)  IBD(I,J,K) = 1
  232        CONTINUE
  230      CONTINUE
           DO 250 K = 1,K21
           DO 250 I = 1,I21
             DO 251 J = 1,J21-1
               IF(IFL(I,J+1,K)-IFL(I,J,K).EQ.-1) IBD(I,J,K) = 1
  251        CONTINUE
             DO 252 J = J21,2,-1
               IF(IFL(I,J,K)-IFL(I,J-1,K).EQ.1)  IBD(I,J,K) = 1
  252        CONTINUE
  250      CONTINUE
           DO 260 J = 1,J21
           DO 260 I = 1,I21
             DO 261 K = 1,K21-1
               IF(IFL(I,J,K+1)-IFL(I,J,K).EQ.-1) IBD(I,J,K) = 1
  261        CONTINUE
             DO 262 K = K21,2,-1
               IF(IFL(I,J,K)-IFL(I,J,K-1).EQ.1) IBD(I,J,K)=1
  262        CONTINUE
  260      CONTINUE
C
           DO 291 K = 1,K21
           DO 291 J = 1,J21
           DO 291 I = 1,I21
             IBD(I,J,K) = IBD(I,J,K)*IFL(I,J,K)
  291      CONTINUE
C
C***  CALCULATION OF WEIGHT FUNCTION FOR PRESSURE B.C.
          DO 270 K = 2,K21-1
          DO 270 J = 2,J21-1
          DO 270 I = 2,I21-1
            IQQ = IBD(I+1,J,K)+IBD(I-1,J,K)+IBD(I,J+1,K)+IBD(I,J-1,K)
     1           +IBD(I,J,K+1)+IBD(I,J,K-1)
            IQQ = IQQ*(1-IFL(I,J,K))
              IF(IQQ.EQ.1) IWT(I,J,K) = 60
              IF(IQQ.EQ.2) IWT(I,J,K) = 30
              IF(IQQ.EQ.3) IWT(I,J,K) = 20
              IF(IQQ.EQ.4) IWT(I,J,K) = 15
              IF(IQQ.EQ.5) IWT(I,J,K) = 12
              IF(IQQ.EQ.6) IWT(I,J,K) = 10
  270     CONTINUE
C
C***  INITIAL CONDITION
          DO 20 K = 1,K21
          DO 20 J = 1,J21
          DO 20 I = 1,I21
            U(I,J,K) = 1.
            V(I,J,K) = 0.
            W(I,J,K) = 0.
            P(I,J,K) = 0.
   20     CONTINUE
C
C***  MAIN LOOP
      DO 30 L = 1,LM
C
C***  BOUNDARY CONDITION FOR VELOCITY
C***  X-Z
        DO 40 K = 1,K21
        DO 40 I = 1,I21
          U(I,J21,K) =  U(I,J19,K)
          V(I,J21,K) = -V(I,J19,K)
          W(I,J21,K) =  W(I,J19,K)
          U(I,1,K)   =  U(I,3,K)
          V(I,1,K)   = -V(I,3,K)
          W(I,1,K)   =  W(I,3,K)
          U(I,J20,K) =  U(I,J19,K)
          V(I,J20,K) =  0.
          W(I,J20,K) =  W(I,J19,K)
          U(I,2,K)   =  U(I,3,K)
          V(I,2,K)   =  0.
          W(I,2,K)   =  W(I,3,K)
   40   CONTINUE
C***  Y-Z
        DO 41 K = 1,K21
        DO 41 J = 1,J21
          U(I21,J,K) =  U(I19,J,K)
          V(I21,J,K) =  V(I19,J,K)
          W(I21,J,K) =  W(I19,J,K)
          U(I20,J,K) =  U(I19,J,K)
          V(I20,J,K) =  V(I19,J,K)
          W(I20,J,K) =  W(I19,J,K)
   41   CONTINUE
C***  X-Y
        DO 43 J = 1,J21
        DO 43 I = 1,I21
          U(I,J,1)   =  U(I,J,K21-3)
          V(I,J,1)   =  V(I,J,K21-3)
          W(I,J,1)   =  W(I,J,K21-3)
          U(I,J,2)   =  U(I,J,K21-2)
          V(I,J,2)   =  V(I,J,K21-2)
          W(I,J,2)   =  W(I,J,K21-2)
          U(I,J,K21) =  U(I,J,4)
          V(I,J,K21) =  V(I,J,4)
          W(I,J,K21) =  W(I,J,4)
          U(I,J,K20) =  U(I,J,3)
          V(I,J,K20) = V(I,J,3)
          W(I,J,K20) = W(I,J,3)
   43   CONTINUE
C
C***  CALCULATION OF RHS OF POISSON EQUATION FOR PRESSURE
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
           Q(I,J,K) = -U1*U1-V2*V2-W3*W3-2.*(U2*V1+V3*W2+W1*U3)
     1                +TD*(U1+V2+W3)
   22    CONTINUE
         DO 15 K = 2,K20
         DO 15 J = 2,J20
         DO 15 I = 2,I20
           D(I,J,K) = .5/(XA(I)+YA(J)+ZA(K))
   15    CONTINUE
C
C***   SOLUTION OF POISSON EQUATION BY G-S METHOD
       DO 23 KK = 1,KM
          G1 = 0.
C
C***   BOUNDARY CONDITION FOR PRESSURE
C***   OBSTACLES
         DO 285 K = 2,K20
         DO 285 J = 2,J20
         DO 285 I = 2,I20
           S(I,J,K) = IBD(I+1,J,K)*P(I+1,J,K)+IBD(I-1,J,K)*P(I-1,J,K)
     1               +IBD(I,J+1,K)*P(I,J+1,K)+IBD(I,J-1,K)*P(I,J-1,K)
     1               +IBD(I,J,K+1)*P(I,J,K+1)+IBD(I,J,K-1)*P(I,J,K-1)
  285    CONTINUE
         DO 286 K = 2,K20
         DO 286 J = 2,J20
         DO 286 I = 2,I20
           P(I,J,K) = S(I,J,K)*(1-IFL(I,J,K))*IWT(I,J,K)/60.
     1               +P(I,J,K)*IFL(I,J,K)
  286    CONTINUE
C
C***  Y-Z
         DO 24 K = 1,K21
         DO 24 J = 1,J21
           P(2,J,K)   = P(3,J,K)
           P(I20,J,K) = P(I19,J,K)
   24    CONTINUE
C
C***  X-Z
         DO 25 K = 1,K21
         DO 25 I = 1,I21
           P(I,2,K)   = P(I,3,K)
           P(I,J20,K) = P(I,J19,K)
   25    CONTINUE
C
C***  X-Y
         DO 29 J = 1,J21
         DO 29 I = 1,I21
           P(I,J,2)   = P(I,J,K19)
           P(I,J,K20) = P(I,J,3)
   29    CONTINUE
C
C***  GAUSS-SEIDEL METHOD
         DO 26  K = 3,K19
         DO 26  J = 3,J19
         DO 216 I = 3,I19
            UL(I) = D(I,J,K)*(XC(I)*P(I+1,J,K)+YC(J)*P(I,J+1,K)
     1             +XD(I)*P(I-1,J,K)+YD(J)*P(I,J-1,K)
     2             +ZC(K)*P(I,J,K+1)+ZD(K)*P(I,J,K-1)-Q(I,J,K))-P(I,J,K)
  216    CONTINUE
         DO 217 I = 3,I19
  217      G1 = G1+UL(I)*UL(I)
         DO 218 I = 3,I19
  218      P(I,J,K) = UL(I) + P(I,J,K)
   26    CONTINUE
         IF( G1.LE..0001) GO TO 31
   23    CONTINUE
   31    WRITE(*,*) L,K,G1
C
C***  CALCULATION OF NEW VELOCITY (N-S EQ.)
         DO 27 K = 3,K19
         DO 27 J = 3,J19
         DO 27 I = 3,I19
             UNX = U(I,J,K)*(U(I+1,J,K)-U(I-1,J,K))/2.
     1            -ABS(U(I,J,K))*(U(I+1,J,K)-2.*U(I,J,K)+U(I-1,J,K))/2.
             UNY = V(I,J,K)*(U(I,J+1,K)-U(I,J-1,K))/2.
     1            -ABS(V(I,J,K))*(U(I,J+1,K)-2.*U(I,J,K)+U(I,J-1,K))/2.
             UNZ = W(I,J,K)*(U(I,J,K+1)-U(I,J,K-1))/2.
     1            -ABS(W(I,J,K))*(U(I,J,K+1)-2.*U(I,J,K)+U(I,J,K-1))/2.
C
             VNX = U(I,J,K)*(V(I+1,J,K)-V(I-1,J,K))/2.
     1            -ABS(U(I,J,K))*(V(I+1,J,K)-2.*V(I,J,K)+V(I-1,J,K))/2.
             VNY = V(I,J,K)*(V(I,J+1,K)-V(I,J-1,K))/2.
     1            -ABS(V(I,J,K))*(V(I,J+1,K)-2.*V(I,J,K)+V(I,J-1,K))/2.
             VNZ = W(I,J,K)*(V(I,J,K+1)-V(I,J,K-1))/2.
     1            -ABS(W(I,J,K))*(V(I,J,K+1)-2.*V(I,J,K)+V(I,J,K-1))/2.
C
             WNX = U(I,J,K)*(W(I+1,J,K)-W(I-1,J,K))/2.
     1            -ABS(U(I,J,K))*(W(I+1,J,K)-2.*W(I,J,K)+W(I-1,J,K))/2.
             WNY = V(I,J,K)*(W(I,J+1,K)-W(I,J-1,K))/2.
     1            -ABS(V(I,J,K))*(W(I,J+1,K)-2.*W(I,J,K)+W(I,J-1,K))/2.
             WNZ = W(I,J,K)*(W(I,J,K+1)-W(I,J,K-1))/2.
     1            -ABS(W(I,J,K))*(W(I,J,K+1)-2.*W(I,J,K)+W(I,J,K-1))/2.
C
             UV  = (U(I+1,J,K)-2.*U(I,J,K)+U(I-1,J,K))*XA(I)
     1            +(U(I,J+1,K)-2.*U(I,J,K)+U(I,J-1,K))*YA(J)
     2            +(U(I,J,K+1)-2.*U(I,J,K)+U(I,J,K-1))*ZA(K)
     3            -(U(I+1,J,K)-U(I-1,J,K))*.5*XB(I)
     4            -(U(I,J+1,K)-U(I,J-1,K))*.5*YB(J)
     5            -(U(I,J,K+1)-U(I,J,K-1))*.5*ZB(K)
             VV  = (V(I+1,J,K)-2.*V(I,J,K)+V(I-1,J,K))*XA(I)
     1            +(V(I,J+1,K)-2.*V(I,J,K)+V(I,J-1,K))*YA(J)
     2            +(V(I,J,K+1)-2.*V(I,J,K)+V(I,J,K-1))*ZA(K)
     3            -(V(I+1,J,K)-V(I-1,J,K))*.5*XB(I)
     4            -(V(I,J+1,K)-V(I,J-1,K))*.5*YB(J)
     5            -(V(I,J,K+1)-V(I,J,K-1))*.5*ZB(K)
             WV  = (W(I+1,J,K)-2.*W(I,J,K)+W(I-1,J,K))*XA(I)
     1            +(W(I,J+1,K)-2.*W(I,J,K)+W(I,J-1,K))*YA(J)
     2            +(W(I,J,K+1)-2.*W(I,J,K)+W(I,J,K-1))*ZA(K)
     3            -(W(I+1,J,K)-W(I-1,J,K))*.5*XB(I)
     4            -(W(I,J+1,K)-W(I,J-1,K))*.5*YB(J)
     5            -(W(I,J,K+1)-W(I,J,K-1))*.5*ZB(K)
C
         D(I,J,K) = U(I,J,K)+DT*(-(UNX*XG(I)+UNY*YG(J)+UNZ*ZG(K))
     1             -(P(I+1,J,K)-P(I-1,J,K))*.5*XG(I)+R1*UV)
         Q(I,J,K) = V(I,J,K)+DT*(-(VNX*XG(I)+VNY*YG(J)+VNZ*ZG(K))
     1             -(P(I,J+1,K)-P(I,J-1,K))*.5*YG(J)+R1*VV)
         S(I,J,K) = W(I,J,K)+DT*(-(WNX*XG(I)+WNY*YG(J)+WNZ*ZG(K))
     1            -(P(I,J,K+1)-P(I,J,K-1))*.5*ZG(K)+R1*WV)
   27    CONTINUE
         DO 245 K = 3,K19
         DO 245 J = 3,J19
         DO 245 I = 3,I19
           U(I,J,K) = D(I,J,K)*IFL(I,J,K)
           V(I,J,K) = Q(I,J,K)*IFL(I,J,K)
           W(I,J,K) = S(I,J,K)*IFL(I,J,K)
  245    CONTINUE
   30  CONTINUE
C
C
C***  CALCULATION OF STREAM FUNCTION
        DO 28 K = 1,K21
        DO 28 I = 1,I21
          R(I,K) = U(I,J11,K)
   28   CONTINUE
          CALL PSI(R,Z,ID,KD,I21,K21)
        WRITE(*,*) 'INPUT ANY NUMBER'
         READ(*,*) AA
        IWRITE=1
          IF(IWRITE.EQ.1) THEN
             WRITE(*,*) '---P---'
             WRITE(*,*) ((P(I,J11,K),I=2,I21,4),K=2,K21,4)
             WRITE(*,*) '---U---'
             WRITE(*,*) ((U(I,J11,K),I=2,I21,4),K=2,K21,4)
             WRITE(*,*) '---V---'
             WRITE(*,*) ((V(I,J11,K),I=2,I21,4),K=2,K21,4)
             WRITE(*,*) '---W---'
             WRITE(*,*) ((W(I,J11,K),I=2,I21,4),J=2,K21,4)
          END IF
C
C*****  SAVE DATA
         WRITE(*,*) 'DO YOU WANT TO SAVE DATA? (Yes=1) '
         READ(*,*) ISV
         IF(ISV.EQ.1) THEN
            CALL SDATA(U,X,Y,Z,ID,JD,KD,I21,J21,K21,1)
            CALL SDATA(V,X,Y,Z,ID,JD,KD,I21,J21,K21,2)
            CALL SDATA(W,X,Y,Z,ID,JD,KD,I21,J21,K21,2)
            CALL SDATA(P,X,Y,Z,ID,JD,KD,I21,J21,K21,2)
         END IF
C
        STOP
      END
C      
      SUBROUTINE SDATA(D,X,Y,Z,ID,JD,KD,I21,J21,K21,IC)
      DIMENSION D(ID,JD,KD),X(ID),Y(JD),Z(KD)
      IF(IC.EQ.1) THEN
        WRITE(3,*) I21,J21,K21
        WRITE(3,*) (X(I),I=1,I21)
        WRITE(3,*) (Y(J),J=1,J21)
        WRITE(3,*) (Z(K),K=1,K21)
      END IF
      WRITE(3,*) (((D(I,J,K),I=1,I21),J=1,J21),K=1,K21)
      RETURN
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
        DO 71 I = 1,I21
        DO 71 J = 2,J21
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
