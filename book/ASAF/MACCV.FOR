C**********************************************************************
C     UNSTEADY FLOW IN CUBIC CAVITY  --- MAC METHOD                   *
C**********************************************************************
      PARAMETER (NX=51,NY=51)
      COMMON U(NX,NY),V(NX,NY),P(NX,NY),R(NX,NY),D(NX,NY)
      WRITE(*,*) 'INPUT NUMBER OF MESH (<51) (20,20)'
       READ(*,*) MX,MY
         I21 = MX + 1
         J21 = MY + 1
      WRITE(*,*) 'INPUT TIME INCLEMENT DT  (0.01)'
       READ(*,*) DT
      WRITE(*,*) 'INPUT REYNOLDS NUMBER  (40)'
       READ(*,*) RE
      WRITE(*,*) 'INPUT NUMBER OF TIME STEP  (300)'
       READ(*,*) LM
      WRITE(*,*) 'INPUT MAXMUM NUMBER OF ITERATION FOR POISSON EQ.(50)'
       READ(*,*) KM
         DX = 1./FLOAT(I21-3)
         DY = 1./FLOAT(J21-3)
         XD = 1./DX
         YD = 1./DY
        XDH = .5*XD
        YDH = .5*YD
         TD = 1./DT
         R1 = 1./RE
        I20 = I21-1
        I19 = I21-2
        J20 = J21-1
        J19 = J21-2
         A1 = .5*DY*DY/(DX*DX+DY*DY)
         A2 = .5*DX*DX/(DX*DX+DY*DY)
         A3 = .5*DY*DY/(1.+DY*DY/(DX*DX))
C
C***  INITIAL CONDITION
        DO 20 J = 1,J21
        DO 20 I = 1,I21
         U(I,J) = 0.
         V(I,J) = 0.
         P(I,J) = 0.
   20   CONTINUE
C
C***   TIME MARCHING
        DO 30 L = 1,LM
C
C***   BOUNDARY CONDITION FOR VELOCITY
          DO 40 I = 1,I21
            U(I,J20) = 1.
            V(I,J21) = V(I,J19)
            V(I,J20) = 0.
            V(I,1)   = V(I,3)
            V(I,2)   = 0.
            U(I,1)   = -U(I,2)
   40     CONTINUE
          DO 41 J = 1,J21
            U(I21,J) = U(I19,J)
            U(I20,J) = 0.
            V(I20,J) = -V(I19,J)
            U(1,J)   = U(3,J)
            U(2,J)   = 0.
            V(1,J)   = -V(2,J)
   41     CONTINUE
C
C***  CALCULATE Right Hand Side OF POISSON EQUATION
          DIVV = 0.
          DO 22 J = 2,J19
          DO 22 I = 2,I19
               U1 = (U(I+1,J)-U(I,J))*XD
C              U2 = (U(I,J+1)-U(I,J-1))*XDH
C              V1 = (V(I+1,J)-V(I-1,J))*YDH
               V2 = (V(I,J+1)-V(I,J))*YD
           D(I,J) = U1+V2
             DIVV = DIVV+ABS(U1+V2)
               UA =.25*(U(I,J)+U(I+1,J)+U(I+1,J+1)+U(I,J+1))
               UB =.25*(U(I,J)+U(I+1,J)+U(I+1,J-1)+U(I,J-1))
               VA =.25*(V(I,J)+V(I,J+1)+V(I+1,J+1)+V(I+1,J))
               VB =.25*(V(I,J)+V(I,J+1)+V(I-1,J+1)+V(I-1,J))
           R(I,J) =-U1*U1-2.*(UA-UB)*(VA-VB)*XD*YD-V2*V2+TD*(U1+V2)
   22     CONTINUE
C
C***  SOLVING POISSON EQ. BY SOR METHOD
          DO 50 K = 1,KM
               G2 = 0.
C
C***  PRESSURE BOUNDARY CONDITION
            DO 51 J = 1,J21
             P(1,J)   = P(2,J)
             P(I20,J) = P(I19,J)
   51       CONTINUE
            DO 52 I = 1,I21
             P(I,1)   = P(I,2)
             P(I,J20) = P(I,J19)
   52       CONTINUE
C
C***  SOR METHOD 
             DO 53 J = 2,J19
             DO 53 I = 2,I19
               ULI = A1*(P(I+1,J)+P(I-1,J))+A2*(P(I,J+1)+P(I,J-1))
     1              -A3*R(I,J)-P(I,J)
                G2 = G2+ULI*ULI
            P(I,J) = ULI+P(I,J)
   53       CONTINUE
         IF( G2.LE..000001) GO TO 54
   50 CONTINUE
   54 continue
       IF(MOD(L,2).EQ.0) WRITE(*,*) L,K,G2
C
C*** TIME INTEGRATION FOR NS EQ.
            DO 60 J = 2,J19
            DO 60 I = 3,I19
                 V1 = .25*(V(I,J)+V(I,J+1)+V(I-1,J+1)+V(I-1,J))
                 UN = U(I,J)*(U(I+1,J)-U(I-1,J))*XDH
     1               +V1*(U(I,J+1)-U(I,J-1))*YDH
                 UV = (U(I+1,J)-2.*U(I,J)+U(I-1,J))*XD*XD
     1               +(U(I,J+1)-2.*U(I,J)+U(I,J-1))*YD*YD
C
             D(I,J) = U(I,J)+DT*(-UN-(P(I,J)-P(I-1,J))*XD+R1*UV)
   60       CONTINUE
            DO 61 J = 3,J19
            DO 61 I = 2,I19
                 U1 = .25*(U(I,J)+U(I+1,J)+U(I+1,J-1)+U(I,J-1))
                 VN = U1*(V(I+1,J)-V(I-1,J))*XDH
     1               +V(I,J)*(V(I,J+1)-V(I,J-1))*YDH
                 VV = (V(I+1,J)-2.*V(I,J)+V(I-1,J))*XD*XD
     1               +(V(I,J+1)-2.*V(I,J)+V(I,J-1))*YD*YD
C
             R(I,J) = V(I,J)+DT*(-VN-(P(I,J)-P(I,J-1))*YD+R1*VV)
   61       CONTINUE
            DO 62 J = 2,J19
            DO 62 I = 3,I19
             U(I,J) = D(I,J)
   62       CONTINUE
            DO 63 J = 3,J19
            DO 63 I = 2,I19
             V(I,J) = R(I,J)
   63       CONTINUE
   30 CONTINUE
      WRITE(*,*) '---U---'
      WRITE(*,600) ((U(I,J),I=2,20,2),J=2,20,2)
      WRITE(*,*) '---V---'
      WRITE(*,600) ((V(I,J),I=2,20,2),J=2,20,2)
      WRITE(*,*) '---P---'
      WRITE(*,600) ((P(I,J),I=2,20,2),J=2,20,2)
  600 FORMAT(1H ,10F6.3)
C
C***  CALCULATING PSI 
        DO 70 I = 1,I21
         R(I,1) = 0.
   70   CONTINUE
        DO 71 J = 2,J21
        DO 71 I = 1,I21
         R(I,J) = R(I,J-1)+DY*U(I,J)
   71   CONTINUE
           CALL OUT(R,NX,NY,I21,J21)
C
      WRITE(*,*) 'Save data? Yes=1, No=0'
       READ(*,*) ISAVE
        IF(ISAVE.EQ.1) THEN
        DO 80 J = 1,J21
        DO 80 I = 1,I21
          WRITE(17,*) U(I,J),V(I,J),P(I,J)
   80   CONTINUE
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
