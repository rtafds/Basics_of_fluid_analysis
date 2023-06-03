C**********************************************************************
C     UNSTEADY FLOW IN THE ROOM (NO HEAT TRANSFER)                    *
C           PSI-OMEGA METHOD     NON-UNIFORM MESH SIZE                *
C**********************************************************************
C
      PARAMETER(MX=51,MY=51,MXP=MX+1,MYP=MY+1)
      DIMENSION PSI(MX,MY),OMG(MX,MY),TMP(MX,MY),X(MXP),Y(MYP)
     1,         A1(MX),B1(MX),C1(MX),A2(MX),B2(MX),C2(MX)
     2,         A3(MY),B3(MY),C3(MY),A4(MY),B4(MY),C4(MY)
C
C***  READ PARAMETERS
  99  WRITE(*,*)'INPUT NUMBER OF MESH FOR X & Y DIRECTION(<51) (20,20)'
       READ(*,*) NA,NB
           NX = NA + 1
           NY = NB + 1
      WRITE(*,*)'INPUT REYNOLDS NUMBER RE & TIME INCREMENT DT(100,.01)'
       READ(*,*) RE,DT
      WRITE(*,*) 'INPUT NUMBER OF TIME STEP (400)'
       READ(*,*) NMAX
      WRITE(*,*) 'INPUT ACCELARATION PARAMETER FOR PSI. (1.0)'
       READ(*,*) CONST
      WRITE(*,*) 'INPUT MAXIMUM ERROR EPS (0.00001)  '
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
           DO 110 I = 1,NX 
               BXA  = BX*FLOAT(I-1)/FLOAT(NX-1)
               X(I) = FA*(EXP(BXA)-1.)/(EXP(BXA)+1.)
  110      CONTINUE
            X00    = 2.*X(1)-X(2)
            X(NX+1) = 2.*X(NX)-X(NX-1)
           DO 120 J = 1,NY 
               BYA  = BY*FLOAT(J-1)/FLOAT(NY-1)
               Y(J) = FB*(EXP(BYA)-1.)/(EXP(BYA)+1.)
  120      CONTINUE
            Y00    = 2.*Y(1)-Y(2)
            Y(NY+1) = 2.*Y(NY)-Y(NY-1)
C
C***  CALCULATE METRICS (ONE-DIMENSION)
           DO 130 I = 1,NX
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
           DO 140 J = 1,NY
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
              PS0   = Y(NY)-Y(JA)
C
C***  INITIAL CONDITION FOR PSI AND OMEGA
           DO 10  J = 1,NY
           DO 10  I = 1,NX
              PSI(I,J)  = 0.0
              OMG(I,J)  = 0.0
   10      CONTINUE
           DO 11 I  = 1,NX
              PSI(I,NY) = PS0
   11      CONTINUE
           DO 12 J  = JB,NY
              PSI(NX,J) = PS0
   12      CONTINUE
C
C***  MAIN LOOP
C
      DO 100 N = 1,NMAX
C
C***  BOUNDARY CONDITION (STEP1)
C***  LEFT AND RIGHT
           DO 20  J = 1,NY
             OMG(1,J)  = (A2(1)*C1(1)/A1(1)-C2(1))*PSI(2,J)
             OMG(NX,J) = (C2(NX)*A1(NX)/C1(NX)-A2(NX))*PSI(NX-1,J)
     1                  +(C2(NX)*B1(NX)/C1(NX)-B2(NX))*PS0
   20      CONTINUE
           DO 21  J = JA+1,NY
             OMG(1,J)  = 0.
             PSI(1,J)  = (Y(J)-Y(JA))/(Y(NY)-Y(JA))*PS0
   21      CONTINUE
           DO 22  J = 1,JB-1
             OMG(NX,J) = OMG(NX-1,J)
             PSI(NX,J) = (Y(J)-Y(1))/(Y(JB)-Y(1))*PS0
   22      CONTINUE
C***  BOTTOM AND TOP
           DO 30  I = 1,NX
             OMG(I,1)  = (A4(1)*C3(1)/A3(1)-C4(1))*PSI(I,2)
             OMG(I,NY) = (C4(NY)*A3(NY)/C3(NY)-A4(NY))*PSI(I,NY-1)
     1                  +(C4(NY)*B3(NY)/C3(NY)-B4(NY))*PS0
   30      CONTINUE
C
C***  CALCULATE NEW PSI (STEP2)
         DO 55  K = 1,100 
           DO 60 J = 2,NY-1
           DO 60 I = 2,NX-1
C
             TMP(I,J)  = PSI(I,J)
C
             RHS       = -(A2(I)*PSI(I-1,J)+C2(I)*PSI(I+1,J)
     1                    +A4(J)*PSI(I,J-1)+C4(J)*PSI(I,J+1)
     2                    +OMG(I,J))/(B2(I)+B4(J)) -PSI(I,J)
             PSI(I,J)  = PSI(I,J)+CONST*RHS
   60      CONTINUE
C
             ERR2      = 0.
           DO 70 J = 2,NY-1
           DO 70 I = 2,NX-1
                BB = ABS(PSI(I,J)-TMP(I,J))
                IF(BB.GE.ERR2) ERR2 = BB
   70      CONTINUE
           IF(ERR2.LT.EPS) GO TO 75
   55    CONTINUE
   75    CONTINUE
                  WRITE(*,*) N,K,ERR2
C
C***  CALCULATE NEW OMEGA (STEP3)
        DO 40 J = 2,NY-1
        DO 40 I = 2,NX-1
C
          TMP(I,J)=OMG(I,J)
C
          RHS=-(A3(J)*PSI(I,J-1)+B3(J)*PSI(I,J)+C3(J)*PSI(I,J+1))
     1        *(A1(I)*OMG(I-1,J)+B1(I)*OMG(I,J)+C1(I)*OMG(I+1,J))
     2        +(A3(J)*OMG(I,J-1)+B3(J)*OMG(I,J)+C3(J)*OMG(I,J+1))
     3        *(A1(I)*PSI(I-1,J)+B1(I)*PSI(I,J)+C1(I)*PSI(I+1,J))
     4        +(A2(I)*OMG(I-1,J)+B2(I)*OMG(I,J)+C2(I)*OMG(I+1,J) 
     5         +A4(J)*OMG(I,J-1)+B4(J)*OMG(I,J)+C4(J)*OMG(I,J+1))/RE
          OMG(I,J) = OMG(I,J)+DT*RHS
   40   CONTINUE
C
  100 CONTINUE
C***  END OF MAIN LOOP
        CALL OUT(PSI,MX,MY,NX,NY)
        WRITE(*,*) 'SAVING IN FILE?  YES=1'
        READ(*,*) ISAVE
          IF(ISAVE.EQ.1) THEN
            WRITE(7,*) NX,NY
            WRITE(7,*) ((PSI(I,J),I=1,NX),J=1,NY),
     1                 ((OMG(I,J),I=1,NX),J=1,NY)
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