C**********************************************************************
C     STEADY FLOW IN CUBIC CAVITY --- PSI-OMEGA METHOD                *
C**********************************************************************
C
      PARAMETER(MX=51,MY=51)
      DIMENSION PSI(MX,MY),OMG(MX,MY),TMP(MX,MY)
C
C***  READ AND CALCULATE PARAMETERS
  99  WRITE(*,*) 'INPUT NUMBER OF MESH (<51)  (20)'
       READ(*,*) NA
         NN = NA + 1
         IF(NN.LE.1.OR.NN.GE.52) GO TO 99
       WRITE(*,*) 'INPUT REYNOLDS NUMBER (40)'
        READ(*,*) RE
       WRITE(*,*) 'INPUT MAXIMUM NUMBER OF ITERATION (500)'
        READ(*,*) NMAX
       WRITE(*,*) 'INPUT ACCELARATION TWO PARAMETERS (1.0,1.0)'
       WRITE(*,*) 'FOR Omega EQ. AND Psi EQ.  '
        READ(*,*) CONST1,CONST2
      WRITE(*,*) 'INPUT MAXIMUM ERROR EPS (0.00001)'
        READ(*,*) EPS
          H = 1./FLOAT(NN-1)
         HI = 1./H
C
C***  INITIAL CONDITION FOR PSI AND OMEGA
        DO 10 J = 1,NN
        DO 10 I = 1,NN
          PSI(I,J)=0.0
          OMG(I,J)=0.0
   10   CONTINUE
C
C***  MAIN LOOP
C
      DO 100 N = 1,NMAX
C
C***  BOUNDARY CONDITION (STEP1)
C***  LEFT AND RIGHT
        DO 20 J = 1,NN
          OMG(1,J)  = -2.*PSI(2,J)*HI*HI
          OMG(NN,J) = -2.*PSI(NN-1,J)*HI*HI
   20   CONTINUE
C***  BOTTOM AND TOP
        DO 30 I = 1,NN
          OMG(I,1)  = -2.*PSI(I,2)*HI*HI
          OMG(I,NN) = -2.*(PSI(I,NN-1)+H)*HI*HI
   30   CONTINUE
C
C***  CALCULATE NEW OMEGA (STEP2)
        DO 40 J = 2,NN-1
        DO 40 I = 2,NN-1
C
          TMP(I,J) = OMG(I,J)
C
          RHS = (OMG(I+1,J)+OMG(I-1,J)+OMG(I,J+1)+OMG(I,J-1))/4.
     1         +((PSI(I+1,J)-PSI(I-1,J))*(OMG(I,J+1)-OMG(I,J-1))
     2         -(PSI(I,J+1)-PSI(I,J-1))*(OMG(I+1,J)-OMG(I-1,J)))*RE/16.
          OMG(I,J) = OMG(I,J)*(1.-CONST1)+RHS*CONST1
   40   CONTINUE
C
           ERR1 = 0.
        DO 50 J = 2,NN-1
        DO 50 I = 2,NN-1
             BB = ABS(OMG(I,J)-TMP(I,J))
             IF(BB.GE.ERR1) ERR1 = BB
   50   CONTINUE
C
C***  CALCULATE NEW PSI (STEP3)
        DO 60 J = 2,NN-1
        DO 60 I = 2,NN-1
C
          TMP(I,J) = PSI(I,J)
C
          RHS = (PSI(I+1,J)+PSI(I-1,J)+PSI(I,J+1)+PSI(I,J-1))/4.
     1          +OMG(I,J)*H*H/4.
          PSI(I,J) = PSI(I,J)*(1.-CONST2)+RHS*CONST2
   60   CONTINUE
C
           ERR2 = 0.
        DO 70 J = 2,NN-1
        DO 70 I = 2,NN-1
             AA = AMAX1(1E-8,ABS(TMP(I,J)))
             BB = ABS(PSI(I,J)-TMP(I,J))/AA
             IF(BB.GE.ERR2) ERR2 = BB
   70   CONTINUE
C
       IF(MOD(N,50).EQ.0)
     1   WRITE(*,*) N,' ERROR(OMG)=',ERR1,'  ERROR(PSI)=',ERR2
        IF(AMAX1(ERR1,ERR2).LE.EPS) GO TO 80
C
  100 CONTINUE
C***  END OF MAIN LOOP
       WRITE(*,*) 'NOT CONVERGE!  DO YOU WANT CONTINUE? (YES=1)'
        READ(*,*) II
          IF(II.EQ.1) GO TO 99
   80   CALL OUT(PSI,MX,MY,NN,NN)
       WRITE(*,*) 'SAVING IN FILE?  YES=1  '
        READ(*,*) ISAVE
          IF(ISAVE.EQ.1) THEN
            WRITE(7,*) NN,NN
            WRITE(7,*) ((PSI(I,J),I=1,NN),J=1,NN),
     1                 ((OMG(I,J),I=1,NN),J=1,NN)
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
        WRITE(*,*) 'INPUT ANY NUMBER (1) '
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