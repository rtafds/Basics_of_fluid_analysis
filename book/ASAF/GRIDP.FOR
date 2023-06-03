C**********************************************************************
C     ELLIPTIC GRID GENERATION                                        *
C**********************************************************************
      PARAMETER(MX=42,MY=16)
      COMMON X(MX,MY),Y(MX,MY),P(MX,MY),Q(MX,MY)
      WRITE(*,*) 'INPUT RADIUS OF INNER AND OUTER CIRCLE (1.0,4.0)'
       READ(*,*) R1,R2
      WRITE(*,*) 'INPUT DISTANCE BETWEEN CENTERS OF CIRCLE (1.0)'
       READ(*,*) X3
C
         PAI = 4.*ATAN(1.)
         DO 10 K = 1,MY-1
         DO 10 I = 1,MX
           R   = R1+(R2-R1-X3-0.1)*FLOAT(K-1)/FLOAT(MY-1)
           TET = 2.0*PAI*FLOAT(I-1)/FLOAT(MX-2)
           X(I,K) = 0.6*R*COS(TET)
           Y(I,K) = R*SIN(TET)
   10    CONTINUE
         DO 20 I = 1,MX
           TET = 2.0*PAI*FLOAT(I-1)/FLOAT(MX-2)
           X(I,MY) = R2*COS(TET)-X3
           Y(I,MY) = R2*SIN(TET)
   20    CONTINUE
C
         CALL ELLIP
C
         WRITE(*,*) 'Save data? Yes=1, No=0'
          READ(*,*) ISAVE
           IF(ISAVE.EQ.1) THEN
            DO 23 K = 1,MY
            DO 23 I = 1,MX
              WRITE(7,*) X(I,K),Y(I,K)
   23       CONTINUE
           END IF
C
C***  PLOT
         WRITE(*,*) 'Plot? Yes=1,No=0'
          READ(*,*) IPLT
          IF(IPLT.EQ.1) THEN
C          
         WRITE(*,*) 'FCT,X11,Y11 (40,6,5)'
          READ(*,*) FCT,X11,Y11
c            CALL PLOTS(0)
c            CALL FACTOR(FCT)
c            CALL PLOT(0.,0.,-3)
            DO 30 K = 1,MY
            DO 30 I = 1,MX
               IC=2
               IF(I.EQ.1) IC=3
               X1 = X(I,K)+X11
               Y1 = Y(I,K)+Y11
c               CALL PLOT(X1,Y1,IC)
   30       CONTINUE
            DO 40 I = 1,MX
            DO 40 K = 1,MY
               IC=2
               IF(K.EQ.1) IC=3
               X1 = X(I,K)+X11
               Y1 = Y(I,K)+Y11
c               CALL PLOT(X1,Y1,IC)
   40       CONTINUE
c               CALL PLOT(0.,0.,-3)
c               CALL PLOT(0.,0.,999)
          END IF
      STOP
      END
C**********************************************************************
C     SUBROUTINE FOR SOLVING ELLIPTIC EQ.                             *
C**********************************************************************
      SUBROUTINE ELLIP
      PARAMETER(IM=42,KM=16)
      COMMON X(IM,KM),Y(IM,KM),P(IM,KM),Q(IM,KM)
        KD = IM
        LD = KM
      WRITE(*,*) 'INPUT NUMBER OF ITERATION (300)'
       READ(*,*) K100
      WRITE(*,*) 'TYPE OF GENERATION (=0 LAPLACE  '
      WRITE(*,*) '=1 CONCENTRATION -POINT  =2 CONCENTRATION -LINE'
       READ(*,*) ITYP
C
C***  CALCULATION OF SOURCE TERM       
        IF(ITYP.EQ.0) THEN
          DO 30  L = 1,LD
          DO 30  K = 1,KD
            P(K,L) = 0.
            Q(K,L) = 0.
   30     CONTINUE
      ELSE IF(ITYP.EQ.1) THEN
          WRITE(*,*)'INPUT CONCENTRATION(A) & DECAY(C) FACTOR (-1,.4)'
           READ(*,*) AA,CC
          WRITE(*,*) 'INPUT MESH POINT (K,L) (15,8)'
           READ(*,*) K4,L4
          DO 40  L = 1,LD
          DO 40  K = 1,KD
            SG     = 1.
            IF(K.GE.K4) SG = -1.
            P(K,L) =-AA*SG*EXP(-CC*SQRT(FLOAT(K-K4)**2+FLOAT(L-L4)**2))
            Q(K,L) = 0.
   40       CONTINUE
      ELSE
          WRITE(*,*)'INPUT CONCENTRATION(A) & DECAY(C) FACTOR (-.5,.5)'
           READ(*,*) AA,CC
          WRITE(*,*) 'INPUT LINE NUMBER (L) (8)'
           READ(*,*) L3
          DO 50  L = 1,LD
          DO 50  K = 1,KD
            P(K,L) = 0.
            SG     = 1.
            IF(L.GT.L3) SG = -1.
            Q(K,L) =-AA*SG*EXP(-CC*ABS(L-L3))
   50     CONTINUE
      END IF
             K30   = K100/3
             K60   = K30*2
             EPS1  = 0.000001
             CONST = 0.05
C
C***  SOLVING ELLIPTIC EQ. BY SOR METHOD
         DO 120 KK = 1,K100
              GOSA = 0.0
              IF(CONST.LE.1..AND.KK.GT.K30) CONST = 0.25
              IF(CONST.LE.1..AND.KK.GT.K60) CONST = 1.0
              IF(ITYP.EQ.0) CONST = 1.0
           DO 121 L = 2,LD-1
           DO 121 K = 2,KD-1
               XXI  = (X(K+1,L)-X(K-1,L))*.5
               YXI  = (Y(K+1,L)-Y(K-1,L))*.5
               XETA = (X(K,L+1)-X(K,L-1))*.5
               YETA = (Y(K,L+1)-Y(K,L-1))*.5
               ALPH = XETA*XETA+YETA*YETA
               BETA = XXI*XETA+YXI*YETA
               GAMA = XXI*XXI+YXI*YXI
               AJAC = XXI*YETA-XETA*YXI
               WX   = P(K,L)*ALPH/(AJAC*AJAC)
               QE   = Q(K,L)*GAMA/(AJAC*AJAC)
               AG   = 0.5/(ALPH+GAMA)
               SX   = (ALPH*(X(K+1,L)+X(K-1,L))-.5*BETA*(X(K+1,L+1)
     1               -X(K-1,L+1)-X(K+1,L-1)+X(K-1,L-1))+GAMA*(X(K,L+1)
     2               +X(K,L-1))+AJAC*AJAC*(WX*XXI+QE*XETA))*AG-X(K,L)
               SY   = (ALPH*(Y(K+1,L)+Y(K-1,L))-.5*BETA*(Y(K+1,L+1)
     1               -Y(K-1,L+1)-Y(K+1,L-1)+Y(K-1,L-1))+GAMA*(Y(K,L+1)
     2               +Y(K,L-1))+AJAC*AJAC*(WX*YXI+QE*YETA))*AG-Y(K,L)
             X(K,L) = X(K,L)+CONST*SX
             Y(K,L) = Y(K,L)+CONST*SY
             GOSA   = GOSA+SX*SX+SY*SY
  121      CONTINUE
C
C***   PERIODIC BOUNDARY CONDITION
           DO 123 L = 1,LD
            X(1,L)  = X(KD-1,L)
            Y(1,L)  = Y(KD-1,L)
            X(KD,L) = X(2,L)
            Y(KD,L) = Y(2,L)
  123      CONTINUE
           IF(GOSA.LT.EPS1) GO TO 122
           WRITE(*,*)KK,GOSA
  120     CONTINUE
  122   CONTINUE
      RETURN
      END
