C**********************************************************************
C     PROGRAM FOR SOLVING POTENTIAL FLOW PLOBLEM                      *
C**********************************************************************
      PROGRAM POTEN
      PARAMETER(JDIM=62,KDIM=31)
      COMMON/VEL/ U(JDIM,KDIM),V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, ITYP
      COMMON/FLWPRM/  EPS, ERR, CONST, ALP
C
           CALL GRID
           CALL DATA
           CALL METJ
      DO 10 KK = 1,3
           CALL INIT(KK)
           DO 50 I = 1,ISTEP0
               ERR = 0.
               CALL BC(KK)
               CALL PSI
               IF(ERR.LT.EPS) GO TO 60
               IF(MOD(I,2).EQ.0) WRITE(*,*) KK,I,ERR
   50      CONTINUE
   60      WRITE(*,*) ERR
           IF(KK.EQ.1) THEN
              DO 20 K = 1,KMAX
              DO 20 J = 1,JMAX
               U(J,K) = P(J,K)
   20         CONTINUE
           END IF
           IF(KK.EQ.2) THEN
              DO 25 K = 1,KMAX
              DO 25 J = 1,JMAX
               V(J,K) = P(J,K)
   25        CONTINUE
           END IF
   10 CONTINUE
           CALL SUP
           WRITE(*,*) 'Save data? Yes=1 No=0  '
            READ(*,*) ISAVE
            IF(ISAVE.EQ.1) CALL OUTP
            CALL OUT2
      STOP
      END
C**********************************************************************
C     SUBROUTINE FOR DETERMING  AND CALCULATING PARAMETERS            *
C**********************************************************************
      SUBROUTINE DATA
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, ITYP
      COMMON/FLWPRM/  EPS, ERR, CONST, ALP
        JM     = JMAX-1
        KM     = KMAX-1
        ISTEP0 = 10
        NSTEPS = 100
        EPS    = 0.001
        CONST  = 1.0
      WRITE(*,*) 'INPUT NUMBER OF ITERATION  (100)  '
       READ(*,*) ISTEP0
      WRITE(*,*) 'INPUT ATTACK ANGLE  (30)  '
       READ(*,*) ALPP
      WRITE(*,*) 'INPUT MAXIMUM ERROR  (0.00001)  '
       READ(*,*) EPS
        PAI    = ATAN(1.)*4.
        ALP    = PAI*ALPP/180.
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR CALCULATING OR READING GRIDS                     *
C**********************************************************************
      SUBROUTINE GRID
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON/GRIDXY/ X(JDIM,KDIM),Y(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0,  ITYP
      DIMENSION RR(KDIM)
C
      WRITE(*,*) 'DO YOU WANT TO READ FILE ? (YES=1,NO=0) (0) '
       READ(*,*) IGR
C
      IF (IGR.EQ.1) THEN
C
         READ(10,*) IA,IB
         READ(10,*) ((X(I,J),I=1,IA),J=1,IB),
     1              ((Y(I,J),I=1,IA),J=1,IB)
         JMAX = IA
         KMAX = IB
         WRITE(*,*) 'Periodic ? (Yes=1, No=0)'
          READ(*,*) ITYP 
C
      ELSE
C
C***  
         ITYP = 1
         JMAX = 62
         KMAX = 31
         PAI  = ATAN(1.)*4.
         AA   = 1.0
         BB   = 0.4
   99    WRITE(*,*) 'RR=A+A*R+A*R**2+A*R**3+...INPUT A & R? (0.1,1.1) '
          READ(*,*) HH,RA
          RR(1)  = 1.
         DO 15 K = 2,KMAX
           RR(K) = RR(K-1)+HH*RA**(K-1)
   15    CONTINUE
           WRITE(*,*) (RR(K),K=1,KMAX)
           WRITE(*,*) 'OK?  YES=1  NO=0  (1)  '
            READ(*,*) IOK
            IF(IOK.NE.1) GO TO 99
         DO 20 K = 1,KMAX
         DO 20 J = 1,JMAX
              TT = PAI*2.*FLOAT(J-2)/FLOAT(JMAX-2)
C             TT = PAI*FLOAT(J-1)/FLOAT(JMAX-1)
              BC = BB+(AA-BB)*FLOAT(K-1)/FLOAT(KMAX-1)
          X(J,K) = AA*RR(K)*COS(TT)
          Y(J,K) = BC*RR(K)*SIN(TT)
   20    CONTINUE
      END IF
C
      RETURN
      END
C
C********************************************************************
C     SUBROUTINE FOR CALCULATING METRICS                            *
C********************************************************************
      SUBROUTINE METJ
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON/GRIDXY/ X(JDIM,KDIM),Y(JDIM,KDIM)
      COMMON/METRIC/
     1 XX(JDIM,KDIM),XY(JDIM,KDIM),YX(JDIM,KDIM),YY(JDIM,KDIM),
     2 C1(JDIM,KDIM),C2(JDIM,KDIM),C3(JDIM,KDIM),
     3 C4(JDIM,KDIM),C5(JDIM,KDIM),AJ(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, ITYP
      DO 125 K=1,KMAX
      DO 125 J=1,JMAX
      IF(K.EQ.1) THEN
        XE=0.5*(-X(J,3)+4.0*X(J,2)-3.0*X(J,1))
        YE=0.5*(-Y(J,3)+4.0*Y(J,2)-3.0*Y(J,1))
      ELSE IF(K.EQ.KMAX) THEN
       XE=0.5*(X(J,KMAX-2)-4.0*X(J,KMAX-1)+3.0*X(J,KMAX))
       YE=0.5*(Y(J,KMAX-2)-4.0*Y(J,KMAX-1)+3.0*Y(J,KMAX))
      ELSE
        XE=0.5*(X(J,K+1)-X(J,K-1))
        YE=0.5*(Y(J,K+1)-Y(J,K-1))
      END IF
      IF(J.EQ.1) THEN
        XXI=0.5*(-X(3,K)+4.0*X(2,K)-3.0*X(1,K))
        YXI=0.5*(-Y(3,K)+4.0*Y(2,K)-3.0*Y(1,K))
        IF(ITYP.EQ.1) THEN
          XXI=0.5*(X(2,K)-X(JMAX-2,K))
          YXI=0.5*(Y(2,K)-Y(JMAX-2,K))
        END IF
      ELSE IF(J.EQ.JMAX) THEN
       XXI=0.5*(X(JMAX-2,K)-4.0*X(JMAX-1,K)+3.0*X(JMAX,K))
       YXI=0.5*(Y(JMAX-2,K)-4.0*Y(JMAX-1,K)+3.0*Y(JMAX,K))
        IF(ITYP.EQ.1) THEN
          XXI=0.5*(X(3,K)-X(JMAX-1,K))
          YXI=0.5*(Y(3,K)-Y(JMAX-1,K))
        END IF
      ELSE
        XXI=0.5*(X(J+1,K)-X(J-1,K))
        YXI=0.5*(Y(J+1,K)-Y(J-1,K))
      END IF
      IF(ITYP.EQ.1.AND.J.EQ.1) THEN
        XXI=0.5*(X(J+1,K)-X(JM-1,K))
        YXI=0.5*(Y(J+1,K)-Y(JM-1,K))
      END IF
      IF(ITYP.EQ.1.AND.J.EQ.JMAX) THEN
        XXI=0.5*(X(3,K)-X(J-1,K))
        YXI=0.5*(Y(3,K)-Y(J-1,K))
      END IF
        AJJ=XXI*YE-XE*YXI
      XX(J,K)=YE/AJJ
      YX(J,K)=-YXI/AJJ
      XY(J,K)=-XE/AJJ
      YY(J,K)=XXI/AJJ
      AJ(J,K)=AJJ
  125 CONTINUE
      DO 10 K=1,KMAX
      DO 10 J=1,JMAX
      C1(J,K)=XX(J,K)**2+XY(J,K)**2
      C3(J,K)=YX(J,K)**2+YY(J,K)**2
      C2(J,K)=2.*(XX(J,K)*YX(J,K)+XY(J,K)*YY(J,K))
   10 CONTINUE
      DO 11 K=2,KM
      DO 11 J=2,JM
      C77=XX(J,K)*(XX(J+1,K)-XX(J-1,K))+YX(J,K)*(XX(J,K+1)
     1   -XX(J,K-1))+XY(J,K)*(XY(J+1,K)-XY(J-1,K))
     2   +YY(J,K)*(XY(J,K+1)-XY(J,K-1))
      C88=XX(J,K)*(YX(J+1,K)-YX(J-1,K))+YX(J,K)*(YX(J,K+1)
     1   -YX(J,K))+XY(J,K)*(YY(J+1,K)-YY(J-1,K))
     2   +YY(J,K)*(YY(J,K+1)-YY(J,K-1))
      C4(J,K)=C77*.5
      C5(J,K)=C88*.5
   11 CONTINUE
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR INITIAL CONDITIONS                               *
C**********************************************************************
      SUBROUTINE INIT(KK)
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON/GRIDXY/ X(JDIM,KDIM),Y(JDIM,KDIM)
      COMMON /VEL/ U(JDIM,KDIM),V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, ITYP
      IF(KK.EQ.1) THEN
         DO 10  K = 1,KMAX
         DO 10  J = 1,JMAX
           P(J,K) = Y(J,K)
   10    CONTINUE
      END IF
      IF(KK.EQ.2) THEN
         DO 20  K = 1,KMAX
         DO 20  J = 1,JMAX
           P(J,K) = -X(J,K)
   20    CONTINUE
      END IF
      IF(KK.EQ.3) THEN
         DO 30  K = 1,KMAX
         DO 30  J = 1,JMAX
           P(J,K) = 0.
   30    CONTINUE
      END IF
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR SOLVING LAPLACE EQ. FOR STREAM FUNCTION          *
C**********************************************************************
      SUBROUTINE PSI
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON /VEL/ U(JDIM,KDIM),V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/METRIC/
     1 XX(JDIM,KDIM),XY(JDIM,KDIM),YX(JDIM,KDIM),YY(JDIM,KDIM),
     2 C1(JDIM,KDIM),C2(JDIM,KDIM),C3(JDIM,KDIM),
     3 C4(JDIM,KDIM),C5(JDIM,KDIM),AJ(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0,  ITYP
      COMMON/FLWPRM/ EPS, ERR, CONST, ALP
C
C***  SOR METHOD
      DO 10 K = 2,KM
      DO 10 J = 2,JM
        CC = .5/(C1(J,K)+C3(J,K))
        PA = C1(J,K)*(P(J+1,K)+P(J-1,K))+C3(J,K)*(P(J,K+1)+P(J,K-1))
     1     +.25*C2(J,K)*(P(J+1,K+1)-P(J-1,K+1)-P(J+1,K-1)+P(J-1,K-1))
     2     +.5*C4(J,K)*(P(J+1,K)-P(J-1,K))
     3     +.5*C5(J,K)*(P(J,K+1)-P(J,K-1))
        PP = PA*CC
        ERR    = ERR+(PP-P(J,K))**2
        P(J,K) = P(J,K)*(1.-CONST)+PP*CONST
   10 CONTINUE
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR BOUNDARY CONDITION                               *
C**********************************************************************
      SUBROUTINE BC(KK)
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON/GRIDXY/X(JDIM,KDIM),Y(JDIM,KDIM)
      COMMON/METRIC/
     1 XX(JDIM,KDIM),XY(JDIM,KDIM),YX(JDIM,KDIM),YY(JDIM,KDIM),
     2 C1(JDIM,KDIM),C2(JDIM,KDIM),C3(JDIM,KDIM),
     3 C4(JDIM,KDIM),C5(JDIM,KDIM),AJ(JDIM,KDIM)
      COMMON /VEL/ U(JDIM,KDIM),V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, ITYP
C
      IF(KK.EQ.1) THEN
         DO 20  J = 1,JMAX
           P(J,1)    = 0.
           P(J,KMAX) = Y(J,KMAX)
   20    CONTINUE
      END IF
C
      IF(KK.EQ.2) THEN
         DO 30  J = 1,JMAX
           P(J,1)    = 0.
           P(J,KMAX) = -X(J,KMAX)
   30 CONTINUE
      END IF
C
      IF(KK.EQ.3) THEN
         DO 40  J = 1,JMAX
           P(J,1)    = 1.
           P(J,KMAX) = 0.
   40    CONTINUE
      END IF
C
      DO 10 K = 1,KMAX
        P(1,K)    = P(JM,K)
        P(JMAX,K) = P(2,K)
   10 CONTINUE
C
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR SUPERPOSITION                                    *
C**********************************************************************
      SUBROUTINE SUP
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON /VEL/ U(JDIM,KDIM),V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/METRIC/
     1 XX(JDIM,KDIM),XY(JDIM,KDIM),YX(JDIM,KDIM),YY(JDIM,KDIM),
     2 C1(JDIM,KDIM),C2(JDIM,KDIM),C3(JDIM,KDIM),
     3 C4(JDIM,KDIM),C5(JDIM,KDIM),AJ(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, ITYP
      COMMON/FLWPRM/  EPS, ERR, CONST, ALP
C        DIMENSION VT(JDIM,2)
C        DO 10 J = 1,JMAX
C          G1 = SQRT(XY(J,1)**2+YY(J,1)**2)
C          G2 = SQRT(XY(J,KMAX)**2+YY(J,KMAX)**2)
C          P1 = U(J,1)*COS(ALP)+V(J,1)*SIN(ALP)+ALMD*P(J,1)
C          P2 = U(J,2)*COS(ALP)+V(J,2)*SIN(ALP)+ALMD*P(J,2)
C          VT(J,1) = G1*(P2-P1)
C          P3 = U(J,KM)*COS(ALP)+V(J,KM)*SIN(ALP)+ALMD*P(J,KM)
C          P4 = U(J,KMAX)*COS(ALP)+V(J,KMAX)*SIN(ALP)+ALMD*P(J,KMAX)
C          VT(J,2) = G2*(P4-P3)
C   10 CONTINUE
       CC = 0.
       G1 = SQRT(XY(3,2)**2+YY(3,2)**2)
       G2 = SQRT(XY(JM-1,2)**2+YY(JM-1,2)**2)
       AA = (G1*U(3,2)-G2*U(JM-2,2))*COS(ALP)
     1     +(G1*V(3,2)-G2*V(JM-2,2))*SIN(ALP)
       BB = G2*P(JM-2,2)-G1*P(3,2)
         IF(BB.NE.0.) CC=AA/BB
         WRITE(*,*) '--- LAMDA = ',CC
         IF(ITYP.EQ.1) CC=0
       DO 20 K = 1,KMAX
       DO 20 J = 1,JMAX
        P(J,K) = U(J,K)*COS(ALP)+V(J,K)*SIN(ALP)+CC*P(J,K)
   20  CONTINUE
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR WRITING DATA INTO FILE                           *
C**********************************************************************
      SUBROUTINE OUTP
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON /VEL/ U(JDIM,KDIM),V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, ITYP
         WRITE(12,*) JMAX,KMAX
         WRITE(12,*) ((P(J,K),J=1,JMAX),K=1,KMAX)
      RETURN
      END
C**********************************************************************
C     SUBPOUTINE FOR OUTPUT                                           *
C**********************************************************************
      SUBROUTINE OUT2
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON/GRIDXY/ X(JDIM,KDIM),Y(JDIM,KDIM)
      COMMON /VEL/ U(JDIM,KDIM),V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, ITYP
      CHARACTER*2 Z(39,21)
      DIMENSION INDEX(39,21)
C
      AREA(P1,Q1,P2,Q2,P3,Q3)=ABS((P2-P1)*(Q3-Q1)-(P3-P1)*(Q2-Q1))
      PMIN=P(1,1)
      DO 10 K=1,KMAX
      DO 10 J=1,JMAX
        IF(P(J,K).LT.PMIN) PMIN=P(J,K)
   10 CONTINUE
      PMAX=P(1,1)
      DO 30 K=1,KMAX
      DO 30 J=1,JMAX
        IF(P(J,K).GT.PMAX) PMAX=P(J,K)
   30 CONTINUE
      WRITE(*,*) 'PMAX=',PMAX,'   PMIN=',PMIN
      JA=1
      KA=1
      DO 40 K2=1,21
        WRITE(*,*) K2
      DO 40 J2=1,39
        INDEX(J2,K2)=0
        XA=FLOAT(J2-19)*.15
        YA=FLOAT(K2-11)*.15*1.2
        DO 60 K=1,KMAX-1
        DO 60 J=1,JMAX-1
          X1=X(J,K)
          Y1=Y(J,K)
          X2=X(J+1,K)
          Y2=Y(J+1,K)
          X3=X(J+1,K+1)
          Y3=Y(J+1,K+1)
          X4=X(J,K+1)
          Y4=Y(J,K+1)
          SA=AREA(X1,Y1,X2,Y2,X3,Y3)
          SB=AREA(X1,Y1,X3,Y3,X4,Y4)
          S1=AREA(X1,Y1,X2,Y2,XA,YA)
          S2=AREA(X2,Y2,X3,Y3,XA,YA)
          S3=AREA(X3,Y3,X4,Y4,XA,YA)
          S4=AREA(X4,Y4,X1,Y1,XA,YA)
          IF(ABS(SA+SB-S1-S2-S3-S4).LE..0001) THEN
            F1=.5*(S4+S1)/(SA+SB)
            F2=.5*(S1+S2)/(SA+SB)
            F3=.5*(S2+S3)/(SA+SB)
            F4=.5*(S3+S4)/(SA+SB)
            PJK=F1*P(J+1,K+1)+F2*P(J,K+1)+F3*P(J,K)+F4*P(J+1,K)
            IND=INT((PJK-PMIN)/(PMAX-PMIN)*100)
            INDEX(J2,K2)=(MOD(IND,9)+1)*11
            GO TO 40
          END IF
   60   CONTINUE
   40 CONTINUE
      DO 50 K=21,1,-1
        WRITE(*,600) (INDEX(J,K),J=1,39)
   50 CONTINUE
      WRITE(*,*) 'INPUT ANY NUMBER (1) '
       READ(*,*) AA
      DO 70 K=1,21
      DO 70 J=1,39
       IND = INDEX(J,K)
       IF(IND.EQ.0) THEN
         Z(J,K)='// '
        ELSE
          IF(MOD(IND,2).EQ.0) Z(J,K)='**'
          IF(MOD(IND,2).EQ.1) Z(J,K)='..'
        END IF
   70  CONTINUE
       DO 75 K=21,1,-1
         WRITE(*,601) (Z(J,K),J=1,39)
   75   CONTINUE
  600   FORMAT(1H ,39I2)
  601   FORMAT(1H ,39A2)
      WRITE(*,*) 'PMAX=',PMAX,'   PMIN=',PMIN
      RETURN
      END
