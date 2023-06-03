C**********************************************************************
C     2 dimensional flow slover (genaralized coordinate sysytem)      *
C        Stream function - Vorticity (Psi-Omega) Method               *
C**********************************************************************
      PROGRAM FLOW2D
      PARAMETER(JDIM=62,KDIM=31)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR
      COMMON/FLWPRM/ RE, REI, DT, DTI, EPS, ERR, CONST
C
        CALL GRID
        CALL DATA
        CALL METJ
C
        CALL INIT
      DO 10 N = 1,NSTEPS
         ERR = 0.
         DO 20 I = 1,ISTEP0
            CALL PSI
            IF(ERR.LT.EPS) GO TO 30
   20    CONTINUE
   30    IF(MOD(N,2).EQ.0) WRITE(*,*) N ,ERR
        CALL OMEGA
        CALL BC
   10 CONTINUE
        WRITE(*,*) 'SAVE DATA? (Yes=1, No=0) '
         READ(*,*) ISAVE
         IF(ISAVE.EQ.1) CALL OUTP
        CALL OUT2
C
      STOP
      END
C********************************************************************
C     SUBROUTINE FOR INPUT DATA                                     *
C********************************************************************
      SUBROUTINE DATA
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR
      COMMON/FLWPRM/ RE, REI, DT, DTI, EPS, ERR, CONST
         JM = JMAX-1
         KM = KMAX-1
         WRITE(*,*) 'INPUT NUMBER OF TIME STEP (400) '
          READ(*,*) NSTEPS
         WRITE(*,*)'INPUT MAX. NUM.OF ITERATION FOR POISSON EQ.(10)'
          READ(*,*) ISTEP0
         WRITE(*,*)' REYNOLDS NUM. & TIME INCREMENT DT? (100,.005) '
          READ(*,*) RE,DT
         WRITE(*,*) 'INPUT MAX. ERROR  (.0001)'
          READ(*,*) EPS
         WRITE(*,*) 'ACCERARATION PARAMETER FOR SOR METHOD? (1)  '
          READ(*,*) CONST
         REI = 1./RE
         DTI = 1./DT
      RETURN
      END
C********************************************************************
C     SUBROUTINE FOR GRID DATA                                      *
C********************************************************************
      SUBROUTINE GRID
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON/GRIDXY/ X(JDIM,KDIM),Y(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR
      DIMENSION RR(KDIM)
C
      WRITE(*,*) 'DO YOU WANT TO READ FILE ? (YES=1,NO=2) (2)'
        READ(*,*) IGR
      WRITE(*,*) 'Periodic ? (Yes=1,No=2) (1)'
        READ(*,*) ITYP
C
      IF (IGR.EQ.1) THEN
C
        READ(38,*) IA,IB,IC
        READ(38,*) (((X(I,J),I=1,IA),J=1,IB),K=1,IC)
        READ(38,*) (((Y(I,J),I=1,IA),J=1,IB),K=1,IC)
        READ(38,*) (((RR(K),I=1,IA),J=1,IB),K=1,IC)
        JMAX = IA
        KMAX = IB
C
      ELSE
C
C***  
        JMAX = 62
        KMAX = 31
        PAI  = ATAN(1.)*4.
        AA   = .4
        BB   = 1.
        HH   = .05
        RA   = 1.125
        RR(1)= 1.
        DO 15 K = 2,KMAX
          RR(K) = RR(K-1)+HH*RA**(K-1)
   15   CONTINUE
        WRITE(*,*) (RR(K),K=1,KMAX)
        DO 20 K = 1,KMAX
        DO 20 J = 1,JMAX
             TT = 2.*PAI*FLOAT(J-2)/FLOAT(JMAX-2)
             BC = BB+(AA-BB)*FLOAT(K-1)/FLOAT(KMAX-1)
         X(J,K) = AA*RR(K)*COS(TT)
         Y(J,K) = BC*RR(K)*SIN(TT)
   20 CONTINUE
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
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR
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
C********************************************************************
C     SUBROUTINE FOR INITIAL CONDITION                              *
C********************************************************************
      SUBROUTINE INIT
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON/GRIDXY/ X(JDIM,KDIM),Y(JDIM,KDIM)
      COMMON /VEL/ V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR
      IF(ITYP.EQ.1) THEN
         DO 10  K = 1,KMAX
         DO 10  J = 1,JMAX
           V(J,K) = 0.
           P(J,K) = Y(J,K)
   10    CONTINUE
      ELSE
         DO 20  K = 1,KMAX
         DO 20  J = 1,JMAX
           V(J,K) = 0.
           P(J,K) = (Y(J,K)-Y(J,1))/(Y(J,KMAX)-Y(J,1))
   20    CONTINUE
       END IF
      RETURN
      END
C********************************************************************
C     SUBROUTINE FOR SOLVING POISSON EQUATION OF STREAM FUNCTION    *
C********************************************************************
      SUBROUTINE PSI
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON /VEL/ V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON /FLW/ D(JDIM,KDIM)
      COMMON/METRIC/
     1 XX(JDIM,KDIM),XY(JDIM,KDIM),YX(JDIM,KDIM),YY(JDIM,KDIM),
     2 C1(JDIM,KDIM),C2(JDIM,KDIM),C3(JDIM,KDIM),
     3 C4(JDIM,KDIM),C5(JDIM,KDIM),AJ(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR
      COMMON/FLWPRM/ RE, REI, DT, DTI, EPS, ERR, CONST
       DO 10 K = 2,KM
       DO 10 J = 2,JM
       CC = .5/(C1(J,K)+C3(J,K))
       PA = C1(J,K)*(P(J+1,K)+P(J-1,K))+C3(J,K)*(P(J,K+1)+P(J,K-1))
     1    +.25*C2(J,K)*(P(J+1,K+1)-P(J-1,K+1)-P(J+1,K-1)+P(J-1,K-1))
     2    +.5*C4(J,K)*(P(J+1,K)-P(J-1,K))
     3    +.5*C5(J,K)*(P(J,K+1)-P(J,K-1))
       PP = (PA+V(J,K))*CC
       ERR= ERR+(PP-P(J,K))**2
       P(J,K) = P(J,K)*(1.-CONST)+PP*CONST
   10  CONTINUE
      RETURN
      END
C********************************************************************
C     SUBROUTINE FOR SOLVING VOLTICITY TRANSPORT EQ.                *
C********************************************************************
      SUBROUTINE OMEGA
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON /VEL/ V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON /FLW/ D(JDIM,KDIM)
      COMMON/METRIC/
     1 XX(JDIM,KDIM),XY(JDIM,KDIM),YX(JDIM,KDIM),YY(JDIM,KDIM),
     2 C1(JDIM,KDIM),C2(JDIM,KDIM),C3(JDIM,KDIM),
     3 C4(JDIM,KDIM),C5(JDIM,KDIM),AJ(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR
      COMMON/FLWPRM/ RE, REI, DT, DTI, EPS, ERR, CONST
       DO 10 K = 2,KM
       DO 10 J = 2,JM
       UNL=(P(J,K+1)-P(J,K-1))*(V(J+1,K)-V(J-1,K))*.25
     2    -(P(J+1,K)-P(J-1,K))*(V(J,K+1)-V(J,K-1))*.25
       VVS=C1(J,K)*(V(J+1,K)-2.*V(J,K)+V(J-1,K))
     1    +C2(J,K)*(V(J+1,K+1)-V(J+1,K-1)-V(J-1,K+1)+V(J-1,K-1))*.25
     2    +C3(J,K)*(V(J,K+1)-2.*V(J,K)+V(J,K-1))
     3+.5*C4(J,K)*(V(J+1,K)-V(J-1,K))+.5*C5(J,K)*(V(J,K+1)-V(J,K-1))
       D(J,K)=V(J,K)+DT*(-UNL/AJ(J,K)+REI*VVS)
   10  CONTINUE
       DO 20 K = 2,KM
       DO 20 J = 2,JM
        V(J,K) = D(J,K)
   20  CONTINUE
      RETURN
      END
C********************************************************************
C     SUBROUTINE FOR BOUNDARY CONDITIONS                            *
C********************************************************************
      SUBROUTINE BC
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON/GRIDXY/X(JDIM,KDIM),Y(JDIM,KDIM)
      COMMON/METRIC/
     1 XX(JDIM,KDIM),XY(JDIM,KDIM),YX(JDIM,KDIM),YY(JDIM,KDIM),
     2 C1(JDIM,KDIM),C2(JDIM,KDIM),C3(JDIM,KDIM),
     3 C4(JDIM,KDIM),C5(JDIM,KDIM),AJ(JDIM,KDIM)
      COMMON /VEL/ V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR
      IF (ITYP.EQ.1) THEN
C***  periodic
        DO 10   K = 1,KMAX
           P(1,K) = P(JM,K)
           V(1,K) = V(JM,K)
         P(JMAX,K)= P(2,K)
         V(JMAX,K)= V(2,K)
   10    CONTINUE
      ELSE
C***  uniform flow at J=1 extapolate at K=KMAX      
        DO 20   K = 1,KMAX
           V(1,K) = 0.
        P(JMAX,K) = P(JM,K)
        V(JMAX,K) = V(JM,K)
   20   CONTINUE
      END IF
C
C      IF(IGR.NE.1) THEN
C***  symmetric
C        DO 30   K = 1,KMAX
C           P(1,K) = 0.
C           V(1,K) = 0.
C         P(JMAX,K)= 0.
C         V(JMAX,K)= 0.
C   30   CONTINUE
C      END IF
C
      IF(ITYP.EQ.1) THEN
C***  no-slip at K=1 and extrapolate  K=KMAX
        DO 40   J = 1,JMAX
          P(J,1) = 0.
          V(J,1) = -2.*C3(J,1)*(P(J,2)-P(J,1))
C         V(J,KMAX)= V(J,KM)
   40   CONTINUE
C   
       ELSE
C
C***  no-slip at K=1 and symmetric K=KMAX 
        DO 50  J = 1,JMAX
          P(J,1) = 0.
          V(J,1) = -2.*C3(J,1)*(P(J,2)-P(J,1))
        V(J,KMAX)= -2.*C3(J,KMAX)*(P(J,KM)-P(J,KMAX))*0.
   50   CONTINUE
C
      END IF
C
      RETURN
      END
C********************************************************************
C     SUBROUTINE FOR WRITING DATA INTO FILE                         * 
C********************************************************************
      SUBROUTINE OUTP
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON/GRIDXY/ X(JDIM,KDIM),Y(JDIM,KDIM)
      COMMON /VEL/ V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR
      WRITE(12,*) JMAX,KMAX
      WRITE(12,*) ((P(J,K),J=1,JMAX),K=1,KMAX),
     1            ((V(J,K),J=1,JMAX),K=1,KMAX)
      RETURN
      END
C**********************************************************************
C     SUBPOUTINE FOR OUTPUT                                           *
C**********************************************************************
      SUBROUTINE OUT2
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON/GRIDXY/ X(JDIM,KDIM),Y(JDIM,KDIM)
      COMMON /VEL/ V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR
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
            IND=INT((PJK-PMIN)/(PMAX-PMIN)*50)
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