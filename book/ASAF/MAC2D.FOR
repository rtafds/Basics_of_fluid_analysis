C**********************************************************************
C     2 dimensional flow slover (genaralized coordinate sysytem)      *
C**********************************************************************
      PROGRAM FLOW2D
      PARAMETER(JDIM=62,KDIM=31)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR, ALP
      COMMON/FLWPRM/ RE, REI, DT, DTI, EPS, ERR, CONST
C
         CALL GRID
         CALL DATA
         CALL METJ
         CALL INIT
C
      DO 10 N = 1,NSTEPS
         CALL RHS
         ERR = 0.
         DO 20 I=1,ISTEP0
            CALL PRS
            CALL PRSBC
            IF(ERR.LT.EPS) GO TO 30
   20    CONTINUE
   30    WRITE(*,*) N ,ERR
         CALL NS
         CALL BC
   10 CONTINUE
C
        WRITE(*,*) 'Save Data? (Yes=1,No=0) '
         READ(*,*) ISAVE
         IF(ISAVE.EQ.1) CALL OUTP
         CALL OUT2
      STOP
      END
C**********************************************************************
C     SUBROUTINE FOR INPUT DATA                                       *
C**********************************************************************
      SUBROUTINE DATA
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR, ALP
      COMMON/FLWPRM/ RE, REI, DT, DTI, EPS, ERR, CONST
         JM = JMAX-1
         KM = KMAX-1
         WRITE(*,*) 'INPUT NUMBER OF TIME STEP  (500)'
          READ(*,*) NSTEPS
         WRITE(*,*)'INPUT MAX. NUMBER OF ITERATION FOR POISSON EQ.(10)'
          READ(*,*) ISTEP0
         WRITE(*,*) 'INPUT REYNOLDS NUM. & TIME INCREMENT DT (100,.01)'
          READ(*,*) RE,DT
         WRITE(*,*) 'INPUT MAX. ERROR  (.0001)'
          READ(*,*) EPS
         WRITE(*,*) 'INPUT ACCERARATION PARAMETER FOR SOR METHOD (1)'
          READ(*,*) CONST
         REI = 1./RE
         DTI = 1./DT
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR GRID DATA                                        *
C**********************************************************************
      SUBROUTINE GRID
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON/GRIDXY/ X(JDIM,KDIM),Y(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR, ALP
      DIMENSION RR(KDIM)
C
      WRITE(*,*) 'DO YOU WANT TO READ GRID FILE ?  (YES=1,NO=0)  (0)'
       READ(*,*) IGR
      WRITE(*,*) 'Is the grid Periodic one ? (Yes=1,No=0)  (1)'
       READ(*,*) ITYP
C
      IF (IGR.EQ.1) THEN
C
         READ(11,*) IA,IB
         READ(11,*) ((X(I,J),I=1,IA),J=1,IB)
         READ(11,*) ((Y(I,J),I=1,IA),J=1,IB)
         JMAX = IA
         KMAX = IB
C
      ELSE
C
C*** GRID AROUND THE ELLIPSOID
         JMAX  = 62
         KMAX  = 31
         PAI   = ATAN(1.)*4.
         AA    = 0.75
         BB    = 1.0
         HH    = 0.025
         RA    = 1.15
         RR(1) = 1.0
         DO 15 K = 2,KMAX
           RR(K) = RR(K-1)+HH*RA**(K-1)
   15    CONTINUE
         WRITE(*,*) (RR(K),K=1,KMAX)
         DO 20 K = 1,KMAX
         DO 20 J = 1,JMAX
           TT = 2.*PAI*FLOAT(J-2)/FLOAT(JMAX-2)
           BC = BB+(AA-BB)*FLOAT(K-1)/FLOAT(KMAX-1)
           X(J,K) = AA*RR(K)*COS(TT)
           Y(J,K) = BC*RR(K)*SIN(TT)
   20   CONTINUE
C         WRITE(15,*) JMAX,KMAX
C         WRITE(15,*) ((X(I,J),I=1,JMAX),J=1,KMAX),
C     1               ((Y(I,J),I=1,JMAX),J=1,KMAX)
C
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
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR, ALP
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
C     SUBROUTINE FOR INITIAL CONDITION                                *
C**********************************************************************
      SUBROUTINE INIT
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON/GRIDXY/ X(JDIM,KDIM),Y(JDIM,KDIM)
      COMMON /VEL/ U(JDIM,KDIM),V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR, ALP
         WRITE(*,*) ' INPUT ATTACK ANGLE?  (30)'
          READ(*,*) ALP
         TTT = ALP*3.1415192/180.
         DO 10  K = 1,KMAX
         DO 10  J = 1,JMAX
           U(J,K) = COS(TTT)
           V(J,K) = SIN(TTT)
           P(J,K) = 0.-.5*(U(J,K)**2+V(J,K)**2)
   10    CONTINUE
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR CALCULATING RHS OF POISSON EQ.                   *
C**********************************************************************
      SUBROUTINE RHS
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON /VEL/ U(JDIM,KDIM),V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON /FLW/ Q(JDIM,KDIM),D(JDIM,KDIM)
      COMMON/METRIC/
     1 XX(JDIM,KDIM),XY(JDIM,KDIM),YX(JDIM,KDIM),YY(JDIM,KDIM),
     2 C1(JDIM,KDIM),C2(JDIM,KDIM),C3(JDIM,KDIM),
     3 C4(JDIM,KDIM),C5(JDIM,KDIM),AJ(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR, ALP
      COMMON/FLWPRM/ RE, REI, DT, DTI, EPS, ERR, CONST
      DO 10 K = 2,KM
      DO 10 J = 2,JM
        UXD = XX(J,K)*(U(J+1,K)-U(J-1,K))+YX(J,K)*(U(J,K+1)-U(J,K-1))
        UYD = XY(J,K)*(U(J+1,K)-U(J-1,K))+YY(J,K)*(U(J,K+1)-U(J,K-1))
        VXD = XX(J,K)*(V(J+1,K)-V(J-1,K))+YX(J,K)*(V(J,K+1)-V(J,K-1))
        VYD = XY(J,K)*(V(J+1,K)-V(J-1,K))+YY(J,K)*(V(J,K+1)-V(J,K-1))
        Q(J,K) = -.25*(UXD*UXD+2.*UYD*VXD+VYD*VYD)+.5*(UXD+VYD)*DTI
   10 CONTINUE
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR PRESSURE BOUNDARY CONDITION                      *
C**********************************************************************
      SUBROUTINE PRSBC
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON /VEL/ U(JDIM,KDIM),V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/METRIC/
     1 XX(JDIM,KDIM),XY(JDIM,KDIM),YX(JDIM,KDIM),YY(JDIM,KDIM),
     2 C1(JDIM,KDIM),C2(JDIM,KDIM),C3(JDIM,KDIM),
     3 C4(JDIM,KDIM),C5(JDIM,KDIM),AJ(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR, ALP
      COMMON/FLWPRM/ RE, REI, DT, DTI, EPS, ERR, CONST
C***  no-slip at K=1   free stream at K=KMAX      
      DO 10 J = 2,JM
        UL        = .5*C2(J,1)*(U(J+1,2)-U(J-1,2))+C5(J,1)*U(J,1)
        VL        = 2.*C3(J,1)*V(J,2)
        P(J,1)    = P(J,2)-REI*AJ(J,1)*(XX(J,1)*VL-XY(J,2)*UL)
        P(J,KMAX) = -.5*(U(J,KMAX)**2+V(J,KMAX)**2)
   10 CONTINUE
C   
      IF(ITYP.EQ.1) THEN
C***  periodic      
         DO 20 K = 1,KMAX
           P(1,K)    = P(JM,K)
           P(JMAX,K) = P(2,K)
   20    CONTINUE
      ELSE
C***  extapolate
         DO 30 K = 1,KMAX
           P(1,K)    = P(2,K)
           P(JMAX,K) = P(JM,K)
   30    CONTINUE
      END IF
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR CALCULATING PRESSURE                             *
C**********************************************************************
      SUBROUTINE PRS
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON /VEL/ U(JDIM,KDIM),V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON /FLW/ Q(JDIM,KDIM),D(JDIM,KDIM)
      COMMON/METRIC/
     1 XX(JDIM,KDIM),XY(JDIM,KDIM),YX(JDIM,KDIM),YY(JDIM,KDIM),
     2 C1(JDIM,KDIM),C2(JDIM,KDIM),C3(JDIM,KDIM),
     3 C4(JDIM,KDIM),C5(JDIM,KDIM),AJ(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR, ALP
      COMMON/FLWPRM/ RE, REI, DT, DTI, EPS, ERR, CONST
      DO 10 K = 2,KM
      DO 10 J = 2,JM
        CC = .5/(C1(J,K)+C3(J,K))
        PA = C1(J,K)*(P(J+1,K)+P(J-1,K))+C3(J,K)*(P(J,K+1)+P(J,K-1))
     1      +.25*C2(J,K)*(P(J+1,K+1)-P(J-1,K+1)-P(J+1,K-1)+P(J-1,K-1))
     2      +.5*C4(J,K)*(P(J+1,K)-P(J-1,K))
     3      +.5*C5(J,K)*(P(J,K+1)-P(J,K-1))
        PP = (PA-Q(J,K))*CC
        ERR = ERR+(PP-P(J,K))**2
        P(J,K) = P(J,K)*(1.-CONST)+PP*CONST
   10 CONTINUE
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR CALCULATING VELOCITY                             *
C**********************************************************************
      SUBROUTINE NS
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON /VEL/ U(JDIM,KDIM),V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON /FLW/ Q(JDIM,KDIM),D(JDIM,KDIM)
      COMMON/METRIC/
     1 XX(JDIM,KDIM),XY(JDIM,KDIM),YX(JDIM,KDIM),YY(JDIM,KDIM),
     2 C1(JDIM,KDIM),C2(JDIM,KDIM),C3(JDIM,KDIM),
     3 C4(JDIM,KDIM),C5(JDIM,KDIM),AJ(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR, ALP
      COMMON/FLWPRM/ RE, REI, DT, DTI, EPS, ERR, CONST
C
      DO 10 K = 2,KM
      DO 10 J = 2,JM
C***  conservative form
        UNL = XX(J,K)*(U(J+1,K)**2-U(J-1,K)**2)*.5
     1       +YX(J,K)*(U(J,K+1)**2-U(J,K-1)**2)*.5
     2       +XY(J,K)*(U(J+1,K)*V(J+1,K)-U(J-1,K)*V(J-1,K))*.5
     3       +YY(J,K)*(U(J,K+1)*V(J,K+1)-U(J,K-1)*V(J,K-1))*.5
        VNL = XX(J,K)*(U(J+1,K)*V(J+1,K)-U(J-1,K)*V(J-1,K))*.5
     1       +YX(J,K)*(U(J,K+1)*V(J,K+1)-U(J,K-1)*V(J,K-1))*.5
     2       +XY(J,K)*(V(J+1,K)**2-V(J-1,K)**2)*.5
     3       +YY(J,K)*(V(J,K+1)**2-V(J,K-1)**2)*.5
C***  non-conservative form
C       UA  = (U(J,K)*XX(J,K)+V(J,K)*XY(J,K))*.5
C       VA  = (U(J,K)*YX(J,K)+V(J,K)*YY(J,K))*.5
C       UNL = UA*(U(J+1,K)-U(J-1,K))+VA*(U(J,K+1)-U(J,K-1))
C       VNL = UA*(V(J+1,K)-V(J-1,K))+VA*(V(J,K+1)-V(J,K-1))
        UVS = C1(J,K)*(U(J+1,K)-2.*U(J,K)+U(J-1,K))
     1       +C2(J,K)*(U(J+1,K+1)-U(J+1,K-1)-U(J-1,K+1)+U(J-1,K-1))*.25
     2       +C3(J,K)*(U(J,K+1)-2.*U(J,K)+U(J,K-1))+.5*C4(J,K)
     3       *(U(J+1,K)-U(J-1,K))+.5*C5(J,K)*(U(J,K+1)-U(J,K-1))
        VVS = C1(J,K)*(V(J+1,K)-2.*V(J,K)+V(J-1,K))
     1       +C2(J,K)*(V(J+1,K+1)-V(J+1,K-1)-V(J-1,K+1)+V(J-1,K-1))*.25
     2       +C3(J,K)*(V(J,K+1)-2.*V(J,K)+V(J,K-1))+.5*C4(J,K)
     3       *(V(J+1,K)-V(J-1,K))+.5*C5(J,K)*(V(J,K+1)-V(J,K-1))
C
        D(J,K) = U(J,K)+DT*(-UNL-.5*XX(J,K)*(P(J+1,K)-P(J-1,K))
     1          -.5*YX(J,K)*(P(J,K+1)-P(J,K-1))+REI*UVS)
        Q(J,K) = V(J,K)+DT*(-VNL-.5*XY(J,K)*(P(J+1,K)-P(J-1,K))
     1          -.5*YY(J,K)*(P(J,K+1)-P(J,K-1))+REI*VVS)
   10 CONTINUE
      DO 20 K = 2,KM
      DO 20 J = 2,JM
        U(J,K) = D(J,K)
        V(J,K) = Q(J,K)
   20 CONTINUE
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR BOUNDARY CONDITION FOR VELOCITY                  *
C**********************************************************************
      SUBROUTINE BC
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON /VEL/ U(JDIM,KDIM),V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR, ALP
C
C***  no-slip at K=1  free stream at K=KMAX
      TTT=ALP*3.141592/180.
      DO 10 J = 1,JMAX
         U(J,1)    = 0.
         V(J,1)    = 0.
         U(J,KMAX) = COS(TTT)
         V(J,KMAX) = SIN(TTT)
   10 CONTINUE
C
      IF(ITYP.EQ.1) THEN
C***  poriodic
          DO 20 K = 1,KMAX
            U(1,K)    = U(JM,K)
            V(1,K)    = V(JM,K)
            U(JMAX,K) = U(2,K)
            V(JMAX,K) = V(2,K)
   20     CONTINUE
      ELSE
C*** symmetric
          DO 30 K = 1,KMAX
            U(1,K)    = U(2,K)
            V(1,K)    = 0.
            U(JMAX,K) = U(JM,K)
            V(JMAX,K) = 0.
   30     CONTINUE
      END IF
      RETURN
      END
C**********************************************************************
C     WRITE DATA INTO FILE                                            *
C**********************************************************************
      SUBROUTINE OUTP
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON /VEL/ U(JDIM,KDIM),V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR, ALP
         WRITE(12,*) JMAX,KMAX
         WRITE(12,*) ((U(J,K),J=1,JMAX),K=1,KMAX),
     1               ((V(J,K),J=1,JMAX),K=1,KMAX),
     1               ((P(J,K),J=1,JMAX),K=1,KMAX)
      RETURN
      END
C**********************************************************************
C     SUBPOUTINE FOR OUTPUT                                           *
C**********************************************************************
      SUBROUTINE OUT2
      PARAMETER ( JDIM=62,KDIM=31)
      COMMON/GRIDXY/ X(JDIM,KDIM),Y(JDIM,KDIM)
      COMMON /VEL/ U(JDIM,KDIM),V(JDIM,KDIM),P(JDIM,KDIM)
      COMMON/METRIC/
     1 XX(JDIM,KDIM),XY(JDIM,KDIM),YX(JDIM,KDIM),YY(JDIM,KDIM),
     2 C1(JDIM,KDIM),C2(JDIM,KDIM),C3(JDIM,KDIM),
     3 C4(JDIM,KDIM),C5(JDIM,KDIM),AJ(JDIM,KDIM)
      COMMON/GRDPRM/ JM, JMAX, KM, KMAX, ISTEP0, NSTEPS, ITYP, IGR, ALP
      COMMON/FLWPRM/ RE, REI, DT, DTI, EPS, ERR, CONST
      CHARACTER*2 Z(39,31)
      DIMENSION INDEX(39,31)
C
      AREA(P1,Q1,P2,Q2,P3,Q3)=ABS((P2-P1)*(Q3-Q1)-(P3-P1)*(Q2-Q1))
C
      DO 20 J=1,JMAX
      P(J,1)=0.
      DO 25 K=2,KMAX
      P(J,K)=P(J,K-1)+.5*(U(J,K)+U(J,K-1))*(Y(J,K)-Y(J,K-1))
     1               -.5*(V(J,K)+V(J,K-1))*(X(J,K)-X(J,K-1))
   25 CONTINUE
   20 CONTINUE
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
      DO 40 K2=1,25
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
            IND=INT((PJK-PMIN)/(PMAX-PMIN)*80)
            INDEX(J2,K2)=(MOD(IND,9)+1)*11
            GO TO 40
          END IF
   60   CONTINUE
   40 CONTINUE
      DO 70 K=1,25
      DO 70 J=1,39
        IND=INDEX(J,K)
        IF(IND.EQ.0) THEN
          Z(J,K)='//'
        ELSE
          IF(MOD(IND,2).EQ.0) Z(J,K)='..'
          IF(MOD(IND,2).EQ.1) Z(J,K)='**'
        END IF
   70 CONTINUE
      DO 50 K=25,1,-1
        WRITE(*,600) (INDEX(J,K),J=1,39)
   50 CONTINUE
  600   FORMAT(1H ,39I2)
      WRITE(*,*) 'INPUT ANY NUMBER (1) '
       READ(*,*) AA
      WRITE(*,*) '  '
      DO 75 K=25,1,-1
        WRITE(*,601) (Z(J,K),J=1,39)
   75 CONTINUE
  601   FORMAT(1H ,39A2)
      WRITE(*,*) 'PMAX=',PMAX,'   PMIN=',PMIN
      RETURN
      END
