C**********************************************************************
C     THREE DIMENSIONAL FLOW SOLVER
C**********************************************************************
      PROGRAM FLOW3D
      PARAMETER (JD=22,KD=16,LD=11)
      COMMON/XYZ/X(JD,KD,LD),Y(JD,KD,LD),Z(JD,KD,LD)
      COMMON/MET1/XX(JD,KD,LD),XY(JD,KD,LD),XZ(JD,KD,LD)
     1           ,YX(JD,KD,LD),YY(JD,KD,LD),YZ(JD,KD,LD)
     2           ,ZX(JD,KD,LD),ZY(JD,KD,LD),ZZ(JD,KD,LD)
      COMMON/MET2/C1(JD,KD,LD),C2(JD,KD,LD),C3(JD,KD,LD)
     1           ,C4(JD,KD,LD),C5(JD,KD,LD),C6(JD,KD,LD)
     2           ,C7(JD,KD,LD),C8(JD,KD,LD),C9(JD,KD,LD)
      COMMON/FLOW1/U(JD,KD,LD),V(JD,KD,LD),W(JD,KD,LD),P(JD,KD,LD)
      COMMON/FLOW2/D(JD,KD,LD),Q(JD,KD,LD),S(JD,KD,LD)
      COMMON/FV1/JM,KM,LM,JM1,KM1,LM1,IL,NL,IW,NC
      COMMON/FV2/DT,DTI,REI,EPS,CONST
C
         CALL EDATA
         CALL GRID
         CALL METJ
         CALL METJ2
         CALL INIT
         CALL BC
C
         DO 10 NC=1,NL
           CALL RHS
           CALL PRS
           CALL NS
           CALL BC
             IF(NC.GT.20.AND.NC-NC/IW*IW.EQ.0) THEN
               CALL OUT
             END IF
   10    CONTINUE
         IWR = 7
           IF(IWR.GE.7) THEN
               WRITE(IWR,*) JD,KD
               DO 15 K=1,KD
               DO 15 J=1,JD
                  WRITE(IWR,*) X(J,K,6),Y(J,K,6),U(J,K,6),V(J,K,6)
   15          CONTINUE
           END IF
         WRITE(*,*) 'DO YOU WANT TO SAVE DATA? '
          READ(*,*) ISV
           IF(ISV.EQ.1) CALL SDATA
         DO 20 LL=1,LD
           WRITE(*,*) 'L  6    (0<L<11)  (=0 STOP)'
           READ(*,*) L
           IF(L.LE.0) STOP
           CALL OUT2(L)
   20    CONTINUE
        STOP
      END
C**********************************************************************
C     SUBROUTINE FOR BOUNDARY CONDITION
C**********************************************************************
      SUBROUTINE BC
      PARAMETER (JD=22,KD=16,LD=11)
      COMMON/FLOW1/U(JD,KD,LD),V(JD,KD,LD),W(JD,KD,LD),P(JD,KD,LD)
      COMMON/FV1/JM,KM,LM,JM1,KM1,LM1,IL,NL,IW,NC
      COMMON/IBND/LFL
C***  X-Y  LFL=0 :NO-SLIP   =1 :FREE-SLIP   =2 :PERIODIC
        IF(LFL.EQ.1) THEN
          DO 50 K = 1,KM
          DO 50 J = 1,JM
            U(J,K,1)  =  U(J,K,2)
            V(J,K,1)  =  V(J,K,2)
            W(J,K,1)  =  W(J,K,2)*0.
            U(J,K,LM) =  U(J,K,LM1)
            V(J,K,LM) =  V(J,K,LM1)
            W(J,K,LM) =  W(J,K,LM1)*0.
   50     CONTINUE
        END IF
        IF(LFL.EQ.0) THEN
          DO 55 K = 1,KM
          DO 55 J = 1,JM
            U(J,K,1)  =  U(J,K,2)*0.
            V(J,K,1)  = -V(J,K,2)*0.
            W(J,K,1)  =  W(J,K,2)*0.
            U(J,K,LM) =  U(J,K,LM1)*0.
            V(J,K,LM) =  V(J,K,LM1)*0.
            W(J,K,LM) =  W(J,K,LM1)*0.
   55     CONTINUE
        END IF
        IF(LFL.EQ.2) THEN
          DO 59 K = 1,KM
          DO 59 J = 1,JM
            U(J,K,1)  =  U(J,K,LM-1)
            V(J,K,1)  =  V(J,K,LM-1)
            W(J,K,1)  =  W(J,K,LM-1)
            U(J,K,LM) =  U(J,K,2)
            V(J,K,LM) =  V(J,K,2)
            W(J,K,LM) =  W(J,K,2)
   59     CONTINUE
        END IF
C*** X-Z
          DO 51 L = 1,LM
          DO 51 J = 1,JM
            U(J,1,L)  =  U(J,3,L)*0.
            V(J,1,L)  = -V(J,3,L)*0.
            W(J,1,L)  =  W(J,3,L)*0.
            U(J,KM,L) =  U(J,KM-2,L)*0.-1.
            V(J,KM,L) = -V(J,KM-2,L)*0.
            W(J,KM,L) =  W(J,KM-2,L)*0.
   51     CONTINUE
C*** X-Y (Periodic)
          DO 52 L = 1,LM
          DO 52 K = 1,KM
            U(1,K,L)  =  U(JM-1,K,L)
            V(1,K,L)  =  V(JM-1,K,L)
            W(1,K,L)  =  W(JM-1,K,L)
            U(JM,K,L) =  U(2,K,L)
            V(JM,K,L) =  V(2,K,L)
            W(JM,K,L) =  W(2,K,L)
   52     CONTINUE
        RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR PREPARING DATA
C**********************************************************************
      SUBROUTINE EDATA
      COMMON/FV1/JM,KM,LM,JM1,KM1,LM1,IL,NL,IW,NC
      COMMON/FV2/DT,DTI,REI,EPS,CONST
      COMMON/IBND/LFL
      WRITE(*,*) 'NUMBER OF GRID (X,Y,Z DIRECTION)? (21,15,10)'
       READ(*,*) JM1,KM1,LM1
          JM  =  JM1 + 1
          KM  =  KM1 + 1
          LM  =  LM1 + 1
          IL  =  10
      WRITE(*,*) 'TIME INCREMENT DT & REYNOLDS NUMBER RE? (.01,100)'
       READ(*,*) DT,RE
      WRITE(*,*) 'NUMBER OF TIME STEPS? (400)' 
       READ(*,*) NL
      WRITE(*,*) 'BOUNDARY TYPE? 0:NO-SLIP 1:FREE-SLIP 2:PERI ?? (1)'
       READ(*,*)LFL
          IW    =  NL+10
          DTI   =  1./DT
          REI   =  1./RE
          EPS   =  0.001
          CONST = 1.
        RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR MAKING GRIDS
C**********************************************************************
      SUBROUTINE GRID
      PARAMETER (JD=22,KD=16,LD=11)
      COMMON/XYZ/X(JD,KD,LD),Y(JD,KD,LD),Z(JD,KD,LD)
      COMMON/FV1/JM,KM,LM,JM1,KM1,LM1,IL,NL,IW,NC
      COMMON/HGT/QS(JD,KD,2),ZS(JD,KD),ZQ(2*LD+1)
      DIMENSION HT(KD),RA(JD)
          ITYP  =  1
          IF(ITYP.EQ.1) THEN
            WRITE(*,*) ' R=1.+r**1+r**2+r**3+... r? (1.1)'
             READ(*,*) R
            WRITE(*,*) 'RADIUS OF OUTER CIRCLE ? (5.0) '
             READ(*,*) W
            WRITE(*,*) 'ANGLE BETWEEN UNIT NORMAL AND AXIS (30.)'
             READ(*,*) TET
              TET  =  90. - TET
            WRITE(*,*) 'HIGHT OF THE CYLINDER? (4.0)'
              READ(*,*) H
              PAI  =  ATAN(1.0)*4.0
              TE   =  TET/180.*PAI
              A    =  (W-1.)*(R-1.)/(R**FLOAT(KM-1)-1.)
              RA(1)=  1.
            DO 220 K = 2,KM
  220         RA(K) = RA(K-1)+A*R**(K-2)
            WRITE(*,*) 'UNIFORM MESH (1) OR NOT (FACTOR) (0,0.9)'
             READ(*,*) IGTY,BB
              MZN  = 2*LM-1
              BD   = 0.5*ALOG((1.0+BB)/(1.0-BB))
            DO 170 K = 1,MZN
              EY   = (-1.+2.*FLOAT(K-1)/FLOAT(MZN-1))*2.*BD
              ZQ(K)= (1.+(EXP(EY)-1.)/(EXP(EY)+1.)
     1              *(EXP(2.*BD)+1.)/(EXP(2.*BD)-1.))*.5
  170       CONTINUE
            WRITE(*,*) (ZQ(I),I=1,LM)
            DO 230 L = 1,LM
            DO 230 K = 1,KM
            DO 230 J = 1,JM
              A  = .25*RA(K)/SIN(TE)
              B  = .25*RA(K)
              AL = PAI*2.*FLOAT(J-2)/FLOAT(JD-2)
              X(J,K,L)  = A*COS(AL)
              Y(J,K,L)  = B*SIN(AL)
              Z(J,K,L)  = H*ZQ(L)
              IF(IGTY.EQ.1) Z(J,K,L)  = H*FLOAT(L-1)/FLOAT(LM-1)
              IF(TET.NE.90.)
     1        X(J,K,L)  = A*COS(AL)+Z(J,K,L)/TAN(TE)
  230       CONTINUE
          ELSE
            READ(11) JB,KB,LB
            READ(11) (((X(J,K,L),J=1,JB),K=1,KB),L=1,LB)
            READ(11) (((Y(J,K,L),J=1,JB),K=1,KB),L=1,LB)
            READ(11) (((Z(J,K,L),J=1,JB),K=1,KB),L=1,LB)
          END IF
        RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR INITIAL CONDITION
C**********************************************************************
      SUBROUTINE INIT
      PARAMETER (JD=22,KD=16,LD=11)
      COMMON/XYZ/X(JD,KD,LD),Y(JD,KD,LD),Z(JD,KD,LD)
      COMMON/MET1/XX(JD,KD,LD),XY(JD,KD,LD),XZ(JD,KD,LD)
     1           ,YX(JD,KD,LD),YY(JD,KD,LD),YZ(JD,KD,LD)
     2           ,ZX(JD,KD,LD),ZY(JD,KD,LD),ZZ(JD,KD,LD)
      COMMON/FLOW1/U(JD,KD,LD),V(JD,KD,LD),W(JD,KD,LD),P(JD,KD,LD)
      COMMON/FV1/JM,KM,LM,JM1,KM1,LM1,IL,NL,IW,NC
      DO 70 L = 1,LM
      DO 70 K = 1,KM
      DO 70 J = 1,JM
      AJ=XX(J,K,L)*YY(J,K,L)*ZZ(J,K,L)+XY(J,K,L)*YZ(J,K,L)*ZX(J,K,L)
     1  +XZ(J,K,L)*YX(J,K,L)*ZY(J,K,L)-XZ(J,K,L)*YY(J,K,L)*ZX(J,K,L)
     2  -XY(J,K,L)*YX(J,K,L)*ZZ(J,K,L)-XX(J,K,L)*YZ(J,K,L)*ZY(J,K,L)
      U(J,K,L)=(YY(J,K,L)*ZZ(J,K,L)-YZ(J,K,L)*ZY(J,K,L))/AJ*0.- 1.
     1        *(Z(J,K,L)/Z(J,K,LM))**(1./5.)
      V(J,K,L)=(YZ(J,K,L)*ZX(J,K,L)-YX(J,K,L)*ZZ(J,K,L))/AJ*0.
      W(J,K,L)=(YX(J,K,L)*ZY(J,K,L)-ZX(J,K,L)*YY(J,K,L))/AJ*0.
   70 CONTINUE
      DO 71 L = 1,LM
      DO 71 K = 1,KM
      DO 71 J = 1,JM
      P(J,K,L)=-.5*(U(J,K,L)**2+V(J,K,L)**2+W(J,K,L)**2)*0.
   71 CONTINUE
      PMIN=P(1,1,1)
      DO 72 L = 1,LM
      DO 72 K = 1,KM
      DO 72 J=  1,JM
      IF(P(J,K,L).LT.PMIN) PMIN=P(J,K,L)
   72 CONTINUE
      DO 74 L = 1,LM
      DO 74 K = 1,KM
      DO 74 J = 1,JM
      P(J,K,L)=P(J,K,L)-PMIN
   74 CONTINUE
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR CALCULATING METRICS 
C**********************************************************************
      SUBROUTINE METJ
      PARAMETER (JD=22,KD=16,LD=11)
      COMMON/XYZ/X(JD,KD,LD),Y(JD,KD,LD),Z(JD,KD,LD)
      COMMON/MET1/XX(JD,KD,LD),XY(JD,KD,LD),XZ(JD,KD,LD)
     1           ,YX(JD,KD,LD),YY(JD,KD,LD),YZ(JD,KD,LD)
     2           ,ZX(JD,KD,LD),ZY(JD,KD,LD),ZZ(JD,KD,LD)
      COMMON/FLOW2/D(JD,KD,LD),Q(JD,KD,LD),S(JD,KD,LD)
      COMMON/FV1/JM,KM,LM,JM1,KM1,LM1,IL,NL,IW,NC
      DIMENSION XJ(JD),YJ(JD),ZJ(JD),XK(JD),YK(JD),ZK(JD)
     1         ,XL(JD),YL(JD),ZL(JD)
      EEP=1.E-30
      DX2        = 2.0
      DY2        = 2.0
      DZ2        = 2.0
C
      DO 51 L = 1,LM
      DO 51 K = 1,KM
      J          = 1
      J1         = J + 1
      J2         = J + 2
      XJ(J)      = -(3.0*X(J,K,L) -4.0*X(J1,K,L) +X(J2,K,L))/DX2
      YJ(J)      = -(3.0*Y(J,K,L) -4.0*Y(J1,K,L) +Y(J2,K,L))/DX2
      ZJ(J)      = -(3.0*Z(J,K,L) -4.0*Z(J1,K,L) +Z(J2,K,L))/DX2
C
      DO 11 J = 2,JM1
      XJ(J)      = (X(J+1,K,L) -X(J-1,K,L))/DX2
      YJ(J)      = (Y(J+1,K,L) -Y(J-1,K,L))/DX2
      ZJ(J)      = (Z(J+1,K,L) -Z(J-1,K,L))/DX2
   11 CONTINUE
C
      J          = JM
      J1         = J - 1
      J2         = J - 2
      XJ(J)      =  (3.0*X(J,K,L) -4.0*X(J1,K,L) +X(J2,K,L))/DX2
      YJ(J)      =  (3.0*Y(J,K,L) -4.0*Y(J1,K,L) +Y(J2,K,L))/DX2
      ZJ(J)      =  (3.0*Z(J,K,L) -4.0*Z(J1,K,L) +Z(J2,K,L))/DX2
C
      IF( K.EQ. 1            ) THEN
          K1         = K + 1
          K2         = K + 2
          DO 21 J = 1,JM
          XK(J)      = -(3.0*X(J,K,L) -4.0*X(J,K1,L) +X(J,K2,L))/DY2
          YK(J)      = -(3.0*Y(J,K,L) -4.0*Y(J,K1,L) +Y(J,K2,L))/DY2
          ZK(J)      = -(3.0*Z(J,K,L) -4.0*Z(J,K1,L) +Z(J,K2,L))/DY2
   21     CONTINUE
C
      ELSE IF( K.EQ.KM          ) THEN
          K1         = K - 1
          K2         = K - 2
          DO 23 J = 1,JM
          XK(J)      =  (3.0*X(J,K,L) -4.0*X(J,K1,L) +X(J,K2,L))/DY2
          YK(J)      =  (3.0*Y(J,K,L) -4.0*Y(J,K1,L) +Y(J,K2,L))/DY2
          ZK(J)      =  (3.0*Z(J,K,L) -4.0*Z(J,K1,L) +Z(J,K2,L))/DY2
   23     CONTINUE
C
      ELSE
          KP = K+1
          KR = K-1
          DO 22 J = 1,JM
          XK(J)      = (X(J,KP,L) -X(J,KR,L))/DY2
          YK(J)      = (Y(J,KP,L) -Y(J,KR,L))/DY2
          ZK(J)      = (Z(J,KP,L) -Z(J,KR,L))/DY2
   22     CONTINUE
C
      END IF
C
      IF( L.EQ.1 ) THEN
          L1         = L + 1
          L2         = L + 2
          DO 31 J = 1,JM
          XL(J)      = -(3.0*X(J,K,L) -4.0*X(J,K,L1) +X(J,K,L2))/DZ2
          YL(J)      = -(3.0*Y(J,K,L) -4.0*Y(J,K,L1) +Y(J,K,L2))/DZ2
          ZL(J)      = -(3.0*Z(J,K,L) -4.0*Z(J,K,L1) +Z(J,K,L2))/DZ2
   31     CONTINUE
C
      ELSE IF( L.EQ.LM   ) THEN
          L1         = L - 1
          L2         = L - 2
          DO 33 J = 1,JM
          XL(J)      =  (3.0*X(J,K,L) -4.0*X(J,K,L1) +X(J,K,L2))/DZ2
          YL(J)      =  (3.0*Y(J,K,L) -4.0*Y(J,K,L1) +Y(J,K,L2))/DZ2
          ZL(J)      =  (3.0*Z(J,K,L) -4.0*Z(J,K,L1) +Z(J,K,L2))/DZ2
   33     CONTINUE
C
      ELSE
          DO 32 J = 1,JM
          XL(J)      = (X(J,K,L+1) -X(J,K,L-1))/DZ2
          YL(J)      = (Y(J,K,L+1) -Y(J,K,L-1))/DZ2
          ZL(J)      = (Z(J,K,L+1) -Z(J,K,L-1))/DZ2
   32     CONTINUE
      END IF
C
      DO 41 J = 1,JM
      Q(J,K,L)   = XJ(J)*YK(J)*ZL(J) +XL(J)*YJ(J)*ZK(J) +
     1             XK(J)*YL(J)*ZJ(J) -XJ(J)*YL(J)*ZK(J) -
     2             XL(J)*YK(J)*ZJ(J) -XK(J)*YJ(J)*ZL(J)
   41 CONTINUE
C
      DO 42 J = 1,JM
      Q(J,K,L)   = 1.0/( Q(J,K,L)   +EEP )
   42 CONTINUE
C
      DO 45 J = 1,JM
      XX(J,K,L)  = (YK(J)*ZL(J) -ZK(J)*YL(J))*Q(J,K,L)
      XY(J,K,L)  = (ZK(J)*XL(J) -XK(J)*ZL(J))*Q(J,K,L)
      XZ(J,K,L)  = (XK(J)*YL(J) -YK(J)*XL(J))*Q(J,K,L)
      YX(J,K,L)  = (ZJ(J)*YL(J) -YJ(J)*ZL(J))*Q(J,K,L)
      YY(J,K,L)  = (XJ(J)*ZL(J) -ZJ(J)*XL(J))*Q(J,K,L)
      YZ(J,K,L)  = (YJ(J)*XL(J) -XJ(J)*YL(J))*Q(J,K,L)
      ZX(J,K,L)  = (YJ(J)*ZK(J) -ZJ(J)*YK(J))*Q(J,K,L)
      ZY(J,K,L)  = (XK(J)*ZJ(J) -XJ(J)*ZK(J))*Q(J,K,L)
      ZZ(J,K,L)  = (XJ(J)*YK(J) -YJ(J)*XK(J))*Q(J,K,L)
   45 CONTINUE
C
   51 CONTINUE
C
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR CALCULATING METRICS (FOR LAPLACIAN)
C**********************************************************************
      SUBROUTINE METJ2
      PARAMETER (JD=22,KD=16,LD=11)
      COMMON/MET1/XX(JD,KD,LD),XY(JD,KD,LD),XZ(JD,KD,LD)
     1           ,YX(JD,KD,LD),YY(JD,KD,LD),YZ(JD,KD,LD)
     2           ,ZX(JD,KD,LD),ZY(JD,KD,LD),ZZ(JD,KD,LD)
      COMMON/MET2/C1(JD,KD,LD),C2(JD,KD,LD),C3(JD,KD,LD)
     1           ,C4(JD,KD,LD),C5(JD,KD,LD),C6(JD,KD,LD)
     2           ,C7(JD,KD,LD),C8(JD,KD,LD),C9(JD,KD,LD)
      COMMON/FV1/JM,KM,LM,JM1,KM1,LM1,IL,NL,IW,NC
      DO 10 L = 1,LM
      DO 10 K = 1,KM
      DO 10 J = 1,JM
      C1(J,K,L)=XX(J,K,L)**2+XY(J,K,L)**2+XZ(J,K,L)**2
      C2(J,K,L)=YX(J,K,L)**2+YY(J,K,L)**2+YZ(J,K,L)**2
      C3(J,K,L)=ZX(J,K,L)**2+ZY(J,K,L)**2+ZZ(J,K,L)**2
      C4(J,K,L)=XX(J,K,L)*YX(J,K,L)+XY(J,K,L)*YY(J,K,L)
     1         +XZ(J,K,L)*YZ(J,K,L)
      C5(J,K,L)=YX(J,K,L)*ZX(J,K,L)+YY(J,K,L)*ZY(J,K,L)
     1         +YZ(J,K,L)*ZZ(J,K,L)
      C6(J,K,L)=ZX(J,K,L)*XX(J,K,L)+ZY(J,K,L)*XY(J,K,L)
     1         +ZZ(J,K,L)*XZ(J,K,L)
   10 CONTINUE
      DO 11 L = 2,LM1
      DO 11 K = 2,KM1
      DO 11 J = 2,JM1
      C77=XX(J,K,L)*(XX(J+1,K,L)-XX(J-1,K,L))+YX(J,K,L)*(XX(J,K+1,L)
     1   -XX(J,K-1,L))+ZX(J,K,L)*(XX(J,K,L+1)-XX(J,K,L-1))
     2   +XY(J,K,L)*(XY(J+1,K,L)-XY(J-1,K,L))+YY(J,K,L)*(XY(J,K+1,L)
     3   -XY(J,K-1,L))+ZY(J,K,L)*(XY(J,K,L+1)-XY(J,K,L-1))
     4   +XZ(J,K,L)*(XZ(J+1,K,L)-XZ(J-1,K,L))+YZ(J,K,L)*(XZ(J,K+1,L)
     5   -XZ(J,K-1,L))+ZZ(J,K,L)*(XZ(J,K,L+1)-XZ(J,K,L-1))
      C88=XX(J,K,L)*(YX(J+1,K,L)-YX(J-1,K,L))+YX(J,K,L)*(YX(J,K+1,L)
     1   -YX(J,K-1,L))+ZX(J,K,L)*(YX(J,K,L+1)-YX(J,K,L-1))
     2   +XY(J,K,L)*(YY(J+1,K,L)-YY(J-1,K,L))+YY(J,K,L)*(YY(J,K+1,L)
     3   -YY(J,K-1,L))+ZY(J,K,L)*(YY(J,K,L+1)-YY(J,K,L-1))
     4   +XZ(J,K,L)*(YZ(J+1,K,L)-YZ(J-1,K,L))+YZ(J,K,L)*(YZ(J,K+1,L)
     5   -YZ(J,K-1,L))+ZZ(J,K,L)*(YZ(J,K,L+1)-YZ(J,K,L-1))
      C99=XX(J,K,L)*(ZX(J+1,K,L)-ZX(J-1,K,L))+YX(J,K,L)*(ZX(J,K+1,L)
     1   -ZX(J,K-1,L))+ZX(J,K,L)*(ZX(J,K,L+1)-ZX(J,K,L-1))
     2   +XY(J,K,L)*(ZY(J+1,K,L)-ZY(J-1,K,L))+YY(J,K,L)*(ZY(J,K+1,L)
     3   -ZY(J,K-1,L))+ZY(J,K,L)*(ZY(J,K,L+1)-ZY(J,K,L-1))
     4   +XZ(J,K,L)*(ZZ(J+1,K,L)-ZZ(J-1,K,L))+YZ(J,K,L)*(ZZ(J,K+1,L)
     5   -ZZ(J,K-1,L))+ZZ(J,K,L)*(ZZ(J,K,L+1)-ZZ(J,K,L-1))
      C7(J,K,L) = C77*.5
      C8(J,K,L) = C88*.5
      C9(J,K,L) = C99*.5
   11 CONTINUE
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR GETTING NEW VELOCITY
C**********************************************************************
      SUBROUTINE NS
      PARAMETER (JD=22,KD=16,LD=11)
      COMMON/MET1/XX(JD,KD,LD),XY(JD,KD,LD),XZ(JD,KD,LD)
     1           ,YX(JD,KD,LD),YY(JD,KD,LD),YZ(JD,KD,LD)
     2           ,ZX(JD,KD,LD),ZY(JD,KD,LD),ZZ(JD,KD,LD)
      COMMON/MET2/C1(JD,KD,LD),C2(JD,KD,LD),C3(JD,KD,LD)
     1           ,C4(JD,KD,LD),C5(JD,KD,LD),C6(JD,KD,LD)
     2           ,C7(JD,KD,LD),C8(JD,KD,LD),C9(JD,KD,LD)
      COMMON/FLOW1/U(JD,KD,LD),V(JD,KD,LD),W(JD,KD,LD),P(JD,KD,LD)
      COMMON/FLOW2/D(JD,KD,LD),Q(JD,KD,LD),S(JD,KD,LD)
      COMMON/FV1/JM,KM,LM,JM1,KM1,LM1,IL,NL,IW,NC
      COMMON/FV2/DT,DTI,REI,EPS,CONST
      DIMENSION US(JD),UA(JD),UN(JD),UL(JD),VS(JD),VA(JD)
     1         ,VN(JD),VL(JD),WS(JD),WA(JD),WN(JD),WL(JD)
      DO 40 L = 2,LM1
      DO 40 K = 2,KM1
      DO 41 J = 2,JM1
        UP = U(J,K,L)*XX(J,K,L)+V(J,K,L)*XY(J,K,L)+W(J,K,L)*XZ(J,K,L)
        VP = U(J,K,L)*YX(J,K,L)+V(J,K,L)*YY(J,K,L)+W(J,K,L)*YZ(J,K,L)
        WP = U(J,K,L)*ZX(J,K,L)+V(J,K,L)*ZY(J,K,L)+W(J,K,L)*ZZ(J,K,L)
        US(J) = UP
        VS(J) = VP
        WS(J) = WP
        UA(J) = ABS(UP)
        VA(J) = ABS(VP)
        WA(J) = ABS(WP)
   41 CONTINUE
      DO 43 J = 2,JM1
      UN(J)=US(J)*(U(J+1,K,L)-U(J-1,K,L))
     1     -UA(J)*(U(J+1,K,L)-2.*U(J,K,L)+U(J-1,K,L))
     3     +VS(J)*(U(J,K+1,L)-U(J,K-1,L))
     4     -VA(J)*(U(J,K+1,L)-2.*U(J,K,L)+U(J,K-1,L))
     6     +WS(J)*(U(J,K,L+1)-U(J,K,L-1))
     7     -WA(J)*(U(J,K,L+1)-2.*U(J,K,L)+U(J,K,L-1))
      VN(J)=US(J)*(V(J+1,K,L)-V(J-1,K,L))
     1     -UA(J)*(V(J+1,K,L)-2.*V(J,K,L)+V(J-1,K,L))
     3     +VS(J)*(V(J,K+1,L)-V(J,K-1,L))
     4     -VA(J)*(V(J,K+1,L)-2.*V(J,K,L)+V(J,K-1,L))
     6     +WS(J)*(V(J,K,L+1)-V(J,K,L-1))
     7     -WA(J)*(V(J,K,L+1)-2.*V(J,K,L)+V(J,K,L-1))
      WN(J)=US(J)*(W(J+1,K,L)-W(J-1,K,L))
     1     -UA(J)*(W(J+1,K,L)-2.*W(J,K,L)+W(J-1,K,L))
     3     +VS(J)*(W(J,K+1,L)-W(J,K-1,L))
     4     -VA(J)*(W(J,K+1,L)-2.*W(J,K,L)+W(J,K-1,L))
     6     +WS(J)*(W(J,K,L+1)-W(J,K,L-1))
     7     -WA(J)*(W(J,K,L+1)-2.*W(J,K,L)+W(J,K,L-1))
   43 CONTINUE
      DO 44 J = 2,JM1
      UL(J)=C1(J,K,L)*(U(J+1,K,L)-2.*U(J,K,L)+U(J-1,K,L))
     1     +C2(J,K,L)*(U(J,K+1,L)-2.*U(J,K,L)+U(J,K-1,L))
     2     +C3(J,K,L)*(U(J,K,L+1)-2.*U(J,K,L)+U(J,K,L-1))
     3+(U(J+1,K+1,L)-U(J-1,K+1,L)-U(J+1,K-1,L)+U(J-1,K-1,L))*C4(J,K,L)/2
     4+(U(J,K+1,L+1)-U(J,K-1,L+1)-U(J,K+1,L-1)+U(J,K-1,L-1))*C5(J,K,L)/2
     5+(U(J+1,K,L+1)-U(J-1,K,L+1)-U(J+1,K,L-1)+U(J-1,K,L-1))*C6(J,K,L)/2
     6+.5*C7(J,K,L)*(U(J+1,K,L)-U(J-1,K,L))+.5*C8(J,K,L)*(U(J,K+1,L)
     7-U(J,K-1,L))+.5*C9(J,K,L)*(U(J,K,L+1)-U(J,K,L-1))
      VL(J)=C1(J,K,L)*(V(J+1,K,L)-2.*V(J,K,L)+V(J-1,K,L))
     1     +C2(J,K,L)*(V(J,K+1,L)-2.*V(J,K,L)+V(J,K-1,L))
     2     +C3(J,K,L)*(V(J,K,L+1)-2.*V(J,K,L)+V(J,K,L-1))
     3+(V(J+1,K+1,L)-V(J-1,K+1,L)-V(J+1,K-1,L)+V(J-1,K-1,L))*C4(J,K,L)/2
     4+(V(J,K+1,L+1)-V(J,K-1,L+1)-V(J,K+1,L-1)+V(J,K-1,L-1))*C5(J,K,L)/2
     5+(V(J+1,K,L+1)-V(J-1,K,L+1)-V(J+1,K,L-1)+V(J-1,K,L-1))*C6(J,K,L)/2
     6+.5*C7(J,K,L)*(V(J+1,K,L)-V(J-1,K,L))+.5*C8(J,K,L)*(V(J,K+1,L)
     7-V(J,K-1,L))+.5*C9(J,K,L)*(V(J,K,L+1)-V(J,K,L-1))
      WL(J)=C1(J,K,L)*(W(J+1,K,L)-2.*W(J,K,L)+W(J-1,K,L))
     1     +C2(J,K,L)*(W(J,K+1,L)-2.*W(J,K,L)+W(J,K-1,L))
     2     +C3(J,K,L)*(W(J,K,L+1)-2.*W(J,K,L)+W(J,K,L-1))
     3+(W(J+1,K+1,L)-W(J-1,K+1,L)-W(J+1,K-1,L)+W(J-1,K-1,L))*C4(J,K,L)/2
     4+(W(J,K+1,L+1)-W(J,K-1,L+1)-W(J,K+1,L-1)+W(J,K-1,L-1))*C5(J,K,L)/2
     5+(W(J+1,K,L+1)-W(J-1,K,L+1)-W(J+1,K,L-1)+W(J-1,K,L-1))*C6(J,K,L)/2
     6+.5*C7(J,K,L)*(W(J+1,K,L)-W(J-1,K,L))+.5*C8(J,K,L)*(W(J,K+1,L)
     7-W(J,K-1,L))+.5*C9(J,K,L)*(W(J,K,L+1)-W(J,K,L-1))
   44 CONTINUE
      DO 45 J = 2,JM1
      PJ=(P(J+1,K,L)-P(J-1,K,L))*.5
      PK=(P(J,K+1,L)-P(J,K-1,L))*.5
      PL=(P(J,K,L+1)-P(J,K,L-1))*.5
      UINC=-UN(J)*.5-XX(J,K,L)*PJ-YX(J,K,L)*PK-ZX(J,K,L)*PL+REI*UL(J)
      VINC=-VN(J)*.5-XY(J,K,L)*PJ-YY(J,K,L)*PK-ZY(J,K,L)*PL+REI*VL(J)
      WINC=-WN(J)*.5-XZ(J,K,L)*PJ-YZ(J,K,L)*PK-ZZ(J,K,L)*PL+REI*WL(J)
      D(J,K,L)=U(J,K,L)+DT*UINC
      Q(J,K,L)=V(J,K,L)+DT*VINC
      S(J,K,L)=W(J,K,L)+DT*WINC
   45 CONTINUE
   40 CONTINUE
      DO 46 L = 2,LM1
      DO 46 K = 2,KM1
      DO 46 J = 2,JM1
      U(J,K,L)=D(J,K,L)
      V(J,K,L)=Q(J,K,L)
      W(J,K,L)=S(J,K,L)
   46 CONTINUE
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR OUTPUT
C**********************************************************************
      SUBROUTINE OUT
      PARAMETER (JD=22,KD=16,LD=11)
      COMMON/FLOW1/U(JD,KD,LD),V(JD,KD,LD),W(JD,KD,LD),P(JD,KD,LD)
      COMMON/FV1/JM,KM,LM,JM1,KM1,LM1,IL,NL,IW,NC
      DIMENSION LCN(10)
        LCN(1) = 2
        LCN(2) = 4
        LCN(3) = 6
        LCN(4) = 8
        LCN(5) = 10
      DO 91 LL = 1,5
        L = LCN(LL)
        WRITE(6,*) 'OUTPUT--- U,V,W,P --- L=',L
        DO 90 J = 2,JM
          WRITE(*,*) 'J=',J
          WRITE(*,600)(U(J,K,L),K=1,11  )
          WRITE(*,600)(V(J,K,L),K=1,11  )
          WRITE(*,600)(W(J,K,L),K=1,11  )
          WRITE(*,600)(P(J,K,L),K=1,11  )
   90   CONTINUE
   91 CONTINUE
  600 FORMAT(1H ,11F8.4)
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR THE SOLUTION OF POISSON EQ. FOR PRESSURE
C**********************************************************************
      SUBROUTINE PRS
      PARAMETER (JD=22,KD=16,LD=11)
      COMMON/MET1/XX(JD,KD,LD),XY(JD,KD,LD),XZ(JD,KD,LD)
     1           ,YX(JD,KD,LD),YY(JD,KD,LD),YZ(JD,KD,LD)
     2           ,ZX(JD,KD,LD),ZY(JD,KD,LD),ZZ(JD,KD,LD)
      COMMON/MET2/C1(JD,KD,LD),C2(JD,KD,LD),C3(JD,KD,LD)
     1           ,C4(JD,KD,LD),C5(JD,KD,LD),C6(JD,KD,LD)
     2           ,C7(JD,KD,LD),C8(JD,KD,LD),C9(JD,KD,LD)
      COMMON/FLOW1/U(JD,KD,LD),V(JD,KD,LD),W(JD,KD,LD),P(JD,KD,LD)
      COMMON/FLOW2/D(JD,KD,LD),Q(JD,KD,LD),S(JD,KD,LD)
      COMMON/FV1/JM,KM,LM,JM1,KM1,LM1,IL,NL,IW,NC
      COMMON/FV2/DT,DTI,REI,EPS,CONST
      COMMON/IBND/LFL
         DO 37 L = 2,LM1
         DO 37 K = 2,KM1
         DO 37 J = 2,JM1
           S(J,K,L) = .5/(C1(J,K,L)+C2(J,K,L)+C3(J,K,L))
   37    CONTINUE
C
C***  ITERATION
      DO 30 NN = 1,IL
C
C***  BOUNDARY CONDITION
C***  Y-Z
        DO 31 L = 1,LM
        DO 31 K = 1,KM
          P(1,K,L)  = P(JM-1,K,L)
          P(JM,K,L) = P(2,K,L)
   31   CONTINUE
C***  Y-Z
        DO 32 L = 1,LM
        DO 32 J = 1,JM
          P(J,1,L)  = P(J,2,L)
   32   CONTINUE
        DO 39 L = 1,LM
        DO 39 J = 1,JM
          P(J,KM,L)  = P(J,KM1,L)
   39 CONTINUE
C***  X-Y
      IF(LFL.EQ.2) THEN
        DO 53 K = 1,KM
        DO 53 J = 1,JM
          P(J,K,1)  = P(J,K,LM1)
          P(J,K,LM) = P(J,K,2)
   53   CONTINUE
      ELSE
        DO 33 K = 1,KM
        DO 33 J = 1,JM
          P(J,K,1)  = P(J,K,2)
          P(J,K,LM) = P(J,K,LM1)
   33   CONTINUE
      END IF
C
C***  SOR METHOD
      GOSA=0.
      DO 36 L = 2,LM1
      DO 36 K = 2,KM1
      DO 36 J = 2,JM1
      PA=((C1(J,K,L)+.5*C7(J,K,L))*P(J+1,K,L)+(C1(J,K,L)-.5*C7(J,K,L))
     1  *P(J-1,K,L)+(C2(J,K,L)+.5*C8(J,K,L))*P(J,K+1,L)+(C2(J,K,L)
     2  -.5*C8(J,K,L))*P(J,K-1,L)+(C3(J,K,L)+.5*C9(J,K,L))*P(J,K,L+1)
     3  +(C3(J,K,L)-.5*C9(J,K,L))*P(J,K,L-1)+C4(J,K,L)*.5
     4  *(P(J+1,K+1,L)-P(J-1,K+1,L)-P(J+1,K-1,L)+P(J-1,K-1,L))+.5
     5*(P(J,K+1,L+1)-P(J,K-1,L+1)-P(J,K+1,L-1)+P(J,K-1,L-1))*C5(J,K,L)
     6+(P(J+1,K,L+1)-P(J-1,K,L+1)-P(J+1,K,L-1)+P(J-1,K,L-1))*C6(J,K,L)
     7*.5 -Q(J,K,L))*S(J,K,L)-P(J,K,L)
      P(J,K,L)=P(J,K,L)+PA*CONST
      GOSA=GOSA+PA*PA
   36 CONTINUE
      IF(GOSA.LT.EPS) GO TO 34
   30 CONTINUE
   34  IF(NC-NC/2*2.EQ.0) WRITE(*,*) NC,NN,GOSA
      RETURN
      END
C**********************************************************************
C     SUBROUTINE FOR CALCULATING RHS OF POISSON EQ.
C**********************************************************************
      SUBROUTINE RHS
      PARAMETER (JD=22,KD=16,LD=11)
      COMMON/MET1/XX(JD,KD,LD),XY(JD,KD,LD),XZ(JD,KD,LD)
     1           ,YX(JD,KD,LD),YY(JD,KD,LD),YZ(JD,KD,LD)
     2           ,ZX(JD,KD,LD),ZY(JD,KD,LD),ZZ(JD,KD,LD)
      COMMON/FLOW1/U(JD,KD,LD),V(JD,KD,LD),W(JD,KD,LD),P(JD,KD,LD)
      COMMON/FLOW2/D(JD,KD,LD),Q(JD,KD,LD),S(JD,KD,LD)
      COMMON/FV1/JM,KM,LM,JM1,KM1,LM1,IL,NL,IW,NC
      COMMON/FV2/DT,DTI,REI,EPS,CONST
        DO 20 L = 2,LM1
        DO 20 K = 2,KM1
        DO 20 J = 2,JM1
          UJ = (U(J+1,K,L)-U(J-1,K,L))*.5
          UK = (U(J,K+1,L)-U(J,K-1,L))*.5
          UL = (U(J,K,L+1)-U(J,K,L-1))*.5
          VJ = (V(J+1,K,L)-V(J-1,K,L))*.5
          VK = (V(J,K+1,L)-V(J,K-1,L))*.5
          VL = (V(J,K,L+1)-V(J,K,L-1))*.5
          WJ = (W(J+1,K,L)-W(J-1,K,L))*.5
          WK = (W(J,K+1,L)-W(J,K-1,L))*.5
          WL = (W(J,K,L+1)-W(J,K,L-1))*.5
          UX = XX(J,K,L)*UJ+YX(J,K,L)*UK+ZX(J,K,L)*UL
          UY = XY(J,K,L)*UJ+YY(J,K,L)*UK+ZY(J,K,L)*UL
          UZ = XZ(J,K,L)*UJ+YZ(J,K,L)*UK+ZZ(J,K,L)*UL
          VX = XX(J,K,L)*VJ+YX(J,K,L)*VK+ZX(J,K,L)*VL
          VY = XY(J,K,L)*VJ+YY(J,K,L)*VK+ZY(J,K,L)*VL
          VZ = XZ(J,K,L)*VJ+YZ(J,K,L)*VK+ZZ(J,K,L)*VL
          WX = XX(J,K,L)*WJ+YX(J,K,L)*WK+ZX(J,K,L)*WL
          WY = XY(J,K,L)*WJ+YY(J,K,L)*WK+ZY(J,K,L)*WL
          WZ = XZ(J,K,L)*WJ+YZ(J,K,L)*WK+ZZ(J,K,L)*WL
C         D(J,K,L) = UX+VY+WZ
          Q(J,K,L) = -UX*UX-VY*VY-WZ*WZ-2.*(VX*UY+WY*VZ+UZ*WX)
     1               +(UX+VY+WZ)*DTI
   20   CONTINUE
      RETURN
      END
C**********************************************************************
C     SUBROTINE FOR SAVING DATA
C**********************************************************************
      SUBROUTINE SDATA
      PARAMETER (JD=22,KD=16,LD=11)
      COMMON/XYZ/X(JD,KD,LD),Y(JD,KD,LD),Z(JD,KD,LD)
      COMMON/FLOW1/U(JD,KD,LD),V(JD,KD,LD),W(JD,KD,LD),P(JD,KD,LD)
      COMMON/FV1/JM,KM,LM,JM1,KM1,LM1,IL,NL,IW,NC
        WRITE(3,*) JM,KM,LM
        DO 10 L = 1,LM
        DO 10 K = 1,KM
        DO 10 J = 1,JM
          WRITE(3,*) X(J,K,L),Y(J,K,L),Z(J,K,L)
     1              ,U(J,K,L),V(J,K,L),W(J,K,L),P(J,K,L)
   10   CONTINUE
        RETURN
      END
C**********************************************************************
C     SUBPOUTINE FOR OUTPUT                                           *
C**********************************************************************
      SUBROUTINE OUT2(L)
      PARAMETER ( JD=22,KD=16,LD=11)
      COMMON/XYZ/X(JD,KD,LD),Y(JD,KD,LD),Z(JD,KD,LD)
      COMMON/FLOW1/U(JD,KD,LD),V(JD,KD,LD),W(JD,KD,LD),P(JD,KD,LD)
      COMMON/FV1/JM,KM,LM,JM1,KM1,LM1,IL,NL,IW,NC
      DIMENSION PP(JD,KD)
      CHARACTER*2 ZA(39,31)
      DIMENSION INDEX(39,31)
C
      AREA(P1,Q1,P2,Q2,P3,Q3)=ABS((P2-P1)*(Q3-Q1)-(P3-P1)*(Q2-Q1))
C
      DO 20 J=1,JM
      PP(J,1)=0.
      DO 25 K=2,KM
      PP(J,K)=PP(J,K-1)+.5*(U(J,K,L)+U(J,K-1,L))*(Y(J,K,L)-Y(J,K-1,L))
     1                 -.5*(V(J,K,L)+V(J,K-1,L))*(X(J,K,L)-X(J,K-1,L))
   25 CONTINUE
   20 CONTINUE
      PMIN=PP(1,1)
      DO 10 K=1,KM
      DO 10 J=1,JM
        IF(PP(J,K).LT.PMIN) PMIN=PP(J,K)
   10 CONTINUE
      PMAX=PP(1,1)
      DO 30 K=1,KM
      DO 30 J=1,JM
        IF(PP(J,K).GT.PMAX) PMAX=PP(J,K)
   30 CONTINUE
      WRITE(*,*) 'PMAX=',PMAX,'   PMIN=',PMIN
      DO 40 K2=1,25
        WRITE(*,*) K2
      DO 40 J2=1,39
        INDEX(J2,K2)=0
        XA=FLOAT(J2-20)*.15*0.25
        YA=FLOAT(K2-13)*.15*1.2*0.25
        DO 60 K=1,KM-1
        DO 60 J=1,JM-1
          X1=X(J,K,L)
          Y1=Y(J,K,L)
          X2=X(J+1,K,L)
          Y2=Y(J+1,K,L)
          X3=X(J+1,K+1,L)
          Y3=Y(J+1,K+1,L)
          X4=X(J,K+1,L)
          Y4=Y(J,K+1,L)
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
            PJK=F1*PP(J+1,K+1)+F2*PP(J,K+1)+F3*PP(J,K)+F4*PP(J+1,K)
c            IND=INT((PJK-PMIN)/(PMAX-PMIN)*80)
            IND=INT((PJK-PMIN)/(PMAX-PMIN)*40)
            INDEX(J2,K2)=(MOD(IND,9)+1)*11
            GO TO 40
          END IF
   60   CONTINUE
   40 CONTINUE
      DO 70 K=1,25
      DO 70 J=1,39
        IND=INDEX(J,K)
        IF(IND.EQ.0) THEN
          ZA(J,K)='//'
        ELSE
          IF(MOD(IND,2).EQ.0) ZA(J,K)='..'
          IF(MOD(IND,2).EQ.1) ZA(J,K)='**'
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
        WRITE(*,601) (ZA(J,K),J=1,39)
   75 CONTINUE
  601   FORMAT(1H ,39A2)
      WRITE(*,*) 'PMAX=',PMAX,'   PMIN=',PMIN
      RETURN
      END
