C**********************************************************************
C     2 dimensional heat equation -- ADI METHOD --                    *
C**********************************************************************
      PARAMETER(NX=51,NY=51)
      DIMENSION A(60),B(60),C(60),D(60),U(NX,NY),UU(NX,NY)
      WRITE(*,*)'Input number of mesh for X & Y direction(<51) (20,20)'
       READ(*,*) KX,KY
         MX = KX + 1
         MY = KY + 1
      WRITE(*,*) 'Input time incremect dT? (0.02) '
       READ(*,*) DT 
      WRITE(*,*) 'Input number of time step? (100)'
       READ(*,*) NLAST
         DX = 1./(MX-1)
         DY = 1./(MY-1)
         R1 = .5*DT/DX**2
         R2 = .5*DT/DY**2
C
C***  Initial Condition
      DO 10 K = 1,MY
        DO 10 J = 1,MX
         U(J,K) = 0.
        UU(J,K) = 0.
   10 CONTINUE
C
C***  Time Marching
      DO 20 N = 1,NLAST
C
C***  Boundary Condition
        DO 30 K = 1,MY
           U(1,K) = 0.5
          U(MX,K) = 0.
          UU(1,K) = 0.5
         UU(MX,K) = 0.
   30   CONTINUE
        DO 35 J = 1,MX
           U(J,1) = 1.
          U(J,MY) = 0.
          UU(J,1) = 1.
         UU(J,MY) = 0.
   35   CONTINUE
C
C***  ADI method (X-direction)
        DO 40 K = 2,MY-1
          DO 45 J = 2,MX-1
             A(J) = -R1
             B(J) = 2.*R1+1.
             C(J) = -R1
             D(J) = U(J,K)+R2*(U(J,K+1)-2.*U(J,K)+U(J,K-1))
   45     CONTINUE
             D(2) = D(2)+R1*U(1,K)
          D(MX-1) = D(MX-1)+R1*U(MX,K)
C
          CALL THOMAS(2,MX-1,A,B,C,D)
C
         DO 50 J = 2,MX-1
         UU(J,K) = D(J)
   50    CONTINUE
   40 CONTINUE
C
C*** ADI method (Y-direction)
        DO 60 J = 2,MX-1
          DO 65 K = 2,MY-1
             A(K) = -R2
             B(K) = 2.*R2+1.
             C(K) = -R2
             D(K) = UU(J,K)+R1*(UU(J+1,K)-2.*UU(J,K)+UU(J-1,K))
   65     CONTINUE
             D(2) = D(2)+R2*UU(J,1)
          D(MY-1) = D(MY-1)+R2*UU(J,MY)
C
          CALL THOMAS(2 ,MY-1,A,B,C,D)
C
          DO 70 K = 2,MY-1
           U(J,K) = D(K)
   70     CONTINUE
   60 CONTINUE
C
      IF(MOD(N,50).EQ.0) CALL OUT(U,NX,NY,MX,MY)
C
   20 CONTINUE
C
      STOP
      END
C
C**********************************************************************
C     SUBROUTINE FOR SOLVING TRI-DIAGONAL MATRIX                      *
C**********************************************************************
      SUBROUTINE THOMAS(IL,IU,A,B,C,D)
      DIMENSION A(1),B(1),C(1),D(1)
C
      IP=IL+1
      DO 10 I = IP,IU
            R = C(I)/B(I-1)
         B(I) = B(I)-R*A(I-1)
         D(I) = D(I)-R*D(I-1)
   10 CONTINUE
C
      D(IU) = D(IU)/B(IU)
      DO 20 I = IP,IU
            J = IU-I+IL
         D(J) = (D(J)-A(J)*D(J+1))/B(J)
   20 CONTINUE
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