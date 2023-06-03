C**********************************************************************
C     2 dimensional heat equation  -- Euler explicit method           *
C**********************************************************************
      PARAMETER(NX=51,NY=51)
      DIMENSION U(NX,NY),UU(NX,NY)
      WRITE(*,*)'Input number of mesh for X & Y direction(<50) (20,20)'
       READ(*,*) KX,KY
        MX = KX+1
        MY = KY+1
      WRITE(*,*) 'Input time increment dT  (.0005)'
       READ(*,*) DT
      WRITE(*,*) 'Input number of time step  (400)'
       READ(*,*) NLAST
        DX = 1./(MX-1)
        DY = 1./(MY-1)
        R1 = DT/DX**2
        R2 = DT/DY**2
C***  Initial Condition
        DO 10 K = 1,MY
          DO 10 J = 1,MX
           U(J,K) = 0.
          UU(J,K) = 0.
   10   CONTINUE
C
C***  Time Marching
      DO 20 N = 1,NLAST
C
C***  Boundary Condition
        DO 30 K = 1,MY
         U(1,K) = 0.5
         U(MX,K)= 0.
   30   CONTINUE
        DO 35 J = 1,MX
         U(J,1) = 1.
         U(J,MY)= 0.
   35   CONTINUE
C
C***  Euler explicit Method
        DO 40 K = 2,MY-1
          DO 40 J = 2,MX-1
           UU(J,K)= U(J,K)+R1*(U(J+1,K)-2.*U(J,K)+U(J-1,K))
     1              +R2*(U(J,K+1)-2.*U(J,K)+U(J,K-1))
   40   CONTINUE
        DO 45 K = 2,MY-1
          DO 45 J = 2,MX-1
           U(J,K) = UU(J,K)
   45   CONTINUE
C   
      IF(MOD(N,50).EQ.0) CALL OUT(U,NX,NY,MX,MY)
C
   20 CONTINUE
C
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