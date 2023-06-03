C**********************************************************************
C     LAPLACE EQUATION                                                *
C**********************************************************************
      PARAMETER(NX=51,NY=51)
      DIMENSION U(NX,NY),UU(NX,NY)
C
C**** INPUT & CALCULATE PARAMETERS
      WRITE(*,*) 'Input number of mesh (KX<51,KY<51). (20,20)'
       READ(*,*) KX,KY
       MX = KX + 1
       MY = KY + 1
      WRITE(*,*) 'Input number of iteration for Laplace equation.(1000)'
       READ(*,*) NN
      WRITE(*,*) 'Input maximum error (EPSIRON). (0.00001)'
       READ(*,*) EPS
C
C**** INITIAL CONDITION
       DO 10 J = 1,MY
       DO 10 I = 1,MX
        U(I,J) = 0.
   10  CONTINUE
C
C**** BOUNDARY CONDITION
       DO 20 I = 1,MX
        U(I,1) = 1.0
        U(I,MX)= 0.0
   20  CONTINUE
       DO 30 J = 1,MY
        U(1,J) = 0.5
        U(MX,J)= 0.0
   30 CONTINUE
C
C**** MAIN LOOP
       DO 40 N = 1,NN
         DO 50 J = 2,MY-1
         DO 50 I = 2,MX-1
         UU(I,J) = U(I,J)
          U(I,J) = .25*(U(I+1,J)+U(I,J+1)+U(I-1,J)+U(I,J-1))
   50    CONTINUE
C
           ERR = 0.
         DO 60 J = 2,MY-1
         DO 60 I = 2,MX-1
              BB = ABS((UU(I,J)-U(I,J))/U(I,J))
              IF(BB.GT.ERR) ERR = BB
   60    CONTINUE
           IF(ERR.LT.EPS) GO TO 70
           IF(MOD(N,50).EQ.0) WRITE(*,*) 'N=',N,'ERROR=',ERR
   40 CONTINUE
C
C**** OUTPUT
   70    CALL OUT(U,NX,NY,MX,MY)
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