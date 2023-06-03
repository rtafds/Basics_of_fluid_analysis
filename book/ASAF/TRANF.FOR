C**********************************************************************
C     TRANSFINITE INTERPOLATION (GRID INSIDE CIRCULAR CYLINDER)       *
C**********************************************************************
      PARAMETER(MX=11,MY=9)
      COMMON X(MX,MY),Y(MX,MY)
         PAI  =  4.*ATAN(1.)
         TET  =  PAI/18.
C
C***  GRID POINTS ON THE BOUNDARY
         DO 10  I = 1,MX
               TE = TET*(I+17)
           X(I,1) = COS(TE)
           Y(I,1) = SIN(TE)
               TF = TET*(11-I)
          X(I,MY) = COS(TF)
          Y(I,MY) = SIN(TF)
   10    CONTINUE
         DO 20  J = 1,MY
               TE = TET*(19-J)
           X(1,J) = COS(TE)
           Y(1,J) = SIN(TE)
               TF = TET*(J+27)
          X(MX,J) = COS(TF)
          Y(MX,J) = SIN(TF)
   20    CONTINUE
C
C***  TRANSFINITE INTERPOLATION
         DO 30  J = 2,MY-1
         DO 30  I = 2,MX-1
                A = FLOAT(I-MX)/FLOAT(1-MX)
                B = FLOAT(J-MY)/FLOAT(1-MY)
           X(I,J) = A*X(1,J)+(1-A)*X(MX,J)+B*X(I,1)+(1-B)*X(I,MY)
     1             -A*B*X(1,1)-A*(1-B)*X(1,MY)-(1-A)*B*X(MX,1)
     2             -(1-A)*(1-B)*X(MX,MY)
           Y(I,J) = A*Y(1,J)+(1-A)*Y(MX,J)+B*Y(I,1)+(1-B)*Y(I,MY)
     1             -A*B*Y(1,1)-A*(1-B)*Y(1,MY)-(1-A)*B*Y(MX,1)
     2             -(1-A)*(1-B)*Y(MX,MY)
   30    CONTINUE
C
         WRITE(*,*) 'Save Data? (Yes=1,No=0) '
         READ(*,*) ISAVE
          IF(ISAVE.EQ.1) THEN
            WRITE(8,*) MX,MY
              DO 35  J = 1,MY
              DO 35  I = 1,MX
                WRITE(8,*) X(I,J),Y(I,J)
   35         CONTINUE
          END IF
        MZ=1
        WRITE(8) MX,MY,MZ
        WRITE(8) (((X(I,J),I=1,MX),J=1,MY),K=1,MZ)
        WRITE(8) (((Y(I,J),I=1,MX),J=1,MY),K=1,MZ)
        WRITE(8) (((0.0,I=1,MX),J=1,MY),K=1,MZ)
C***  PROGRAM FOR PLOTTING
c      CALL PLOTS(0)
      WRITE(*,*) 'FCT,X11,Y11'
       READ(*,*) FCT,X11,Y11
c      CALL FACTOR(FCT)
c      CALL PLOT(0.,0.,-3)
      DO 40 J = 1,MY
      DO 40 I = 1,MX
           IC = 2
           IF(I.EQ.1) IC=3
           XA = X(I,J)+X11
           YA = Y(I,J)+Y11
c           CALL PLOT(XA,YA,IC)
   40 CONTINUE
      DO 50 I = 1,MX
      DO 50 J = 1,MY
           IC = 2
           IF(J.EQ.1) IC=3
           XA = X(I,J)+X11
           YA = Y(I,J)+Y11
c           CALL PLOT(XA,YA,IC)
   50 CONTINUE
c        CALL PLOT(0.,0.,-3)
c        CALL PLOT(0.,0.,999)
      STOP
      END
