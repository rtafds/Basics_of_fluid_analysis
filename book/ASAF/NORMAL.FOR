C***********************************************************************
C     SUBROUTINE FOR MAKING NEARLY ORTHOGONAL GRID                     *
C***********************************************************************
      SUBROUTINE NORMAL(X,Y,T,MX,MY,NX,NY)
      DIMENSION X(MX,MY),Y(MX,MY),T(MX,2)
       DO 10 J = 2,NY
         DO 20 I = 2,NX-1
              T1 =-1000.
              IF(ABS(Y(I+1,J)-Y(I-1,J)).GE.0.00001) 
     1        T1 =-(X(I+1,J)-X(I-1,J))/(Y(I+1,J)-Y(I-1,J))
           DO 30 KK = 1,NX-1
             II = I
             DO 31 K1 = 1,2
               T2 = 1000.
               IF(II.LE.0.OR.II.GE.MX) GO TO 31
               IF(ABS(X(II+1,J)-X(II,J)).GE.0.00001)
     2         T2 = (Y(II+1,J)-Y(II,J))/(X(II+1,J)-X(II,J))
               A1 = (T2*X(II,J)-T1*X(I,J-1)+Y(I,J-1)-Y(II,J))/(T2-T1)
               B1 = T2*(A1-X(II,J))+Y(II,J)
               QQ = SQRT((X(II+1,J)-X(II,J))**2+(Y(II+1,J)-Y(II,J))**2)
               Q1 = SQRT((A1-X(II,J))**2+(B1-Y(II,J))**2)
               Q2 = SQRT((X(II+1,J)-A1)**2+(Y(II+1,J)-B1)**2)
               IF(ABS(QQ-Q1-Q2).LE.0.0001) GO TO 40
               IF(K1.EQ.1) II = I-KK
               IF(K1.EQ.2) II = I+KK
   31        CONTINUE
   30      CONTINUE
   40      CONTINUE
           T(I,1) = .5*(X(I,J)+A1)
           T(I,2) = .5*(Y(I,J)+B1)
   20    CONTINUE
         DO 50 II = 2,MX-1
           X(II,J) = T(II,1)
           Y(II,J) = T(II,2)
   50    CONTINUE
   10  CONTINUE
      RETURN
      END
