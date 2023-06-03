C**********************************************************************
C     SUBROUTINE FOR REARRANGING GRID                                 *
C**********************************************************************
      SUBROUTINE ARRANG(X,Y,T,S,Q,MX,MY,NX,NY,ITYP)
      DIMENSION X(MX,MY),Y(MX,MY),T(MY,2),S(MY),Q(MY)
        ID   =  NX
        JD   =  NY
      DO 20 I = 1,ID
        S(1)  = 0.0
        DO 21 J = 2,JD
         S(J) = S(J-1)+SQRT((X(I,J)-X(I,J-1))**2+(Y(I,J)-Y(I,J-1))**2)
   21   CONTINUE
        T(JD,1) = X(I,JD)
        T(JD,2) = Y(I,JD)
        T(1,1)  = X(I,1)
        T(1,2)  = Y(I,1)
C
        IF(ITYP.EQ.1) THEN
C
          FCTA  = 1.15
          AA    = S(JD)*(FCTA-1.)/(FCTA**(JD-1)-1.)
          Q(1)  = 0.
          DO 30 J = 2,JD
          Q(J)  = Q(J-1)+AA*FCTA**FLOAT(J-2)
   30     CONTINUE
C
        ELSE
C
          DO 35 J = 2,JD
          Q(J)  = S(JD)*FLOAT(J-1)/FLOAT(JD-1)
   35     CONTINUE
C
        END IF
C
          DO 22 JJ = 2,JD-1
                JP = JD
             DO 23 J = 2,JD
             IF(Q(JJ).GT.S(J)) GO TO 23
                JP = J
             GO TO 24
   23        CONTINUE
   24     RTO     = (Q(JJ)-S(JP-1))/(S(JP)-S(JP-1))
          T(JJ,1) = (1.-RTO)*X(I,JP-1)+RTO*X(I,JP)
          T(JJ,2) = (1.-RTO)*Y(I,JP-1)+RTO*Y(I,JP)
   22     CONTINUE
          DO 40 J = 1,JD
           X(I,J) = T(J,1)
           Y(I,J) = T(J,2)
   40     CONTINUE
   20 CONTINUE
      RETURN
      END
