C
C     TRI-DIAGONAL MATRIX SOLVER
C
      SUBROUTINE THOMAS(IL,IU,A,B,C,D)
      DIMENSION A(1),B(1),C(1),D(1)
C
         IP   = IL+1
      DO 10 I = IP,IU
         R    = C(I)/B(I-1)
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
