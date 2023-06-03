C
C     BOUNDARY VALUE PROBLEM FOR ORDINALY DIFFERENCIAL EQUATIONS
C
      DIMENSION A(100),B(100),C(100),D(100)
C      
      WRITE(*,*) 'Input the number of mesh. (20)'
      READ(*,*) N
            H = 1./REAL(N)
C
      DO 10 I = 2,N
         A(I) = 1.0
         B(I) = H**2-2.0
         C(I) = 1.0
         D(I) =-REAL(I-1)*H**3
   10 CONTINUE
C
         CALL THOMAS(2,N,A,B,C,D)
C
      WRITE(*,600)
  600 FORMAT(1H ,7X,'X',7X,'APPROX.',5X,'EXACT',5X,'ERROR')
C
      DO 20 I = 2,N
            X = H*REAL(I-1)
            U = SIN(X)/SIN(1.0)-X
          ERR = (D(I)-U)/U*100.0
        WRITE(*,601) X,D(I),U,ERR
   20 CONTINUE
C
  601 FORMAT(1H ,4F11.6)   
C
      STOP
      END
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