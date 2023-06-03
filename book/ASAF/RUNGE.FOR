C**** FOURTH ORDER RUNGE-KUTTA METHOD
        F(X,Y) = .5*(1.+X)*Y*Y
C
        WRITE(*,*) 'INPUT H & XLAST ? (0.01,1.0)'
         READ(*,*)  H,XLAST
         X = 0.
         Y = 1.
C
    1    YEXACT = 4./(4.-2.*X-X*X)
            ERR = (YEXACT-Y)/YEXACT*100.
            WRITE(*,600) X,Y,YEXACT,ERR
              DEL1 = H*F(X,Y)
              DEL2 = H*F(X+H/2.,Y+DEL1/2.)
              DEL3 = H*F(X+H/2.,Y+DEL2/2.)
              DEL4 = H*F(X+H,Y+DEL3)
           Y = Y+(DEL1+2.*DEL2+2.*DEL3+DEL4)/6.
           X = X+H
        IF(X.LE.XLAST) GO TO 1
  600   FORMAT(1H ,F7.2,3F14.7)
       STOP
      END
