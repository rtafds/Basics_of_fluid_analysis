#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/***********************************************************************
 *     UNSTEADY FLOW IN CUBIC CAVITY
 *           PSI-OMEGA METHOD
 ***********************************************************************/
#define	MX	21
#define	MY	21
#define MX1     MX+1
#define MY1     MY+1

main()
{
	int	i, isave, j, k, n, nmax, nn, nx, ny, ja, jb;
	float	a1[MX1], a2[MX1], a3[MY1], a4[MY1], ar, b1[MX1], b2[MX1], b3[MY1],
		b4[MY1], bb, bx, c1[MX1], c2[MX1], c3[MY1], c4[MY1], const1, const2,
		dt, eps, ff, g1, h, hh, hi, omg[MX1][MY1], pr, psi[MX1][MY1], re,
		rhs, t[MX1][MY1], tv, tx, ty, uu, vv, x[MX1], x1, x2, x3, x4, y[MY1],
		y1, y2, y3, y4, ps0;
	FILE	*fp17;

L_99:
	printf( "INPUT NUMBER OF MESH NX,NY  (21,21)  \n" );
	scanf( "%d,%d", &nx, &ny );
	printf( "INPUT REYNOLDS NUMBER RE & DELTA T && AR? (40,.01,.2) \n" );
	scanf( "%f,%f,%f", &re, &dt, &ar );
	printf( "INPUT MAXIMUM NUMBER OF ITERATION NMAX (400) \n" );
	scanf( "%d", &nmax );
	printf( "INPUT ACCELARATION PARAMETER CONST1 & CONST2 (1,1) \n" );
	scanf( "%f,%f", &const1, &const2 );
	printf( "INPUT MAXIMUM ERROR EPS (0.00001) \n" );
	scanf( "%f", &eps );
	printf( "INPUT CONCENTRATION PARAMETER? BB? (0.98) \n" );
	scanf( "%f", &bb );
	printf( "INPUT JA,JB --ENTRANCE(JA-NX) EXIT(1-JB) (16,5) \n" );
	scanf( "%d,%d", &ja, &jb );
	pr = 0.71;
	h = 1./(float)( nx - 1 );
	hi = 1./h;
	hh = h*h;

	/****  GRID */
	ff = (exp( bb ) + 1.)/(exp( bb ) - 1.);
	for( i = 1; i <= nx; i++ ){
		bx = bb*(float)( i - 1 )/(float)( nx - 1 );
		x[i] = ff*(exp( bx ) - 1.)/(exp( bx ) + 1.);
	}
	for( j = 1; j <= ny; j++ ){
		bx = bb*(float)( j - 1 )/(float)( ny - 1 );
		y[j] = ff*(exp( bx ) - 1.)/(exp( bx ) + 1.);
	}
	h = x[2] - x[1];
	hi = 1./h;
	ps0 = y[ny]-y[ja];
	for( i = 2; i <= (nx - 1); i++ ){
		x1 = x[i+1] - x[i];
		x2 = x[i] - x[i-1];
		x3 = x[i+1] - x[i-1];
		x4 = x[i+1] - 2.*x[i] + x[i-1];
		a1[i] = -x1/(x2*x3);
		b1[i] = x4/(x1*x2);
		c1[i] = x2/(x1*x3);
		a2[i] = 2./(x2*x3);
		b2[i] = -2./(x1*x2);
		c2[i] = 2./(x1*x3);
	}
	for( j = 2; j <= (ny - 1); j++ ){
		y1 = y[j+1] - y[j];
		y2 = y[j] - y[j-1];
		y3 = y[j+1] - y[j-1];
		y4 = y[j+1] - 2.*y[j] + y[j-1];
		a3[j] = -y1/(y2*y3);
		b3[j] = y4/(y1*y2);
		c3[j] = y2/(y1*y3);
		a4[j] = 2./(y2*y3);
		b4[j] = -2./(y1*y2);
		c4[j] = 2./(y1*y3);
	}

	/****  CALCULATE METRICS
	 *
	 ****  INITIAL CONDITION FOR PSI AND OMEGA */
	for( j = 1; j <= ny; j++ ){
		for( i = 1; i <= nx; i++ ){
			psi[i][j] = 0.0;
			omg[i][j] = 0.0;
			t[i][j] = 0.0;
		}
	}
	for( i = 1; i <= nx; i++ ){
			psi[i][ny] = ps0;
		}
	for( j = jb; j <= ny; j++ ){
			psi[nx][j] = ps0;
		}

	/****  MAIN LOOP
	 * */
	for( n = 1; n <= nmax; n++ ){

		/****  BOUNDARY CONDITION (STEP1)
		 ****  LEFT AND RIGHT */
		for( j = 1; j <= ny; j++ ){
			omg[j][1] = (a2[1]*c1[1]/a1[1]-c2[1])*psi[2][j];
			omg[nx][j] = (c2[nx]*a1[nx]/c1[nx]-a2[nx])*psi[nx-1][j]
				    +(c2[nx]*b1[nx]/c1[nx]-b2[nx])*ps0;
			psi[nx][j] = ps0;
			t[1][j] = t[2][j];
			t[nx][j] = t[nx-1][j];
		}
		for( j = ja; j <= ny; j++ ){
			omg[1][j] = 0.;
			psi[1][j] = (y[j]-y[ja])/(y[ny]-y[ja])*ps0;
	    	}
		for( j = 1; j <= jb; j++ ){
		    	omg[nx][j] = omg[nx-1][j];
	    		psi[nx][j] = (y[j]-y[1])/(y[jb]-y[1])*ps0;
		}
		/****  BOTTOM AND TOP */
		for( i = 1; i <= nx; i++ ){
			omg[i][1] = (a4[1]*c3[1]/a3[1]-c4[1])*psi[i][2];
			omg[i][ny] = (c4[ny]*a3[ny]/c3[ny]-a4[ny])*psi[i][ny-1]
			            +(c4[ny]*b3[ny]/c3[ny]-b4[ny])*ps0;
			t[i][1] = 1.;
			t[i][ny] = 0.;
		}

		/****  CALCULATE NEW PSI (STEP2) */
		for( k = 1; k <= 100; k++ ){
			g1 = 0.;
			for( j = 2; j <= (ny - 1); j++ ){
				for( i = 2; i <= (nx - 1); i++ ){
					rhs = -(a2[i]*psi[i-1][j] + c2[i]*psi[i+1][j] +
					      a4[j]*psi[i][j-1] + c4[j]*psi[i][j+1] +
					      omg[i][j])/(b2[i] + b4[j]) - psi[i][j];
					g1 = g1 + rhs*rhs;
					psi[i][j] = psi[i][j] + const2*rhs;
				}
			}

			if( g1 < eps )
				goto L_75;
		}
L_75:
		;
		if( (n%20) == 0 ){
			printf( "%d %d %f \n", n, k, g1 );
		}

		/****  CALCULATE NEW OMEGA (STEP3) */
		for( j = 2; j <= (ny - 1); j++ ){
			for( i = 2; i <= (nx - 1); i++ ){

				rhs = -(a3[j]*psi[i][j-1] + b3[j]*psi[i][j] +
				      c3[j]*psi[i][j+1])*(a1[i]*omg[i-1][j] +
				      b1[i]*omg[i][j] + c1[i]*omg[i+1][j]) + (a3[j]*
				      omg[i][j-1] + b3[j]*omg[i][j] + c3[j]*omg[i][j+1])*
				      (a1[i]*psi[i-1][j] + b1[i]*psi[i][j] + c1[i]*
				      psi[i+1][j]) + (a2[i]*omg[i-1][j] + b2[i]*
				      omg[i][j] + c2[i]*omg[i+1][j] + a4[j]*omg[i][j-1] +
				      b4[j]*omg[i][j] + c4[j]*omg[i][j+1])/re -
				      ar*(a1[i]*t[i-1][j] + b1[i]*t[i][j] + c1[i]*
				      t[i+1][j]);
				omg[i][j] = omg[i][j] + dt*rhs;
			}
		}

		/****  CALCULATE NEW TEMPERATURE (STEP4) */
		for( j = 2; j <= (ny - 1); j++ ){
			for( i = 2; i <= (nx - 1); i++ ){

				uu = a3[j]*psi[i][j-1] + b3[j]*psi[i][j] +
				     c3[j]*psi[i][j+1];
				vv = -(a1[i]*psi[i-1][j] + b1[i]*psi[i][j] +
				     c1[i]*psi[i+1][j]);
				tx = a1[i]*t[i-1][j] + b1[i]*t[i][j] + c1[i]*
				     t[i+1][j];
				ty = a3[j]*t[i][j-1] + b3[j]*t[i][j] + c3[j]*
				     t[i][j+1];
				tv = (a2[i]*t[i-1][j] + b2[i]*t[i][j] + c2[i]*
				     t[i+1][j] + a4[j]*t[i][j-1] + b4[j]*t[i][j] +
				     c4[j]*t[i][j+1])/(re*pr);
				t[i][j] = t[i][j] + dt*(-uu*tx - vv*ty + tv);
			}
		}


	}
	/****  END OF MAIN LOOP */
	out1( t, &nx, &ny );
	printf( "\n" );
	out1( psi, &nx, &ny );
	printf( "SAVING IN FILE?  YES=1\n" );
	scanf( "%d", &isave );
	if( isave == 1 ){
		fp17 = fopen( "file.17", "w" );
		fprintf( fp17, "%d %d \n", nx, ny );
		for( j = 1; j <= ny; j++ ){
			for( i = 1; i <= nx; i++ ){
				fprintf( fp17, "%f", psi[i][j] );
			}
		}
		for( j = 1; j <= ny; j++ ){
			for( i = 1; i <= nx; i++ ){
				fprintf( fp17, "%f", omg[i][j] );
			}
		}
		fprintf( fp17, "\n" );
		fclose( fp17 );
	}
	exit(0);
}

out1(a, nx, ny)
float a[MX1][MY1];
int *nx, *ny;
{
	int	i, iw[81], j;
	float	amax, amin;

	amin = a[1][1];
	for( j = 1; j <= *ny; j++ ){
		for( i = 1; i <= *nx; i++ ){
			if( a[i][j] < amin )
				amin = a[i][j];
		}
	}
	for( j = 1; j <= *ny; j++ ){
		for( i = 1; i <= *nx; i++ ){
			a[i][j] = a[i][j] - amin;
		}
	}
	amax = a[1][1];
	for( j = 1; j <= *ny; j++ ){
		for( i = 1; i <= *nx; i++ ){
			if( a[i][j] > amax )
				amax = a[i][j];
		}
	}
	for( j = *ny; j >= 1; j-- ){
		for( i = 1; i <= *nx; i++ ){
			iw[i] = ( a[i][j]*9.99999/amax );
			iw[i] = iw[i]*11;
		}
		printf( " " );
		for( i = 1L; i <= *nx; i++ ){
			printf( "%2d", iw[i] );
		}
		printf( "\n" );
	}
	return;
}
