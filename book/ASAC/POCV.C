#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/***********************************************************************
 *     UNSTEADY FLOW IN CUBIC CAVITY
 *           PSI-OMEGA METHOD
 ***********************************************************************/
#define	MX	51
#define	MY	51
#define MX1     MX+1
#define MY1     MY+1

main()
{
	int	i, isave, j, k, n, nmax, nn;
	float	aa, bb, const1, const2, dt, eps, err2, h, hh, hi, omg[MX1][MY1],
		psi[MX1][MY1], re, rhs, tmp[MX1][MY1];
	FILE	*fp17;

L_99:
	printf( "INPUT NUMBER OF MESH NN (<52) (20) \n" );
	scanf( "%d", &nn );
	if( nn <= 1 || nn >= 52 )
		goto L_99;
	printf( "INPUT REYNOLDS NUMBER RE & DELTA T (40,0.01) \n" );
	scanf( "%f,%f", &re, &dt );
	printf( "INPUT MAXIMUM NUMBER OF ITERATION NMAX (100) \n" );
	scanf( "%d", &nmax );
	printf( "INPUT ACCELARATION PARAMETER CONST1 & CONST2 (1,1) \n" );
	scanf( "%f,%f", &const1, &const2 );
	printf( "INPUT MAXIMUM ERROR EPS (0.00001) \n" );
	scanf( "%f", &eps );
	h = 1./(float)( nn - 1 );
	hi = 1./h;
	hh = h*h;

	/****  INITIAL CONDITION FOR PSI AND OMEGA */
	for( j = 1; j <= nn; j++ ){
		for( i = 1; i <= nn; i++ ){
			psi[i][j] = 0.0;
			omg[i][j] = 0.0;
		}
	}

	/****  MAIN LOOP
	 * */
	for( n = 1; n <= nmax; n++ ){

		/****  BOUNDARY CONDITION (STEP1)
		 ****  LEFT AND RIGHT */
		for( j = 1; j <= nn; j++ ){
			omg[1][j] = -2.*psi[2][j]*hi*hi;
			omg[nn][j] = -2.*psi[nn-1][j]*hi*hi;
		}
		/****  BOTTOM AND TOP */
		for( i = 1; i <= nn; i++ ){
			omg[i][1] = -2.*psi[i][2]*hi*hi;
			omg[i][nn] = -2.*(psi[i][nn-1] + h)*hi*hi;
		}

		/****  CALCULATE NEW PSI (STEP2) */
		for( k = 1; k <= 100; k++ ){
			for( j = 2; j <= (nn - 1); j++ ){
				for( i = 2; i <= (nn - 1); i++ ){

					tmp[i][j] = psi[i][j];

					rhs = (psi[i+1][j] + psi[i-1][j] + psi[i][j+1] +
					 psi[i][j-1])/4. + omg[i][j]*h*h/4. - psi[i][j];
					psi[i][j] = psi[i][j] + const2*rhs;
				}
			}

			err2 = 0.;
			for( j = 2; j <= (nn - 1); j++ ){
				for( i = 2; i <= (nn - 1); i++ ){
			/*  aa = max( 1e-8, fabs( tmp[i][j] ) ); */
					bb = fabs( psi[i][j] - tmp[i][j] );
					if( bb >= err2 )
						err2 = bb;
				}
			}
			if( err2 < eps )
				goto L_75;
		}
L_75:
		;
		fprintf( stdout, "%d %d %f \n", n, k, err2 );

		/****  CALCULATE NEW OMEGA (STEP3) */
		for( j = 2; j <= (nn - 1); j++ ){
			for( i = 2; i <= (nn - 1); i++ ){

				tmp[i][j] = omg[i][j];

				rhs = ((psi[i+1][j] - psi[i-1][j])*(omg[i][j+1] - 
				      omg[i][j-1]) - (psi[i][j+1] - psi[i][j-1])*
				      (omg[i+1][j] - omg[i-1][j]))/4. + (omg[i+1][j] + 
				      omg[i-1][j] + omg[i][j+1] + omg[i][j-1] - 
				      4.*omg[i][j])/re;
				omg[i][j] = omg[i][j] + dt*rhs/hh;
			}
		}

	}
	/****  END OF MAIN LOOP */
	out1( psi, &nn );
	printf( "SAVING IN FILE?  YES=1\n" );
	scanf( "%d", &isave );
	if( isave == 1 ){
		fp17 = fopen( "file.17", "w" );
		fprintf( fp17, "%d %d \n", nn, nn );
		for( j = 1; j <= nn; j++ ){
			for( i = 1; i <= nn; i++ ){
				fprintf( fp17, "%f ", psi[i][j] );
			}
		}
		for( j = 1; j <= nn; j++ ){
			for( i = 1; i <= nn; i++ ){
				fprintf( fp17, "%f ", omg[i][j] );
			}
		}
		fprintf( fp17, "\n" );
		fclose ( fp17 );
	}
	exit(0);
}

out1(a, nn)
float a[MX1][MY1];
int *nn;
{
	int	i, iw[81], j;
	float	amax, amin;

	amin = a[1][1];
	for( j = 1; j <= *nn; j++ ){
		for( i = 1; i <= *nn; i++ ){
			if( a[i][j] < amin )
				amin = a[i][j];
		}
	}
	for( j = 1; j <= *nn; j++ ){
		for( i = 1; i <= *nn; i++ ){
			a[i][j] = a[i][j] - amin;
		}
	}
	amax = a[1][1];
	for( j = 1; j <= *nn; j++ ){
		for( i = 1; i <= *nn; i++ ){
			if( a[i][j] > amax )
				amax = a[i][j];
		}
	}
	for( j = *nn; j >= 1; j-- ){
		for( i = 1; i <= *nn; i++ ){
			iw[i] = (long)( a[i][j]*9.99999/amax );
		}
		printf( " " );
		for( i = 1; i <= *nn; i++ ){
			printf( "%1d", iw[i] );
		}
		printf( "\n" );
	}
	return;
}
