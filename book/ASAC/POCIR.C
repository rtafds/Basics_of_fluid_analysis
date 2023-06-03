#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/***********************************************************************
 *     UNSTEADY FLOW AROUND CIRCULAR CYLINDER
 *           PSI-OMEGA METHOD
 ***********************************************************************/
#define	MX	41
#define	MY	41
#define MX1     MX+1
#define MY1     MY+1

main()
{
	int	i, ii, j, k, kk, n, nmax, nx, ny;
	float	bb, const1, dt, dx, dx2, dxi, dy, dy2, dyi, eps, err, err1,
		fct, fff, omg[MX1][MY1], pai, psi[MX1][MY1], re, rei, rhs, tmp[MX1][MY1];
	FILE	*fp8;
L_99:
	printf( "INPUT NUMBER OF MESH NX  & NY (41,41) \n" );
	scanf( "%d,%d", &nx, &ny );
	printf( "INPUT REYNOLDS NUMBER RE (80) \n" );
	scanf( "%f", &re );
	printf( "INPUT TIME AND SPACE(RADIAL) INCREMENT DT&DY (.01,.1) \n" );
	scanf( "%f,%f", &dt, &dy );
	printf( "INPUT MAX. NUMBERS OF ITERATION FOR MAIN LOOP AND POISSON EQ.NMAX & KK\(400,40) n" );
	scanf( "%d,%d", &nmax, &kk );
	printf( "INPUT ACCELARATION PARAMETER CONST1 & MAX. ERROR EPS(1.0,.002) \n" );
	scanf( "%f,%f", &const1, &eps );

	pai = atan( 1. )*4.;
	dx = pai/(float)( nx - 1 );
	dxi = 1./dx;
	dyi = 1./dy;
	rei = 1./re;
	dx2 = dxi*dxi;
	dy2 = dyi*dyi;
	fct = 1./(2.*dx2 + 2.*dy2);

	/****  INITIAL CONDITION FOR PSI AND OMEGA */
	for( j = 1; j <= ny; j++ ){
		for( i = 1; i <= nx; i++ ){
			psi[i][j] = exp( (j - 1)*dy )*sin( dx*(i - 1) );
			omg[i][j] = 0.0;
		}
	}

	/****  MAIN LOOP
	 * */
	for( n = 1; n <= nmax; n++ ){
		fff = (n - 1)/30.;
		if( fff >= 1 )
			fff = 1.;

		/****  BOUNDARY CONDITION (STEP1)
		 ****  ON THE CYLINDER */
		for( i = 1; i <= nx; i++ ){
			omg[i][1] = -2.*psi[i][2]*dyi*dyi*fff;
			psi[i][1] = 0.;
		}
		/****  ON THE FAR BOUNDARY */
		for( i = 1; i <= nx; i++ ){
			psi[i][ny] = exp( (ny - 1)*dy )*sin( dx*(i - 1) );
			omg[i][ny] = 0.;
		}
		for( j = 1; j <= ny; j++ ){
			psi[1][j] = 0.;
			omg[1][j] = 0.;
			psi[MX][j] = 0.;
			omg[MX][j] = 0.;
		}
		/****  ALONG THE SYMMETRY LINE
		 *
		 ****  SOLVE POISSON EQUATION FOR PSI (STEP2) */
		fct = 1./(2.*dx2 + 2.*dy2);
		for( k = 1; k <= kk; k++ ){
			err = 0.;
			for( j = 2; j <= (ny - 1L); j++ ){
				for( i = 2; i <= (nx - 1); i++ ){
					rhs = ((psi[i+1][j] + psi[i-1][j])*dx2 + 
					 (psi[i][j+1] + psi[i][j-1])*dy2 + omg[i][j]*
					 exp( 2.*(j - 1)*dy ))*fct;
					err = err + (rhs - psi[i][j])*(rhs - psi[i][j]);
					psi[i][j] = psi[i][j]*(1. - const1) + rhs*
					 const1;
				}
			}
			if( err < 0.00001 )
				goto L_65;
		}
L_65:
		printf( "ITERATION NO. =%d    ERROR(L2) =%f \n", k, err );

		/****  CALCULATE NEW OMEGA (STEP3) */
		for( j = 2; j <= (ny - 1); j++ ){
			for( i = 2; i <= (nx - 1L); i++ ){

				tmp[i][j] = omg[i][j];

				rhs = ((omg[i+1][j] - 2.*omg[i][j] + omg[i-1][j])*
				 dx2 + (omg[i][j+1] - 2.*omg[i][j] + omg[i][j-1])*
				 dy2)*rei + ((psi[i+1][j] - psi[i-1][j])*
				 (omg[i][j+1] - omg[i][j-1]) - (psi[i][j+1] - 
				 psi[i][j-1])*(omg[i+1][j] - omg[i-1][j]))*
				 dxi*dyi/4.;
				omg[i][j] = omg[i][j] + dt*rhs*exp( -2.*(j - 1)*dy );
			}
		}

		err1 = 0.;
		for( j = 2; j <= (ny - 1); j++ ){
			for( i = 2; i <= (nx - 1); i++ ){
				bb = fabs( omg[i][j] - tmp[i][j] );
				if( bb >= err1 )
					err1 = bb;
			}
		}

		printf( "%d  *** ERROR(OMG)=%f   ***\n", n, err1 );
		if( n > 10 && err1 <= eps )
			goto L_90;

		}
	/****  END OF MAIN LOOP
	 * */
	printf( "NOT CONVERGE!  DO YOU WANT CONTINUE? (YES=1)\n" );
	scanf( "%d", &ii );
	if( ii == 1 )
		goto L_99;
L_90:
	out2( psi, &nx, &ny, &dy );
	fp8 = fopen( "file.8", "w" );
	for( j = 1; j <= ny; j++ ){
		for( i = 1; i <= nx; i++ ){
			fprintf( fp8, "%f %f \n", psi[i][j], omg[i][j] );
		}
	}
	fclose( fp8 );
	exit(0);
}

out2(a, nx, ny, dy)
float a[MX1][MY1];
int *nx, *ny;
float *dy;
{
	int	i, ii, ind, index[15][39], j, jj;
	float	aa, amax, amin, dx, pai, rr, rt, tet;

	pai = 4.*atan( 1. );
	dx = pai/(float)( *nx - 1 );

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

	for( j = 1; j <= 15; j++ ){
		for( i = 1; i <= 39; i++ ){
			ind = 0;
			if( i != 25 )
				rt = (float)( j - 1 )/fabs( (float)( i - 25 ) );
			tet = pai/2.;
			if( i <= 24 )
				tet = pai - atan( rt );
			if( i >= 26 )
				tet = atan( rt );
			rr = sqrt( (float)( (i - 25)*(i - 25) + (j - 1)*(j - 1) ) )/3.5;
			if( rr != 0. )
				jj = log( rr )/ *dy + 1;
			ii = tet/dx + 1.5;
			if( (ii >= 1 && ii <= *nx) && (jj >= 1 && jj <= *ny) ){
				aa = a[ii][jj]*100./amax;
				ind = aa + 2;
				if( aa < 0. )
					ind = 8;
				}
			index[i][j] = (ind%10)*11;
		}
	}
	for( j = 15; j >= 1; j-- ){
		printf( " " );
		for( i = 39; i >= 1; i-- ){
			printf( "%2d", index[i][j] );
		}
		printf( "\n" );
	}
	for( j = 2; j <= 15; j++ ){
		printf( " " );
		for( i = 39; i >= 1; i-- ){
			printf( "%2d", index[i][j] );
		}
		printf( "\n" );
	}
	return;
}
