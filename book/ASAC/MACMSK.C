#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define	ID	21
#define	JD	21
#define ID1     ID+1
#define JD1     JD+1

float	u[ID1][JD1], v[ID1][JD1], p[ID1][JD1], q[ID1][JD1], d[ID1][JD1],
s[ID1][JD1];
int	ifl[ID1][JD1], ibd[ID1][JD1], iwt[ID1][JD1];
float	ul[ID1], x[ID1], xa[ID1], xb[ID1], xc[ID1], xd[ID1], xg[ID1], y[JD1],
	ya[JD1], yb[JD1], yc[JD1], yd[JD1], yg[JD1];

main()
{
	int	i, i11, i19, i20, i21, i35, iqq, ir, itest, iwrite, j,
		j11, j19, j20, j21, j35, k21, kk,  km, l, lm;
	float	dt, g1, ok, r1, re, stl, td, u1, u2, unx, uny, uv, v1,
		v2, vnx, vny, vv;
	FILE	*fp35, *fp11, *fp10;

	printf( "DT,RE,LM,IR (0.005,40,400,0) \n" );
	scanf( "%f,%f,%d,%d", &dt, &re, &lm, &ir );
	i21 = ID;
	j21 = JD;
	i11 = (i21 + 1)/2;
	j11 = (j21 + 1)/2;
	td = 1./dt;
	r1 = 1./re;
	i20 = i21 - 1;
	i19 = i21 - 2;
	j20 = j21 - 1;
	j19 = j21 - 2;
	km = 10;
	itest = 0;

	if( itest == 0 ){

		for( i = 1; i <= i21; i++ ){
			x[i] = (float)( i - i11 )/(float)( i21 - i11 );
		}
		for( j = 1; j <= j21; j++ ){
			y[j] = (float)( j - j11 )/(float)( j21 - j11 );
		}

		for( j = 1; j <= j21; j++ ){
			for( i = 1; i <= i21; i++ ){
				stl = sqrt( x[i]*x[i] + y[j]*y[j] );
				ifl[i][j] = 1;
				if( stl < .5 )
					ifl[i][j] = 0;
				ibd[i][j] = 0;
				iwt[i][j] = 0;
			}
		}

	}
	else{

		fp35 = fopen( "file.35", "r" );
		fscanf( fp35, "%i %i ", &i35, &j35 );
		fscanf( fp35, "" );
		for( j = 1; j <= j35; j++ ){
			for( i = 1; i <= i35; i++ ){
				fscanf( fp35, "%+ %i ", &ifl[i][j] );
			}
		}
		fscanf( fp35, "" );
		for( i = 1; i <= i35; i++ ){
			fscanf( fp35, "%+ %i ", &x[i] );
		}
		fscanf( fp35, "" );
		for( j = 1; j <= j35; j++ ){
			fscanf( fp35, "%+ %i ", &y[i] );
		}
		for( j = 1; j <= j21; j++ ){
			for( i = 1; i <= i21; i++ ){
				ibd[i][j] = 0;
				iwt[i][j] = 0;
			}
		}
		fclose( fp35 );
	}

	for( i = 2; i <= i20; i++ ){
		xg[i] = 2./(x[i+1] - x[i-1]);
		xa[i] = xg[i]*xg[i];
		xb[i] = (x[i+1] - 2.*x[i] + x[i-1])*xg[i]*xg[i]*xg[i] ;
		xc[i] = xa[i] - xb[i]*.5;
		xd[i] = xa[i] + xb[i]*.5;
	}
	for( j = 2; j <= j20; j++ ){
		yg[j] = 2./(y[j+1] - y[j-1]);
		ya[j] = yg[j]*yg[j];
		yb[j] = (y[j+1] - 2.*y[j] + y[j-1])*yg[j]*yg[j]*yg[j] ;
		yc[j] = ya[j] - yb[j]*.5;
		yd[j] = ya[j] + yb[j]*.5;
	}

	for( j = 1; j <= j21; j++ ){
		for( i = 1; i <= (i21 - 1); i++ ){
			if( ifl[i+1][j] - ifl[i][j] == -1 )
				ibd[i][j] = 1;
			}
		for( i = i21; i >= 2; i-- ){
			if( ifl[i][j] - ifl[i-1][j] == 1 )
				ibd[i][j] = 1;
			}
		}
	for( i = 1; i <= i21; i++ ){
		for( j = 1; j <= (j21 - 1); j++ ){
			if( ifl[i][j+1] - ifl[i][j] == -1 )
				ibd[i][j] = 1;
			}
		for( j = j21; j >= 2; j-- ){
			if( ifl[i][j] - ifl[i][j-1] == 1 )
				ibd[i][j] = 1;
			}
		}

	for( j = 1; j <= j21; j++ ){
		for( i = 1; i <= i21; i++ ){
			ibd[i][j] = ibd[i][j]*ifl[i][j];
			}
		}
	for( j = 1; j <= j21; j++ ){
		printf( " " );
		for( i = 1; i <= i21; i++ ){
			printf( "%d", ifl[i][j] );
		}
		printf( "  " );
		for( i = 1; i <= i21; i++ ){
			printf( "%d", ibd[i][j] );
		}
		printf( "\n" );
	}
	printf( "OK?\n" );
	scanf( "%d", &ok );

	for( j = 2; j <= (j21 - 1); j++ ){
		for( i = 2; i <= (i21 - 1); i++ ){
			iqq = ibd[i+1][j] + ibd[i-1][j] + ibd[i][j+1] + ibd[i][j-1];
			iqq = iqq*(1 - ifl[i][j]);
			if( iqq == 1 )
				iwt[i][j] = 12;
			if( iqq == 2 )
				iwt[i][j] = 6;
			if( iqq == 3 )
				iwt[i][j] = 4;
			if( iqq == 4 )
				iwt[i][j] = 3;
			}
		}

	for( j = 1; j <= j21; j++ ){
		printf( "   " );
		for( i = 1L; i <= i21; i++ ){
			printf( "%2d", iwt[i][j] );
		}
		printf( "\n" );
	}
	printf( "OK?\n" );
	scanf( "%d", &ok );

	for( j = 1; j <= j21; j++ ){
		for( i = 1; i <= i21; i++ ){
			u[i][j] = 1.;
			v[i][j] = 0.;
			p[i][j] = 0.;
		}
	}
	if( ir == 1 ){
		fp11 = fopen( "file.11", "r" );
		for( j = 1; j <= JD; j++ ){
			for( i = 1; i <= ID; i++ ){
				fscanf( fp11, "%f %f %f ", &u[i][j], &v[i][j], &p[i][j] );
			}
		}
		fclose( fp11 );
	}
	for( l = 1; l <= lm; l++ ){
		for( i = 1; i <= i21; i++ ){
			u[i][j21] = u[i][j19];
			v[i][j21] = -v[i][j19];
			u[i][1] = u[i][3];
			v[i][1] = -v[i][3];
			u[i][j20] = u[i][j19];
			v[i][j20] = 0.;
			u[i][2] = u[i][3];
			v[i][2] = 0.;
		}
		for( j = 1; j <= j21; j++ ){
			u[i21][j] = u[i19][j];
			v[i21][j] = v[i19][j];
			u[i20][j] = u[i19][j];
			v[i20][j] = v[i19][j];
		}
		for( j = 2; j <= j20; j++ ){
			for( i = 2; i <= i20; i++ ){
				u1 = (u[i+1][j] - u[i-1][j])*.5*xg[i];
				u2 = (u[i][j+1] - u[i][j-1])*.5*yg[j];
				v1 = (v[i+1][j] - v[i-1][j])*.5*xg[i];
				v2 = (v[i][j+1] - v[i][j-1])*.5*yg[j];
				q[i][j] = -u1*u1 - v2*v2 - 2.*u2*v1 + td*(u1 + v2);
			}
		}
		for( j = 2; j <= j20; j++ ){
			for( i = 2; i <= i20; i++ ){
				d[i][j] = .5/(xa[i] + ya[j]);
			}
		}
		for( kk = 1; kk <= km; kk++ ){
			g1 = 0.;

			for( j = 2; j <= (j21 - 1); j++ ){
				for( i = 2; i <= (i21 - 1); i++ ){
					s[i][j] = ibd[i+1][j]*p[i+1][j] + 
						  ibd[i-1][j]*p[i-1][j] + ibd[i][j+1]*
						  p[i][j+1] + ibd[i][j-1]*p[i][j-1];
				}
			}
			for( j = 2; j <= (j21 - 1); j++ ){
				for( i = 2; i <= (i21 - 1); i++ ){
					p[i][j] = s[i][j]*(1 - ifl[i][j])*
					 iwt[i][j]/12. + p[i][j]*ifl[i][j];
				}
			}

			for( j = 1; j <= j21; j++ ){
				p[2][j] = p[3][j];
				p[i20][j] = p[i19][j];
			}
			for( i = 1; i <= i21; i++ ){
				p[i][2] = p[i][3];
				p[i][j20] = p[i][j19];
			}

			for( j = 3; j <= j19; j++ ){
				for( i = 3; i <= i19; i++ ){
					ul[i] = d[i][j]*(xc[i]*p[i+1][j] + 
					 yc[j]*p[i][j+1] + xd[i]*p[i-1][j] + 
					 yd[j]*p[i][j-1] - q[i][j]) - 
					 p[i][j];
				}
				for( i = 3; i <= i19; i++ ){
					g1 = g1 + ul[i]*ul[i]*ifl[i][j];
				}
				for( i = 3; i <= i19; i++ ){
					p[i][j] = ul[i] + p[i][j]*ifl[i][j];
				}
			}
			if( g1 <= .0001 )
				goto L_31;
			}
L_31:
		printf( "%12d%12d%12f \n", l, kk, g1 );

		for( j = 3; j <= j19; j++ ){
			for( i = 3; i <= i19; i++ ){
				unx = u[i][j]*(u[i+1][j] - u[i-1][j])/
				      2. - fabs( u[i][j] )*(u[i+1][j] - 2.*
				      u[i][j] + u[i-1][j])/2.;
				uny = v[i][j]*(u[i][j+1] - u[i][j-1])/
				      2. - fabs( v[i][j] )*(u[i][j+1] - 2.*
				      u[i][j] + u[i][j-1])/2.;

				vnx = u[i][j]*(v[i+1][j] - v[i-1][j])/
				      2. - fabs( u[i][j] )*(v[i+1][j] - 2.*
				      v[i][j] + v[i-1][j])/2.;
				vny = v[i][j]*(v[i][j+1] - v[i][j-1])/
				      2. - fabs( v[i][j] )*(v[i][j+1] - 2.*
				      v[i][j] + v[i][j-1])/2.;

				uv = (u[i+1][j] - 2.*u[i][j] + u[i-1][j])*
				     xa[i] + (u[i][j+1] - 2.*u[i][j] + 
				     u[i][j-1])*ya[j] - (u[i+1][j] - 
				     u[i-1][j])*.5*xb[i] - (u[i][j+1] - 
				     u[i][j-1])*.5*yb[j];
				vv = (v[i+1][j] - 2.*v[i][j] + v[i-1][j])*
				     xa[i] + (v[i][j+1] - 2.*v[i][j] + 
				     v[i][j-1])*ya[j] - (v[i+1][j] - 
				     v[i-1][j])*.5*xb[i] - (v[i][j+1] - 
				     v[i][j-1])*.5*yb[j];

				d[i][j] = u[i][j] + dt*(-(unx*xg[i] + 
					  uny*yg[j]) - (p[i+1][j] - p[i-1][j])*
					  .5*xg[i] + r1*uv);
				q[i][j] = v[i][j] + dt*(-(vnx*xg[i] + 
					  vny*yg[j]) - (p[i][j+1] - p[i][j-1])*
					  .5*yg[j] + r1*vv);
			}
		}
		for( j = 3; j <= j19; j++ ){
			for( i = 3; i <= i19; i++ ){
				u[i][j] = d[i][j]*ifl[i][j];
				v[i][j] = q[i][j]*ifl[i][j];
			}
		}
	}
	iwrite = 0;
	if( iwrite == 1 ){
		printf( "---U---\n" );
		for( j = 2; j <= j21; j += 4 ){
			for( i = 2; i <= i21; i += 4 ){
				printf( "%f ", u[i][j] );
			}
		}
		printf( "\n" );
		printf( "---V---\n" );
		for( j = 2; j <= j21; j += 4 ){
			for( i = 2; i <= i21; i += 4 ){
				printf( "%f ", v[i][j] );
			}
		}
		printf( "\n" );
		printf( "---P---\n" );
		for( j = 2; j <= j21; j += 4 ){
			for( i = 2; i <= i21; i += 4 ){
				printf( "%f ", p[i][j] );
				}
			}
		printf( "\n" );
		}
	fp10 = fopen( "file.10", "w" );
	for( j = 1; j <= JD; j++ ){
		for( i = 1; i <= ID; i++ ){
			fprintf( fp10, "%f %f %f \n", u[i][j], v[i][j], p[i][j] );
			}
		}
	fclose( fp10 );
}
