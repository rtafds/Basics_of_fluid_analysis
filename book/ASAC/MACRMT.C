#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define	MD	21
#define	ND	21
#define MD1     MD+1
#define ND1     ND+1


float	u[MD1][ND1], v[MD1][ND1], p[MD1][ND1], q[MD1][ND1], t[MD1][ND1],
r[MD1][ND1],
	ul[MD1], h[MD1], x[MD1], y[ND1], a1[MD1], b1[MD1], c1[MD1], a2[MD1],
	b2[MD1], c2[MD1], a3[ND1], b3[ND1], c3[ND1], a4[ND1], b4[ND1], c4[ND1];

main()
{
	int	i, i19, i20, i21, j, j19, j20, j21, k, km, l, lm, ja, jb;
	float	ar, bb, bx, divv, dt, dy, ff, g2, prn, px, py, r1, r2, re, td,
		tn, tv, ty, u1, u2, uli, un, uv, v1, v2, vn, vv, x1, x2,
		x3, x4, y1, y2, y3, y4;
	FILE	*fp17;

	printf( "DT,RE,AR,LM,BB (0.01,40,20,400,0.98) \n" );
	scanf( "%f,%f,%f,%d,%f", &dt, &re, &ar, &lm, &bb );
	printf( "JA,JB (17,5) \n");
	scanf( "%d,%d", &ja, &jb );
	i21 = MD;
	j21 = ND;
	prn = 0.71;
	dy = 1./(float)( j21 - 3 );
	td = 1./dt;
	r1 = 1./re;
	r2 = 1./(re*prn);

	/****  GRID */
	ff = (exp( bb ) + 1.)/(exp( bb ) - 1.);
	for( i = 1; i <= i21; i++ ){
		bx = bb*(float)( i - 1 )/(float)( i21 - 1 );
		x[i] = ff*(exp( bx ) - 1.)/(exp( bx ) + 1.);
	}
	for( j = 1; j <= j21; j++ ){
		bx = bb*(float)( j - 1 )/(float)( j21 - 1 );
		y[j] = ff*(exp( bx ) - 1.)/(exp( bx ) + 1.);
	}
	for( i = 1; i <= i21; i++ ){
		printf( "%f ", x[i] );
	}
	printf( "\n" );
	for( j = 1; j <= j21; j++ ){
		printf( "%f ", y[j] );
	}
	printf( "\n" );
	for( i = 2; i <= (i21 - 1); i++ ){
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
	for( j = 2; j <= (j21 - 1); j++ ){
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
	i20 = i21 - 1;
	i19 = i21 - 2;
	j20 = j21 - 1;
	j19 = j21 - 2;
	km = 20;
	printf( "KM? (20) \n" );
	scanf( "%d", &km );
	for( j = 1; j <= j21; j++ ){
		for( i = 1; i <= i21; i++ ){
			u[i][j] = 0.;
			v[i][j] = 0.;
			p[i][j] = 0.;
			t[i][j] = 0.;
		}
	}

	for( l = 1; l <= lm; l++ ){
		for( i = 1; i <= i21; i++ ){
			u[i][j21] = 1.;
			v[i][j21] = 0.;
			v[i][1] = 0.;
			u[i][1] = 0.;
			t[i][j21] = 0.;
			t[i][1] = 1.;
		}
		for( j = jb; j <= j21; j++ ){
			u[i21][j] = 0.;
			u[i21][j] = 0.;
			t[i21][j] = t[i20][j];
		}
		for( j = 1; j <= ja; j++ ){
			u[1][j] = 0.;
			v[1][j] = 0.;
			t[1][j] = t[2][j];
		}
		for( j = ja+1; j <= j21; j++ ){
			u[1][j] = 1.;
			v[1][j] = 0.;
			t[1][j] = t[2][j];
		}
		for( j = 1; j<=jb; j++ ){
			u[i21][j] = u[i20][j]*0.+1.;
			v[i21][j] = v[i20][j]*0.;
			t[i21][j] = t[i20][j];
		}

		divv = 0.;
		for( j = 2; j <= j20; j++ ){
			for( i = 2; i <= i20; i++ ){
				u1 = a1[i]*u[i-1][j] + b1[i]*u[i][j] +
				     c1[i]*u[i+1][j];
				u2 = a3[j]*u[i][j-1] + b3[j]*u[i][j] +
				     c3[j]*u[i][j+1];
				v1 = a1[i]*v[i-1][j] + b1[i]*v[i][j] +
				     c1[i]*v[i+1][j];
				v2 = a3[j]*v[i][j-1] + b3[j]*v[i][j] +
				     c3[j]*v[i][j+1];
				ty = a3[j]*t[i][j-1] + b3[j]*t[i][j] +
				     c3[j]*t[i][j+1];
				r[i][j] = -u1*u1 - 2.*u2*v1 - v2*v2 + td*(u1 +
				     v2) - ar*ty;
			}
		}

		for( k = 1; k <= km; k++ ){
			g2 = 0.;
			for( j = 1; j <= j21; j++ ){
				p[1][j] = p[2][j];
				p[i21][j] = p[i20][j];
			}
			for( i = 1; i <= i21; i++ ){
				p[i][1] = p[i][2];
				p[i][j21] = p[i][j20];
			}
			for( j = 2; j <= j20; j++ ){
				for( i = 2; i <= i20; i++ ){
					uli = -(a2[i]*p[i-1][j] + c2[i]*
					      p[i+1][j] + a4[j]*p[i][j-1] +
					      c4[j]*p[i][j+1] - r[i][j])/(b2[i] +
					      b4[j]) - p[i][j];
					g2 = g2 + uli*uli;
					p[i][j] = uli + p[i][j];
				}
			}
			if( g2 <= .000001 )
				goto L_331;
		}
L_331:
		;
		if( (l%20) == 0 ){
			printf( "%d %d %f \n", l, k, g2 );
		}

		for( j = 2; j <= j20; j++ ){
			for( i = 2; i <= i20; i++ ){

				un = u[i][j]*(a1[i]*u[i-1][j] + b1[i]*
				     u[i][j] + c1[i]*u[i+1][j]) + v[i][j]*
				     (a3[j]*u[i][j-1] + b3[j]*u[i][j] +
				     c3[j]*u[i][j+1]);
				uv = a2[i]*u[i-1][j] + b2[i]*u[i][j] +
				     c2[i]*u[i+1][j] + a4[j]*u[i][j-1] +
				     b4[j]*u[i][j] + c4[j]*u[i][j+1];
				px = a1[i]*p[i-1][j] + b1[i]*p[i][j] +
				     c1[i]*p[i+1][j];

				u[i][j] = u[i][j] + dt*(-un - px + r1*uv);
			}
		}
		for( j = 2; j <= j20; j++ ){
			for( i = 2; i <= i20; i++ ){
				vn = u[i][j]*(a1[i]*v[i-1][j] + b1[i]*
				     v[i][j] + c1[i]*v[i+1][j]) + v[i][j]*
				     (a3[j]*v[i][j-1] + b3[j]*v[i][j] +
				     c3[j]*v[i][j+1]);
				vv = a2[i]*v[i-1][j] + b2[i]*v[i][j] +
				     c2[i]*v[i+1][j] + a4[j]*v[i][j-1] +
				     b4[j]*v[i][j] + c4[j]*v[i][j+1];
				py = a3[j]*p[i][j-1] + b3[j]*p[i][j] +
				     c3[j]*p[i][j+1];

				v[i][j] = v[i][j] + dt*(-vn - py + r1*vv -
				     ar*t[i][j]);
			}
		}
		for( j = 2; j <= j20; j++ ){
			for( i = 2; i <= i20; i++ ){
				tn = u[i][j]*(a1[i]*t[i-1][j] + b1[i]*
				     t[i][j] + c1[i]*t[i+1][j]) + v[i][j]*
				     (a3[j]*t[i][j-1] + b3[j]*t[i][j] +
				     c3[j]*t[i][j+1]);
				tv = a2[i]*t[i-1][j] + b2[i]*t[i][j] +
				     c2[i]*t[i+1][j] + a4[j]*t[i][j-1] +
				     b4[j]*t[i][j] + c4[j]*t[i][j+1];

				t[i][j] = t[i][j] + dt*(-tn + r2*tv);
			}
		}
	}
	printf( "---U---\n" );
	printf( " " );
	for( j = 2; j <= 20; j += 2 ){
		for( i = 2; i <= 20; i += 2 ){
			printf( "%6.3f", u[i][j] );
		}
	}
	printf( "\n" );
	printf( "---V---\n" );
	printf( " " );
	for( j = 2; j <= 20; j += 2 ){
		for( i = 2; i <= 20; i += 2 ){
			printf( "%6.3f", v[i][j] );
		}
	}
	printf( "\n" );
	printf( "---P---\n" );
	printf( " " );
	for( j = 2; j <= 20; j += 2 ){
		for( i = 2; i <= 20; i += 2 ){
			printf( "%6.3f", p[i][j] );
		}
	}
	printf( "\n" );
	for( i = 1; i <= i21; i++ ){
		q[i][1] = 0.;
	}
	for( j = 2; j <= j21; j++ ){
		for( i = 1; i <= i21; i++ ){
			q[i][j] = q[i][j-1] + dy*u[i][j];
		}
	}
	out( t, &i21, &j21 );
	printf( "\n" );
	out( q, &i21, &j21 );
	fp17 = fopen( "file.17", "w" );
	for( j = 1; j <= j21; j++ ){
		for( i = 1; i <= i21; i++ ){
			fprintf( fp17, "%f %f %f ", u[i][j], v[i][j], q[i][j] );
		}
	}
	fclose( fp17 );
	exit(0);
}

out(q, i21, j21)
float	q[MD1][ND1];
int	*i21, *j21;
{
	int	i, ind[100], j;
	float	qma, qmi;

	qma = q[1][1];
	qmi = q[1][1];
	for( j = 1; j <= *j21; j++ ){
		for( i = 1; i <= *i21; i++ ){
			if( qma < q[i][j] )
				qma = q[i][j];
			if( qmi > q[i][j] )
				qmi = q[i][j];
		}
	}
	for( j = *j21 - 1; j >= 2; j-- ){
		for( i = 2; i <= (*i21 - 1); i++ ){
			ind[i] = (q[i][j] - qmi)/(qma - qmi)*9.9999;
			ind[i] = ind[i]*11;
		}
		printf( " " );
		for( i = 2; i <= (*i21 - 1); i++ ){
			printf( "%2d ", ind[i] );
		}
		printf( "\n" );
	}
	return;
}
