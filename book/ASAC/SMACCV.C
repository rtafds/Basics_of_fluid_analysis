#include <stdio.h>
#include <math.h>

#define	MD	21
#define	ND	21
#define MD1	MD+1
#define ND1	ND+1

	float	u[MD1][ND1], v[MD1][ND1], p[MD1][ND1], q[MD1][ND1], d[MD1][ND1],
r[MD1][ND1],
		ul[MD1], h[MD1];

main()
{
	int	i, i19, i20, i21,j, j19, j20, j21, k, km, l, lm;
	float	a1, a2, a3, afct, div2, divv, dt, dx, dy, g1, g2, r1, re, td,
		u1, ua, ub, uli, un, uv, v1, v2, va, vb, vn, vv, xd, xdh, yd,
		ydh;
	FILE	*fp17;

	printf( "DT,RE,LM,AFCT (0.01,40,300,0.5) \n" );
	scanf( "%f,%f,%d,%f", &dt, &re, &lm, &afct );
	i21 = MD;
	j21 = ND;
	dx = 1./(float)( i21 - 3 );
	dy = 1./(float)( j21 - 3 );
	xd = 1./dx;
	yd = 1./dy;
	xdh = .5*xd;
	ydh = .5*yd;
	td = 1./dt;
	r1 = 1./re;
	i20 = i21 - 1;
	i19 = i21 - 2;
	j20 = j21 - 1;
	j19 = j21 - 2;
	km = 20;
	printf( "KM? (50) \n" );
	scanf( "%d", &km );
	a1 = .5*dy*dy/(dx*dx + dy*dy);
	a2 = .5*dx*dx/(dx*dx + dy*dy);
	a3 = .5*dy*dy/(1. + dy*dy/(dx*dx));
	for( j = 1; j <= j21; j++ ){
		for( i = 1; i <= i21; i++ ){
			u[i][j] = 0.;
			v[i][j] = 0.;
			p[i][j] = 0.;
			q[i][j] = 0.;
		}
	}

	for( l = 1; l <= lm; l++ ){
		for( i = 1; i <= i21; i++ ){
			u[i][j20] = 1.;
			v[i][j21] = v[i][j19];
			v[i][j20] = 0.;
			v[i][1] = v[i][3];
			v[i][2] = 0.;
			u[i][1] = -u[i][2];
		}
		for( j = 1; j <= j21; j++ ){
			u[i21][j] = u[i19][j];
			u[i20][j] = 0.;
			v[i20][j] = -v[i19][j];
			u[1][j] = u[3][j];
			u[2][j] = 0.;
			v[1][j] = -v[2][j];
		}

		divv = 0.;
		for( j = 2; j <= j19; j++ ){
			for( i = 2; i <= i19; i++ ){
				u1 = (u[i+1][j] - u[i][j])*xd;
				v2 = (v[i][j+1] - v[i][j])*yd;
				d[i][j] = u1 + v2;
				divv = divv + fabs( u1 + v2 );
				ua = .25*(u[i][j] + u[i+1][j] + u[i+1][j+1] + u[i][j+1]);
				ub = .25*(u[i][j] + u[i+1][j] + u[i+1][j-1] + u[i][j-1]);
				va = .25*(v[i][j] + v[i][j+1] + v[i+1][j+1] + v[i+1][j]);
				vb = .25*(v[i][j] + v[i][j+1] + v[i-1][j+1] + v[i-1][j]);
				r[i][j] = -u1*u1 - 2.*(ua - ub)*(va - vb)*xd*yd
					  - v2*v2 + td*(u1 + v2);
			}
		}
		for( k = 1; k <= km; k++ ){
			g2 = 0.;
			for( j = 1; j <= j21; j++ ){
				p[1][j] = p[2][j];
				p[i20][j] = p[i19][j];
			}
			for( i = 1; i <= i21; i++ ){
				p[i][1] = p[i][2];
				p[i][j20] = p[i][j19];
			}
			for( j = 2; j <= j19; j++ ){
				for( i = 2; i <= i19; i++ ){
					uli = a1*(p[i+1][j] + p[i-1][j]) +
					      a2*(p[i][j+1] + p[i][j-1]) -
					      a3*r[i][j] - p[i][j];
					g2 = g2 + uli*uli;
					p[i][j] = uli + p[i][j];
				}
			}
			if( g2 <= .00001 )
				goto L_331;
			}
L_331:
		;
		for( k = 1; k <= km; k++ ){
			g1 = 0.;
			for( j = 1; j <= j21; j++ ){
				q[1][j] = q[2][j];
				q[i20][j] = q[i19][j];
			}
			for( i = 1; i <= i21; i++ ){
				q[i][1] = q[i][2];
				q[i][j20] = q[i][j19];
			}
			for( j = 2; j <= j19; j++ ){
				for( i = 2; i <= i19; i++ ){
					uli = a1*(q[i+1][j] + q[i-1][j]) +
					      a2*(q[i][j+1] + q[i][j-1]) +
					      a3*d[i][j] - q[i][j];
					g1 = g1 + uli*uli;
					q[i][j] = uli + q[i][j];
					}
				}
			if( g1 <= .000001 )
				goto L_31;
			}
L_31:
		;
		for( j = 2; j <= j19; j++ ){
			for( i = 3; i <= i19; i++ ){
				u[i][j] = u[i][j] + (q[i][j] - q[i-1][j])*xd*afct;
			}
		}
		for( j = 3; j <= j19; j++ ){
			for( i = 2; i <= i19; i++ ){
				v[i][j] = v[i][j] + (q[i][j] - q[i][j-1])*yd*afct;
			}
		}
		div2 = 0.;
		for( j = 3; j <= j19; j++ ){
			for( i = 3; i <= i19; i++ ){
				u1 = (u[i+1][j] - u[i][j])*xd;
				v2 = (v[i][j+1] - v[i][j])*yd;
				div2 = div2 + fabs( u1 + v2 );
			}
		}
		if( (l%20) == 0 ){
			printf( "%d %d %g %g %g \n", l, k, g1, div2, u[3][11] );
		}
		for( j = 2; j <= j19; j++ ){
			for( i = 2; i <= i19; i++ ){
				;
			}
		}
		for( j = 2; j <= j19; j++ ){
			p[1][j] = p[2][j];
			p[i20][j] = p[i19][j];
		}
		for( i = 2; i <= i19; i++ ){
			p[i][1] = p[i][2];
			p[i][j20] = p[i][j19];
		}
		for( j = 2; j <= j19; j++ ){
			for( i = 3; i <= i19; i++ ){
				v1 = .25*(v[i][j] + v[i][j+1] + v[i-1][j+1] + v[i-1][j]);
				un = u[i][j]*(u[i+1][j] - u[i-1][j])*
				     xdh + v1*(u[i][j+1] - u[i][j-1])*ydh;
				uv = (u[i+1][j] - 2.*u[i][j] + u[i-1][j])*
				     xd*xd + (u[i][j+1] - 2.*u[i][j] + u[i][j-1])*yd*yd;

				u[i][j] = u[i][j] + dt*(-un - (p[i][j] - 
				          p[i-1][j])*xd + r1*uv);
			}
		}
		for( j = 3; j <= j19; j++ ){
			for( i = 2; i <= i19; i++ ){
				u1 = .25*(u[i][j] + u[i+1][j] + u[i+1][j-1] + u[i][j-1]);
				vn = u1*(v[i+1][j] - v[i-1][j])*xdh + 
				     v[i][j]*(v[i][j+1] - v[i][j-1])*ydh;
				vv = (v[i+1][j] - 2.*v[i][j] + v[i-1][j])*
				     xd*xd + (v[i][j+1] - 2.*v[i][j] + v[i][j-1])*yd*yd;

				v[i][j] = v[i][j] + dt*(-vn - (p[i][j] - 
					  p[i][j-1])*yd + r1*vv);
			}
		}
	}
	printf( "---U---\n" );
	for( j = 2; j <= 20; j += 2 ){
		printf( " " );	
		for( i = 2; i <= 20; i += 2 ){
			printf( "%6.3f", u[i][j] );
		}
		printf( "\n" );	
	}
	printf( "---V---\n" );
	for( j = 2; j <= 20; j += 2 ){
		printf( " " );	
		for( i = 2; i <= 20; i += 2 ){
			printf( "%6.3f", v[i][j] );
			}
		printf( "\n" );	
		}
	printf( "---P---\n" );
	for( j = 2; j <= 20; j += 2 ){
		printf( " " );	
		for( i = 2; i <= 20; i += 2 ){
			printf( "%6.3f", p[i][j] );
			}
		printf( "\n" );	
		}
	for( i = 1; i <= i21; i++ ){
		q[i][1] = 0.;
		}
	for( j = 2; j <= j21; j++ ){
		for( i = 1; i <= i21; i++ ){
			q[i][j] = q[i][j-1] + dy*u[i][j];
			}
		}
	out( i21, j21 );
	fp17 = fopen( "file.17", "w" );
	for( j = 1; j <= j21; j++ ){
		for( i = 1; i <= i21; i++ ){
			fprintf( fp17, "%f %f %f ", u[i][j], v[i][j], q[i][j] );
			}
		}
	fclose( fp17 );
	exit(0);
}

out( i21, j21)
int i21, j21;
{
	int	i, ind[100], j;
	float	qma, qmi;

	qma = q[1][1];
	qmi = q[1][1];
	for( j = 1; j <= j21; j++ ){
		for( i = 1; i <= i21; i++ ){
			if( qma < q[i][j] )
				qma = q[i][j];
			if( qmi > q[i][j] )
				qmi = q[i][j];
			}
		}
        for( j=j21-1; j>=2; j+=-1 ){
                for( i=2; i<=i21-1; i++ ){
                        ind[i]=(int)((q[i][j]-qmi)/(qma-qmi)*9.9999);
                        ind[i]=ind[i]*11;
                        }
                for( i=2; i<=i21-1; i++ ){
                        printf("%2d",ind[i]);
                        }
                printf("\n");
        }
}
