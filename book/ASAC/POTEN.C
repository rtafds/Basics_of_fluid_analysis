#include <stdio.h>
#include <math.h>
#define	JDIM	62
#define	KDIM	31
#define	JDIM1	JDIM+1
#define	KDIM1	KDIM+1
float u[JDIM1][KDIM1], v[JDIM1][KDIM1], p[JDIM1][KDIM1];
float x[JDIM1][KDIM1], y[JDIM1][KDIM1];
float xx[JDIM1][KDIM1], xy[JDIM1][KDIM1], yx[JDIM1][KDIM1], yy[JDIM1][KDIM1], 
      c1[JDIM1][KDIM1], c2[JDIM1][KDIM1], c3[JDIM1][KDIM1], c4[JDIM1][KDIM1], 
      c5[JDIM1][KDIM1], aj[JDIM1][KDIM1];
int jm, jmax, km, kmax, istep0, ityp;
float eps, err, const_, alp;
FILE *fp10,*fp15,*fp12;

main()
{
int i, j, k, kk;

	grid();
	data();
	met1();
	for( kk = 1; kk <= 3; kk++ ){
		init( kk );
		for( i = 1; i <= istep0; i++ ){
			err = 0.;
			bc( kk );
			psi();
			if( err < eps ){
				goto L_60;
			}
			if( i - i/20*20 == 0 ){
				printf( "%d %d %f \n", kk, i, err );
			}
		}
L_60:
		printf( "%f \n", err );
		if( kk == 1 ){
			for( k = 1; k <= kmax; k++ ){
				for( j = 1; j <= jmax; j++ ){
					u[j][k] = p[j][k];
				}
			}
		}
		if( kk == 2 ){
			for( k = 1; k <= kmax; k++ ){
				for( j = 1; j <= jmax; j++ ){
					v[j][k] = p[j][k];
				}
			}
		}
	}
	sup();
	outp();
	exit(0);
}

data()
{
int nsteps;
float alpp, pai;

	jm = jmax - 1;
	km = kmax - 1;
	istep0 = 10;
	nsteps = 100;
	eps = .001;
	const_ = 1.;
	printf("ISTEP0,ALP,EPS ?\n" );
	scanf("%d,%f,%f", &istep0, &alpp, &eps );
	pai = atan( 1. )*4.;
		printf("pai== %e\n",pai);
	alp = pai*alpp/180.;
		printf("alp== %e\n",alp);
	return;
}

grid()
{
int i, ia, ib, igr, iok, j, k;
float aa, bb, bc, hh, pai, ra, rr[KDIM], tt;


	printf("DO YOU WANT TO READ FILE ?  YES(=1),NO(=2)\n" );
	scanf("%d", &igr );

	if( igr == 1 ){

		fp10=fopen( "file.10", "w");
		scanf(fp10, "%d %d", &ia, &ib );
		for( j = 1; j <= ib; j++ ){
			for( i = 1; i <= ia; i++ ){
				scanf( fp10, "%f", &x[i][j] );
			}
		}
		for( j = 1; j <= ib; j++ ){
			for( i = 1; i <= ia; i++ ){
				scanf( fp10, "%f", &y[i][j] );
			}
		}
		jmax = ia;
		kmax = ib;
		ityp = 0;
        	fclose( fp10 );

		}
	else{

		/****   */
		jmax = 62;
		kmax = 31;
		pai = atan( 1. )*4.;
		aa = 1.;
		bb = .3;
L_99:
		printf( "RR=A+A*R+A*R**2+A*R**3+...  INPUT A & R?\n" );
		scanf( "%f,%f", &hh, &ra );
		rr[1] = 1.;
		for( k = 2; k <= kmax; k++ ){
			rr[k] = rr[k - 1] + hh*powi(ra, k - 1);
		}
		for( k = 1; k <= kmax; k++ ){
			printf( "%e ", rr[k] );
			}
		printf( "\n" );
		printf( "OK?  YES=1  NO=2\n" );
		scanf( "%d", &iok );
		if( iok != 1 ){
			goto L_99;
		}
		for( k = 1; k <= kmax; k++ ){
			for( j = 1; j <= jmax; j++ ){
				tt = pai*2.*(float)( j - 2 )/(float)( jmax - 2 );
				bc = bb + (aa - bb)*(float)( k - 1 )/(float)( kmax - 1 );
				x[j][k] = aa*rr[k]*cos( tt );
				y[j][k] = bc*rr[k]*sin( tt );
				}
			}
		ityp = 1;

		fp15=fopen( "file.15", "w");
		fprintf( fp15, "%d %d \n", jmax, kmax );
		for( j = 1; j <= kmax; j++ ){
			for( i = 1; i <= jmax; i++ ){
				fprintf( fp15, "%f ", x[i][j] );
			}
		fprintf( fp15, "\n" );
		}
		for( j = 1; j <= kmax; j++ ){
			for( i = 1; i <= jmax; i++ ){
				fprintf( fp15, "%f ", y[i][j] );
			}
		fprintf( fp15, "\n" );
		}
		fprintf( fp15, "\n" );
		}
        fclose( fp15 );

	return;
}

met1()
{
int j, k;
float ajj, c77, c88, xe, xxi, ye, yxi;

for( k = 1; k <= kmax; k++ ){
	for( j = 1; j <= jmax; j++ ){
		if( k == 1 ){
			xe = 0.5*(-x[j][3] + 4.0*x[j][2] - 3.0*x[j][1]);
			ye = 0.5*(-y[j][3] + 4.0*y[j][2] - 3.0*y[j][1]);
			}
		else if( k == kmax ){
			xe = 0.5*(x[j][kmax-2] - 4.0*x[j][kmax-1] + 3.0*x[j][kmax]);
			ye = 0.5*(y[j][kmax-2] - 4.0*y[j][kmax-1] + 3.0*y[j][kmax]);
			}
		else{
			xe = 0.5*(x[j][k+1] - x[j][k-1]);
			ye = 0.5*(y[j][k+1] - y[j][k-1]);
			}
		if( j == 1 ){
			xxi = 0.5*(-x[3][k] + 4.0*x[2][k] - 3.0*x[1][k]);
			yxi = 0.5*(-y[3][k] + 4.0*y[2][k] - 3.0*y[1][k]);
			if( ityp == 1 ){
				xxi = 0.5*(x[2][k] - x[jmax-2][k]);
				yxi = 0.5*(y[2][k] - y[jmax-2][k]);
				}
			}
		else if( j == jmax ){
			xxi = 0.5*(x[jmax-2][k] - 4.0*x[jmax-1][k] + 3.0*x[jmax][k]);
			yxi = 0.5*(y[jmax-2][k] - 4.0*y[jmax-1][k] + 3.0*y[jmax][k]);
			if( ityp == 1 ){
				xxi = 0.5*(x[3][k] - x[jmax-1][k]);
				yxi = 0.5*(y[3][k] - y[jmax-1][k]);
				}
			}
		else{
			xxi = 0.5*(x[j+1][k] - x[j-1][k]);
			yxi = 0.5*(y[j+1][k] - y[j-1][k]);
			}
		ajj = xxi*ye - xe*yxi;
		xx[j][k] = ye/ajj;
		yx[j][k] = -yxi/ajj;
		xy[j][k] = -xe/ajj;
		yy[j][k] = xxi/ajj;
		aj[j][k] = ajj;
		}
		}
	for( k = 1; k <= kmax; k++ ){
		for( j = 1; j <= jmax; j++ ){
			c1[j][k] = xx[j][k]*xx[j][k]+xy[j][k]*xy[j][k];
			c3[j][k] = yx[j][k]*yx[j][k]+yy[j][k]*yy[j][k];
			c2[j][k] = 2.*(xx[j][k]*yx[j][k]+xy[j][k]*yy[j][k]);
			}
		}
	for( k = 2; k <= km; k++ ){
		for( j = 2; j <= jm; j++ ){
			c77 = xx[j][k]*(xx[j+1][k] - xx[j-1][k]) 
                            + yx[j][k]*(xx[j][k+1] - xx[j][k-1]) 
                            + xy[j][k]*(xy[j+1][k] - xy[j-1][k]) 
			    + yy[j][k]*(xy[j][k+1] - xy[j][k-1]);
			c88 = xx[j][k]*(yx[j+1][k] - yx[j-1][k])
			    + yx[j][k]*(yx[j][k+1] - yx[j][k-1])
			    + xy[j][k]*(yy[j+1][k] - yy[j-1][k])
			    + yy[j][k]*(yy[j][k+1] - yy[j][k-1]);
			c4[j][k] = c77*.5;
			c5[j][k] = c88*.5;
			}
		}
	return;
}

init(kk)
int kk;
{
int j, k;

	if( kk == 1 ){
		for( k = 1; k <= kmax; k++ ){
			for( j = 1; j <= jmax; j++ ){
				p[j][k] = y[j][k];
				}
			}
		}
	if( kk == 2 ){
		for( k = 1; k <= kmax; k++ ){
			for( j = 1; j <= jmax; j++ ){
				p[j][k] = -x[j][k];
				}
			}
		}
	if( kk == 3 ){
		for( k = 1; k <= kmax; k++ ){
			for( j = 1; j <= jmax; j++ ){
				p[j][k] = 0.;
				}
			}
		}
	return;
}

psi()
{
int j, k;
float cc, pa, pp;

 for( k = 2; k <= km; k++ ){
  for( j = 2; j <= jm; j++ ){
	cc = .5/(c1[j][k]+c3[j][k]);
	pa = c1[j][k]*(p[j+1][k]+p[j-1][k])+c3[j][k]*(p[j][k+1]+p[j][k-1])
           + .25*c2[j][k]*(p[j+1][k+1]-p[j-1][k+1]-p[j+1][k-1]+p[j-1][k-1])
           + .5*c4[j][k]*(p[j+1][k]-p[j-1][k]) 
           + .5*c5[j][k]*(p[j][k+1]-p[j][k-1]);
	pp = pa*cc;
	err = err+(pp - p[j][k])*(pp-p[j][k]);
	p[j][k]=p[j][k]*(1.-const_)+pp*const_;
  }
 }
	return;
}

bc(kk)
int kk;
{
int j, k;
float ot[JDIM1], vt[JDIM1];


	if( kk == 1 ){
		for( j = 1; j <= jmax; j++ ){
			p[j][1] = 0.;
			p[j][kmax] = y[j][kmax];
			}
		}

	if( kk == 2 ){
		for( j = 1; j <= jmax; j++ ){
			p[j][1] = 0.;
			p[j][kmax] = -x[j][kmax];
			}
		}

	if( kk == 3 ){
		for( j = 1; j <= jmax; j++ ){
			p[j][1] = 1.;
			p[j][kmax] = 0.;
			}
		}

	for( k = 1; k <= kmax; k++ ){
		p[1][k] = p[jm][k];
		p[jmax][k] = p[2][k];
		}

	return;
}

sup()
{
int j, k;
float aa, bb, cc, g1, g2;

	cc = 0.;
	g1 = sqrt( xy[3][2]*xy[3][2] + yy[3][2]*yy[3][2] );
	g2 = sqrt( xy[jm-1][2]*xy[jm-1][2] + yy[jm-1][2]*yy[jm-1][2] );
	aa = (g1*u[3][2] - g2*u[jm-2][2])*cos( alp ) 
           + (g1*v[3][2] - g2*v[jm-2][2])*sin( alp );
	bb = g2*p[jm-2][2] - g1*p[3][2];
	if( bb != 0. )
		cc = aa/bb;
	printf( "--- LAMDA = %f \n", cc );
	for( k = 1; k <= kmax; k++ ){
		for( j = 1; j <= jmax; j++ ){
			p[j][k] = u[j][k]*cos( alp ) + v[j][k]*
			 sin( alp ) + cc*p[j][k];
			}
		}
	return;
}

outp()
{
int j, k;

	fp12=fopen( "file.12", "w");
	fprintf( fp12, "%d %d \n", jmax, kmax );
	for( k = 1; k <= kmax; k++ ){
		for( j = 1; j <= jmax; j++ ){
			fprintf( fp12, "%f ", p[j][k] );
			}
		}
	fprintf( fp12, "\n" );
        fclose( fp12 );
	return;
}
