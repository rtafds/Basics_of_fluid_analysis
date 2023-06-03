#include <stdio.h>
#include <math.h>

#define MX 42
#define MY 16

float x[MX][MY],y[MX][MY],p[MX][MY],q[MX][MY];

main()
{
	int	ity,k,i,j,ic,i1,k1;
	float	r1,r2,x3,pai,r,tet,fct,x11,y11,dp,x1,y1;
	FILE	*fp11, *fp10;

		printf( "r1,r2,x3\n" );
		scanf( "%f%f%f", &r1,&r2,&x3 );
		pai = 4. * atan( 1. );

		for ( k=1; k<=MY-1; k++ ){
			for ( i=1; i<=MX; i++ ){
				r=r1+(r2-r1-x3-.1)*(float)(k-1)/(float)(MY-1);
				tet=pai*2.*(i-1)/(float)(MX-2);
				x[i][k]=r*cos(tet);
				y[i][k]=.3*r*sin(tet);
			}
		}

		for ( i=1; i<=MX; i++ ){
			tet=pai*2.*(i-1)/(float)(MX-2);
			x[i][MY]=r2*cos(tet)-x3;
			y[i][MY]=r2*sin(tet);
		}
		ellip2();

		fp11 = fopen( "file.11", "w" );
		for ( k=1; k<=MY; k++ ){
			for ( j=1; j<=MX; j++ ){
				fprintf( fp11, "%e  %e\n", x[j][k], y[j][k] );
			}
		}
		fclose( fp11 );

	}


ellip2()
{
	int	kd, ld, k100, ityp, l, k, k4, l4, l3, k30, k60, kk;
	float	aa, cc, sg, pk, qk, const1, xxi, yxi;
	float	xeta, yeta, alph, beta, gama, ajac, wxij, qeta, ag, sx, sy, gosa, eps1;

	kd = MX;
	ld = MY;

	printf( "k100,ityp\n" );
	scanf( "%d%d", &k100, &ityp );

	if ( ityp == 0 ){
		for ( l=1; l<=ld; l++ ){
			for ( k=1; k<=kd; k++ ){
				p[k][l] = 0.;
				q[k][l] = 0.;
			}
		}
	}
	else if ( ityp == 1 ){
		printf( "aa,cc,k4,l4\n" );
		scanf( "%f%f%f%f", &aa, &cc, &k4, &l4 );

		for ( l=1; l<=ld; l++ ){
			for ( k=1; k<=kd; k++ ){
				sg = 1;
				if ( k == k4 )	sg = -1;

				p[k][l]=-aa*sg*exp(-cc*sqrt(pow(((float)(k-k4)),2)
					+pow(((float)(l-l4)),2)));
				q[k][l]=0.;
			}
		}
	}
	else{
		printf( "aa,cc,l3\n" );
		scanf( "%f%f%f", &aa, &cc, &l3 );

		for ( l=1; l<=ld; l++ ){
			p[k][l] = 0.;
			sg = 1.;
			if ( l > l3 )	sg = -1.;
			q[k][l]=-aa*sg*exp(-cc*abs(l-l3));
		}
	}

	k30 = k100 / 3;
	k60 = k30 * 2;
	eps1 = 0.000001;
	pk = 0.;
	qk = 0.;
	const1 = .05;

	for ( kk=1; kk<=k100; kk++ ){
		gosa = 0.0;
		if ( const1 == 1. && kk > k30 )	const1 = .25;
		if ( const1 == 1. && kk > k60 )	const1 = 1.;
		if ( ityp == 0 )	const1 = 1.;

		for ( l=2; l<=ld-1; l++ ){
			for ( k=2; k<=kd-1; k++ ){
				xxi=(x[k+1][l]-x[k-1][l])*.5;
				yxi=(y[k+1][l]-y[k-1][l])*.5;
				xeta=(x[k][l+1]-x[k][l-1])*.5;
				yeta=(y[k][l+1]-y[k][l-1])*.5;
				alph=xeta*xeta+yeta*yeta;
				beta=xxi*xeta+yxi*yeta;
				gama=xxi*xxi+yxi*yxi;
				ajac=xxi*yeta-xeta*yxi;
				wxij=p[k][l]*alph/(ajac*ajac);
				qeta=q[k][l]*gama/(ajac*ajac);
				ag=0.5/(alph+gama);
				sx=(alph*(x[k+1][l]+x[k-1][l])-.5*beta*(x[k+1][l+1]
				  -x[k-1][l+1]-x[k+1][l-1]+x[k-1][l-1])+gama*(x[k][l+1]
				  +x[k][l-1])+ajac*ajac*(wxij*xxi+qeta*xeta))*ag-x[k][l];
				sy=(alph*(y[k+1][l]+y[k-1][l])-.5*beta*(y[k+1][l+1]
				  -y[k-1][l+1]-y[k+1][l-1]+y[k-1][l-1])+gama*(y[k][l+1]
				  +y[k][l-1])+ajac*ajac*(wxij*yxi+qeta*yeta))*ag-y[k][l];
				x[k][l]=x[k][l]+const1*sx;
				y[k][l]=y[k][l]+const1*sy;
				gosa=gosa+sx*sx+sy*sy;
			}
		}

		for ( l=1; l<=ld; l++ ){
			x[1][l]=x[kd-1][l];
			y[1][l]=y[kd-1][l];
			x[kd][l]=x[2][l];
			y[kd][l]=y[2][l];
		}

		if ( gosa < eps1 )	goto r122;
		printf( "%12d  %e\n", kk, gosa );
	}

	r122:
	printf( "%12d  %e\n", kk, gosa );
}
