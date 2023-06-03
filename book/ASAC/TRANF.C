#include <stdio.h>
#include <math.h>

#define MX 11
#define MY  9

float	x[MX][MY], y[MX][MY];

main()
{
	int	i, j, ic;
	float	pai, tet, te, tf, a, b, fct, x11, y11, xa, ya;

	pai = 4. * atan( 1. );
	tet = pai / 18.;

	for ( i=1; i<=MX; i++ ){
		te = tet * ( i + 17 );
		x[i][1] = cos( te );
		y[i][1] = sin( te );
		tf = tet * ( 11 - i );
		x[i][MY] = cos( tf );
		y[i][MY] = sin( tf );
	}

	for ( j=1; j<=MY; j++ ){
		te = tet * ( 19 - j );
		x[1][j] = cos( te );
		y[1][j] = sin( te );
		tf = tet * ( j + 27 );
		x[MX][j] = cos( tf );
		y[MX][j] = sin( tf );
	}

	for ( j=2; j<=MY-1; j++ ){
		for ( i=2; i<=MX-1; i++ ){
			a=(float)(i-MX)/(float)(1-MX);
			b=(float)(j-MY)/(float)(1-MY);
			x[i][j]=a*x[1][j]+(1-a)*x[MX][j]+b*x[i][1]+(1-b)*x[i][MY]-a*b*x[1][1]
			      -a*(1-b)*x[1][MY]-(1-a)*b*x[MX][1]-(1-a)*(1-b)*x[MX][MY];
			y[i][j]=a*y[1][j]+(1-a)*y[MX][j]+b*y[i][1]+(1-b)*y[i][MY]-a*b*y[1][1]
			      -a*(1-b)*y[1][MY]-(1-a)*b*y[MX][1]-(1-a)*(1-b)*y[MX][MY];
		}
	}

/*	plots();			*/
	printf( "fct,x11,y11\n" );
	scanf( "%f%f%f", &fct, &x11, &y11 );
/*	factor( fct );
	plot( 0., 0., -3 );		*/

	for ( j=1; j<=MY; j++ ){
		for ( i=1; i<=MX; i++ ){
			ic = 2;
			if ( i == 1 )	ic = 3;
			xa = x[i][j] + x11;
			ya = y[i][j] + y11;
/*			plot( xa, ya, ic );	*/
		}
	}

	for ( i=1; i<=MX; i++ ){
		for ( j=1; j<=MY; j++ ){
			ic = 2;
			if ( j == 1 )	ic = 3;
			xa = x[i][j] + x11;
			ya = y[i][j] + y11;
/*			plot( xa, ya, ic );	*/
		}
	}

/*	plot( 0., 0., -3 );
	plot( 0., 0., 999 );			*/

}
