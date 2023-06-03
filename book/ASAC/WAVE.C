#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/***********************************************************************
 *     WAVE EQUATION    EXPLICIT METHOD                                *
 *********************************************************************** */
#define	NX	51
#define NX1     NX+1

main()
{
	int	i, i5, k, km, mx;
	float	dt, dx, pai, r, u[NX1], uu[NX1], v[NX1], vv[NX1], x;

	/***** INPUT & CALCULATE PARAMETERS
	 * */
	printf( "MESH POINTS(<52) ?  NUMBER OF TIME STEP ?? (20,1000) \n" );
	scanf( "%d,%d", &mx, &km );
L_99:
	printf( "DELTA T: DT ? (0.001) \n" );
	scanf( "%f", &dt );
	dx = 1./(float)( mx - 1 );
	r = dt/dx;
	if( r > 1 )
		goto L_99;
	i5 = .125/dt;
	if( i5 == 0 )
		i5 = 1;
	pai = atan( 1. )*4.;
	output( u, mx, 1 );

	/***** INITIAL CONDITION */
	for( i = 1; i <= mx; i++ ){
		x = (float)( i - 1 )/(float)( mx - 1 );
		u[i] = 0.25*sin( pai*x );
		v[i] = 0.25*sin( pai*x );
	}

	for( k = 1; k <= km; k++ ){
		u[1] = 0.;
		u[mx] = 0.;

		for( i = 1; i <= mx; i++ ){
			vv[i] = u[i] + .3;
		}

		if( (k%i5) == 1 )
			output( vv, mx, 2 );

		for( i = 2; i <= (mx - 1); i++ ){
			uu[i] = r*r*(u[i-1] - 2.*u[i] + u[i+1]) + 2.*
			 u[i] - v[i];
		}

		for( i = 2; i <= (mx - 1); i++ ){
			v[i] = u[i];
			u[i] = uu[i];
		}

	}

	output( vv, mx, 3 );
	exit(0);
}
output(u, mx, mm)
float u[NX1];
int mx, mm;
{
        int     z[61][24];
        int     i, j;

        if( mm == 1 ){
                for( j = 1; j <= 23; j++ ){
                        for( i = 2; i <= (mx - 1); i++ ){
                                z[i][j] = 10;
                        }
                }
                for( j = 1; j <= 23; j++ ){
                        z[1][j] = 11;
                        z[mx][j] = 11;
                }
                for( i = 1; i <= mx; i++ ){
                        z[i][1] = 12;
                        z[i][23] = 12;
                }
        }

        if( mm == 2 ){
                for( i = 1; i <= mx; i++ ){
                        j = u[i]*40 + .001;
                        if( j > 23 )
                                goto L_40;
                        z[i][j] = 13;
L_40:
                        ;
                }
        }

        if( mm == 3 ){
                for( j = 23; j >= 1; j-- ){
                        printf( " " );
                        for( i = 1; i <= 60; i++ ){
                                if ( z[i][j] == 10 )
                                        printf( " " );
                                else if ( z[i][j] == 11 )
                                        printf( ":" );
                                else if ( z[i][j] == 12 )
                                        printf( "-" );
                                else if ( z[i][j] == 13 )
                                        printf( "*" );
                        }
                        printf( "\n" );
                }
        }
        return;
}
/*
output(u, mx, mm)
float u[NX1];
int mx, mm;
{
	char	z[24][61];
	int	i, j;

	if( mm == 1 ){
		for( j = 1; j <= 23; j++ ){
			for( i = 2; i <= (mx - 1); i++ ){
				z[i][j] = ' ';
			}
		}
		for( j = 1; j <= 23; j++ ){
			z[1][j] = ':';
			z[mx][j] = ':';
		}
		for( i = 1; i <= mx; i++ ){
			z[i][1] = '-';
			z[i][23] = '-';
		}
	}

	if( mm == 2 ){
		for( i = 1; i <= mx; i++ ){
			j = u[i]*40 + .001;
			if( j > 23 )
				goto L_40;
			z[i][j] = '*';
L_40:
			;
		}
	}

	if( mm == 3 ){
		for( j = 23; j >= 1; j-- ){
			printf( " " );
			for( i = 1; i <= 60; i++ ){
				printf( "%c", z[i][j] );
			}
			printf( "\n" );
		}
	}
	return;
}
*/
