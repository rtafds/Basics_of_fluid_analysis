#include <stdio.h>
#include <math.h>
/***********************************************************************
 *     DIFFUSION EQUATION   EULER IMPLICIT METHOD                      *
 *********************************************************************** */
#define	NX	51
#define	NX1	NX+1

float u[NX1], uu[NX1],a[NX1],b[NX1],c[NX1],d[NX1],s[NX1],p[NX1];
int z[61][24];
int mx;

main()
{
int i, i5, ih, k, km, ii;
float dt, dx, r, x;
	/***** INPUT & CALCULATE PARAMETERS */
	printf("MESH POINTS(<52) ?  NUMBER OF TIME STEP ?? (40,150) \n" );
	scanf("%d,%d", &mx, &km );
	printf("DELTA T: DT ? (0.002) \n" );
	scanf("%f", &dt );
	dx = 1./(float)( mx - 1 );
	r = dt/(dx*dx);
	ih = (mx + 1)/2;
	i5 = .05/dt;
	if( i5 == 0 ){
	  i5 = 1;
        }
	output( 1 );

	/***** INITIAL CONDITION */
	for( i = 1; i <= mx; i++ ){
		x = (float)( i - 1 )/(float)( mx - 1 );
		if( i <= ih ){
			u[i] = x;
			}
		else{
			u[i] = 1. - x;
			}
		}

	/***** IMPLICIT METHOD */
	for( k = 1; k <= km; k++ ){
		u[1] = 0.;
		u[mx] = 0.;
		if( (k%i5) == 1 ){
			output( 2 );
		}
		for( i=1; i<=mx; i++){
		 a[i]=r;
		 b[i]=-2.*r-1.;
		 c[i]=r;
		 d[i]=-u[i];
		}
		p[1]=b[1];
		s[1]=d[1];
		for( i=2; i<=mx; i++){
		 p[i]=b[i]-a[i]*c[i-1]/p[i-1];
		 s[i]=d[i]-a[i]*s[i-1]/p[i-1];
		}
		uu[mx]=s[mx]/p[mx];
		for( i = 2; i <= (mx - 1); i++ ){
			ii=mx+1-i;
			uu[ii] = (s[ii]-c[ii]*uu[ii+1])/p[ii];
		}
		for( i = 2; i <= (mx - 1); i++ ){
			u[i] = uu[i];
		}

	}
	output( 3 );
}

output(mm)
int mm;
{
int i, j;


        if( mm == 1 ){
        for( j = 1; j <= 23; j++ ){
                for( i = 1; i <= 60; i++ ){
                 z[i][j]=0;
                }
        }
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
