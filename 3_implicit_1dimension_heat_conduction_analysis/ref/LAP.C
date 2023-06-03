#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/***********************************************************************
 *     LAPLACE EQUATION                                                *
 ***********************************************************************/
#define	NX	51
#define	NY	51
#define	NX1	NX+1
#define	NY1	NY+1
float u[NX][NY], uu[NX][NY];
void out(int, int);

int main(void)
{
int i, j, mx, my, n, nn;
float bb, eps, gosa;

	/***** INPUT & CALCULATE PARAMETERS */
	printf("MX,MY (<52 :MESH POINTS) ? (21,21) \n" );
	scanf("%d,%d", &mx, &my );
	printf("NN (NUMBER OF ITERATION) ? (1000)  \n" );
	scanf("%d", &nn );
	printf("EPSIRON (0.00001)? \n" );
	scanf("%f", &eps );

	/***** INITIAL CONDITION */
	for( j = 1; j <= my; j++ ){
		for( i = 1; i <= mx; i++ ){
			u[i][j] = 0.;
			}
		}

	/***** BOUNDARY CONDITION */
	for( i = 1; i <= mx; i++ ){
		u[i][1] = 1.0;
		u[i][mx] = 0.0;
		}
	for( j = 1; j <= my; j++ ){
		u[1][j] = 0.5;
		u[mx][j] = 0.0;
		}

	/***** MAIN LOOP */
	for( n = 1; n <= nn; n++ ){
 	 for( j = 2; j <= (my - 1); j++ ){
	  for( i = 2; i <= (mx - 1); i++ ){
		uu[i][j] = u[i][j];
		u[i][j] = .25*(u[i+1][j]+u[i][j+1]+u[i-1][j]+u[i][j-1]);
	  }
	 }

		gosa = 0.;
		for( j = 2; j <= (my - 1); j++ ){
			for( i = 2; i <= (mx - 1); i++ ){
				bb = fabs( (uu[i][j] - u[i][j])/u[i][j] );
				if( bb > gosa )
					gosa = bb;
				}
			}
		if( gosa < eps )
			goto L_70;
		if( (n%50) == 0 )
			{
			printf("N=%d GOSA=%g \n", n, gosa );
			}
		}

	/***** OUTPUT */
L_70:
	out(mx,my );

	exit(0);
}

void out(int mx, int my)
{
int i, index[NX1][NY1], j;
float umax, umin;

	umin = u[i][j];
	umax = u[i][j];
	for( j = 1; j <= my; j++ ){
		for( i = 1; i <= mx; i++ ){
			if( umax < u[i][j] )
				umax = u[i][j];
			if( umin > u[i][j] )
				umin = u[i][j];
			}
		}
	for( j = 1; j <= my; j++ ){
		for( i = 1; i <= mx; i++ ){
			index[i][j]=(int)((u[i][j]-umin)/(umax-umin)*9.9999)*11;
			}
		}
	for( j = my; j >= 1; j-- ){
		for( i = 1; i <= 21; i++ ){
			printf("%2d", index[i][j] );
			}
		printf("\n" );
		}

	return;
}