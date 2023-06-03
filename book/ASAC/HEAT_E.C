
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>

/***********************************************************************
 *     DIFFUSION EQUATION   EULER EXPLICIT METHOD                      *
 *********************************************************************** */
#define	NX	51
#define	NX1	NX+1

float u[NX1], uu[NX1];
int z[61][24];
int mx;

void output(int);
void output_anime(float *uw, int mx);

int main(void)
{
    int i, i5, ih, k, km;
    float dt, dx, r, x;
    /***** INPUT & CALCULATE PARAMETERS */
    //printf("MESH POINTS (<52) ? NUMBER OF TIME STEP ?? (20,250)\n" );
    //scanf("%d,%d", &mx, &km );
    //printf("DELTA T: DT ? (0.001) \n" );
    //scanf("%f", &dt );
    mx=20; km=250;
    dt=0.001;
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

    /***** MAIN LOOP */
    for( k = 1; k <= km; k++ ){
        u[1] = 0.;
        u[mx] = 0.;
        if( (k%i5) == 1 ){
                output( 2 );
        }
        //output_anime(*u, mx)
        for( i = 2; i <= (mx - 1); i++ ){
                uu[i] = r*u[i - 1] + (1. - 2*r)*u[i] + r*u[i + 1];
                }
        for( i = 2; i <= (mx - 1); i++ ){
                u[i] = uu[i];
                }

        if( fabs( u[ih] ) >= 10000. ){
                printf( "DIVERGE!\n" );

    			exit(1);
        }
    }
    output( 3 );
    output_anime(uu, mx);
}

void output_anime(float *uw, int mx)
{
    int i;
    char out_dir[256] = "./anime";
    char out_path[256] = "./anime/test.dat";
    
    mkdir(out_dir, 0777); // 作成成功の場合は0が返る
    
    
    FILE *file;
    file = fopen(out_path, "w");
    for(i=1; i<mx; i++){
        fprintf(file, "\t%d\t%f\n",i,uw[i]);
    }
    fclose(file);

}


void output(int mm)
{
    int i, j;
    if( mm == 1 ){
        for(j = 1;j <= 23;j++ ){
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