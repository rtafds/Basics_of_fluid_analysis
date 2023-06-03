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

void output_anime(const float *uw, int mx, int step);

int main(void)
{
    int i, ih, k, km, kop, counter;
    float dt, dx, r, x, alpha;
    /***** SETTING PARAMETERS */
    // ここでパラメータを設定します。

    mx=20; km=250;  // mx:メッシュ数. km:ステップ数
    alpha = 1.;  // 熱伝達率。そのそも時間も長さも単位を定義していないので、単位はあまり考えない。
    dt=0.001;  // Δt
    dx = 1./(float)( mx - 1 );  // メッシュ幅
    r = alpha * dt/(dx*dx);
    ih = (mx + 1)/2;  // 初期条件に使うパラメータ
    kop = 5;  // 書き出すステップ間隔
    counter = 1;  // ステップ間隔で書き出すのたるいので使う変数

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
        u[mx] = 2*r*u[mx - 1] + (1. - 2*r)*u[mx];
        
        if( (k%kop) == 1 ){
            output_anime(u, mx, counter);  // animeの出力
            counter++;
        }
        for( i = 2; i <= (mx - 1); i++ ){
            uu[i] = r*u[i - 1] + (1. - 2*r)*u[i] + r*u[i + 1];
        }
        for( i = 1; i <= mx; i++ ){
            u[i] = uu[i];
        }
        if( fabs( u[ih] ) >= 10000. ){
            printf( "DIVERGE!\n" );

    		exit(1);
        }
    }
}

void output_anime(const float *uw, int mx, int step)
{
    int i;
    char out_dir[256] = "./anime";
    char out_path[256];

    sprintf(out_path,"./anime/temp_%05d.dat",step);
    
    mkdir(out_dir, 0777); // sys.stat.hがないと機能しない。その場合は手動でanimeのフォルダを作ってコメントアウトする。
    
    FILE *file;
    file = fopen(out_path, "w");
    for(i=1; i<=mx; i++){
        fprintf(file, "\t%d\t%f\n",i,uw[i]);
    }
    fclose(file);
}
