#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
/***********************************************************************
 *     DIFFUSION EQUATION   EULER IMPLICIT METHOD                      *
 *********************************************************************** */
#define	NX	51
#define	NX1	NX+1

float u[NX1], uu[NX1],a[NX1],b[NX1],c[NX1],d[NX1],s[NX1],p[NX1];
int z[61][24];
int mx;

void output_anime(const float *uw, int mx, int step);

int main(void)
{
int i, ih, k, km, ii, kop, counter;
float dt, dx, r, x, alpha;
	/***** INPUT & CALCULATE PARAMETERS */

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

	/***** IMPLICIT METHOD */
	for( k = 1; k <= km; k++ ){
		u[1] = 0.;
		u[mx] = 0.;

        if( (k%kop) == 1 ){
            output_anime(u, mx, counter);  // animeの出力
            counter++;
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
