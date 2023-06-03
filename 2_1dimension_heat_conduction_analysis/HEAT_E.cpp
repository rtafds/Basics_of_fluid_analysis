#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sys/stat.h>

using namespace std;

/***********************************************************************
 *     DIFFUSION EQUATION   EULER EXPLICIT METHOD                      *
 *********************************************************************** */

void output_anime(const vector<double>& uw, int mx, int step);

int main(void)
{
    int i, ih, k, mx,km, kop, counter;
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

    vector<double> u(mx+1,0.0),uu(mx+1,0.0);  // 0から始まり、u[mx]はmx+1個目の値

    /***** INITIAL CONDITION */
    for( i = 1; i <= mx; i++ ){
        x = (double)( i - 1 )/(double)( mx - 1 );
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
        
        if( (k%kop) == 1 ){
            output_anime(u, mx, counter);  // animeの出力
            counter++;
        }
        for( i = 2; i <= (mx - 1); i++ ){
            uu[i] = r*u[i - 1] + (1. - 2*r)*u[i] + r*u[i + 1];
        }
        u = uu;  // 元の配列にコピー

        if( fabs( u[ih] ) >= 10000. ){
            cout<< "DIVERGE!"<<endl;

    		exit(1);
        }
    }
}

// この辺はC仕様の方が便利なので、C仕様
void output_anime(const vector<double>& uw, int mx, int step)
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

