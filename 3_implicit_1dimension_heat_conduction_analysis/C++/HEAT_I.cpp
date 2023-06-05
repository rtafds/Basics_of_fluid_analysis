#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sys/stat.h>

using namespace std;

/***********************************************************************
 *     DIFFUSION EQUATION   EULER IMPLICIT METHOD                      *
 *********************************************************************** */

// よりc++っぽいプログラムにします。

void output_anime(const vector<double>& uw, int step);

int JacobiIteration(vector< vector<double> > &A, int n, int &max_iter, double &eps);

int main(void)
{
	int i, ih, k, mx, km, ii, kop, counter, end_of_cal, max_iter;
	double dt, dx, r, x, alpha, eps;
	/***** INPUT & CALCULATE PARAMETERS */

    mx=20; km=250;  // mx:メッシュ数. km:ステップ数
    alpha = 1.;  // 熱伝達率。そのそも時間も長さも単位を定義していないので、単位はあまり考えない。
    dt=0.001;  // Δt
    dx = 1./(double)( mx - 1 );  // メッシュ幅
    r = alpha * dt/(dx*dx);
    ih = (mx + 1)/2;  // 初期条件に使うパラメータ
    kop = 5;  // 書き出すステップ間隔
    counter = 1;  // ステップ間隔で書き出すのたるいので使う変数


	// 0から始める仕様で書きます。
	vector<double> u(mx,0.0), uu(mx,0.0);  // 値と次の値を入れる。
    
    // A n×nの係数行列とn×1の定数項(b)を併せたn×(n+1)の拡大行列．
	vector< vector<double> > A(mx, vector<double>(mx+1, 0.0));

	/***** INITIAL CONDITION */
	for( i = 0; i < mx; i++ ){
		x = (double)( i )/(double)( mx - 1 );
		if( i < ih ){
			u[i] = x;
		}
		else{
			u[i] = 1. - x;
		}
	}

	/***** MAIN LOOP *****/
	for( k = 1; k <= km; k++ ){
		u[0] = 0.;
		u[mx-1] = 0.;

        if( (k%kop) == 1 ){
            output_anime(u, counter);
            counter++;
        }
		
        // 係数の設定
        // 本来は境界まで係数行列に入れなくてもいいが、色々めんどいため入れた。
        // そのため、境界はさっき代入したのが、k+1の値でそのままになって欲しいので、対角成分にだけ1を入れる。
        A[0][0] = 1;
        A[mx-1][mx-1] = 1;
        A[0][mx] = u[0];  // 右辺
        A[mx-1][mx] = u[mx-1];
		for(i=1; i<mx-1; i++){
			A[i][i-1]= -r;
			A[i][i]= 1.+ 2.*r;
			A[i][i+1]= -r;
			A[i][mx]= u[i];  // 右辺
		}

		max_iter = 1000;
		eps = 0.0001;
        
        // ここで陰解法でとく。
		JacobiIteration(A, mx, max_iter, eps);

        // 解を入れる。Aの一番右が解になっている。
		for(i=0; i<mx; i++){
			u[i] = A[i][mx];
		}
	}
}

/*!
 * ヤコビ反復法(Jacobi iterative method)
 *  - 解が収束するのは
 *      ・対角有利(diagonal dominant, 対角要素の絶対値>その行の他の要素の絶対値の和)
 *      ・係数行列が対称(symmetric)かつ正定(positive definite)
 *    のどちらかの場合
 * @param[inout] A n×nの係数行列とn×1の定数項(b)を併せたn×(n+1)の拡大行列．n+1列目に解が入る．
 * @param[in] n n元連立一次方程式
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す) 
 * @return 1:成功,0:失敗
 */

int JacobiIteration(vector< vector<double> > &A, int n, int &max_iter, double &eps)
{
    vector<double> x(n, 0.0);    // 初期値はすべて0とする
    vector<double> y(n, 0.0);
 
    double e = 0.0;
    int k;
    for(k = 0; k < max_iter; ++k){
        // 現在の値を代入して，次の解候補を計算
        for(int i = 0; i < n; ++i){
            y[i] = A[i][n];
            for(int j = 0; j < n; ++j){
                //y[i] -= (j != i ? A[i][j]*x[j] : 0.0); こうやっても書ける
                if (i != j){
                    y[i] -= A[i][j]*x[j];
                }
            }
            y[i] /= A[i][i];
        }

        // 収束判定
        e = 0.0;
        for(int i = 0; i < n; ++i){
            e += fabs(y[i]-x[i]);    // 絶対誤差の場合
            //e += fabs((y[i]-x[i])/y[i]);    // 相対誤差の場合
        }
        if(e <= eps){
            break;
        }
        swap(x, y);
    }
    max_iter = k;
    eps = e;
    for(int i = 0; i < n; ++i){
        A[i][n] = y[i];
    }
    return 1;
}


void output_anime(const vector<double>& uw ,int step)
{
    int i;
    char out_dir[256] = "./anime";
    char out_path[256];

    sprintf(out_path,"./anime/temp_%05d.dat",step);
    
    mkdir(out_dir, 0777); // sys.stat.hがないと機能しない。その場合は手動でanimeのフォルダを作ってコメントアウトする。
    
    FILE *file;
    file = fopen(out_path, "w");
    for(i=0; i<(int)uw.size(); i++){
        fprintf(file, "\t%d\t%f\n",i,uw[i]);
    }
    fclose(file);
}
