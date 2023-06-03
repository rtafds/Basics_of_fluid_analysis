/*!
 * ガウス-ザイデル反復法(Gauss Seidel iterative method)
 *  - 解が収束するのは
 *      ・対角有利(diagonal dominant, 対角要素の絶対値>その行の他の要素の絶対値の和)
 *      ・係数行列が対称(symmetric)かつ正定(positive definite)
 *      ・Σ_j |a_ij/a_ii| < 1 (i = 1～n, j != i) 
 * @param[inout] A n×nの係数行列とn×1の定数項(b)を併せたn×(n+1)の行列．n+1列目に解が入る．
 * @param[in] n n元連立一次方程式
 * @param[inout] max_iter 最大反復数(反復終了後,実際の反復数を返す)
 * @param[inout] eps 許容誤差(反復終了後,実際の誤差を返す) 
 * @return 1:成功,0:失敗
 */
int GaussSeidel(vector< vector<double> > &A, int n, int &max_iter, double &eps)
{
    vector<double> x(n, 0.0);    // 初期値はすべて0とする
    double tmp;
 
    double e = 0.0;
    int k;
    for(k = 0; k < max_iter; ++k){
        // 現在の値を代入して，次の解候補を計算
        e = 0.0;
        for(int i = 0; i < n; ++i){
            tmp = x[i];
            x[i] = A[i][n];
            for(int j = 0; j < n; ++j){
                x[i] -= (j != i ? A[i][j]*x[j] : 0.0);
            }
            x[i] /= A[i][i];
 
            e += fabs(tmp-x[i]);    // 絶対誤差の場合
            //e += fabs((tmp-x[i])/tmp);    // 相対誤差の場合
        }
 
        // 収束判定
        if(e <= eps){
            break;
        }
    }
 
    max_iter = k;
    eps = e;
 
    for(int i = 0; i < n; ++i){
        A[i][n] = x[i];
    }
 
    return 1;
}