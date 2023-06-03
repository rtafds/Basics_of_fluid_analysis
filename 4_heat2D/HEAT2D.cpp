#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <sys/stat.h>

using namespace std;

void output_anime(const vector< vector<double> > &uw, int mx, int my, int step);

int main(void)
{
    int   nlast,i,j,n,nop, mx, my;
    double dt,dx,dy,r1,r2,alpha;

    mx = 21; my=21;  // mesh number of x and y

    dt=0.0005;
    nlast=400;  // ループ回数
    alpha=1.;  // 熱伝達汁

    dx=1./(double)(mx-1);
    dy=1./(double)(my-1);
    r1=dt/(dx*dx);
    r2=dt/(dy*dy);
    nop = 5;  // 書き出すステップ間隔

    // 配列を定義します。
    vector< vector<double> > u(mx, vector<double>(my, 0.0));
    vector< vector<double> > uu(mx, vector<double>(my, 0.0));

    for(i=0;i<mx;i++){
        for(j=0;j<my;j++){
            u[i][j]=0.;
            uu[i][j]=0.;
        }
    }

    // main loop
    for(n=0;n<=nlast;n++){

        // 初期条件を決める。今回は0で行くのでなし。
        // 書き込む場合はメッシュに相当する配列の部分に値を与えればいい。

        if( n%nop == 0){
            //output_anime(u, mx, my, n);
        }

        // 境界条件の設定
        for(i=0;i<mx;i++){
            u[i][0]=0.5; u[i][mx-1]=0.;
        }
        for(j=0;j<my;j++){
            u[0][j]=1.0; u[my-1][j]=0.;
        }
        for(i=1;i<mx-1;i++){
            for(j=1;j<my-1;j++){
                uu[i][j]=u[i][j]+r1*(u[i][j+1]-2.*u[i][j]+u[i][j-1])
                        +r2*(u[i+1][j]-2.*u[i][j]+u[i-1][j]);
            }
        }
        for(i=1;i<mx-1;i++){
            for(j=1;j<my-1;j++){
                u[i][j]=uu[i][j];
            }
        }
    }
    return 0;
}


void output_anime(const vector< vector<double> > &uw, int mx, int my, int step)
{
    int i,j;
    char out_dir[256] = "./anime";
    char out_path[256];

    sprintf(out_path,"./anime/temp_%05d.dat",step);
    
    mkdir(out_dir, 0777); // sys.stat.hがないと機能しない。その場合は手動でanimeのフォルダを作ってコメントアウトする。
    
    FILE *file;
    file = fopen(out_path, "w");
    for(i=0; i<mx; i++){
        for(j=0; j<my; j++){
            fprintf(file, "\t%d\t%d\t%f\n",i,j,uw[i][j]);
        }
    }
    fclose(file);
}
