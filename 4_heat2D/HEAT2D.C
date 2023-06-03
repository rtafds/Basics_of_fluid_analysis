#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>

#define MX  21  // mesh number of x
#define MY  21  // mesh number of y
#define MX1  MX+1
#define MY1  MY+1
float u[MX1][MY1],uu[MX1][MY1];

void output_anime(const float uw[MX1][MY1], int mx, int my, int step);

int main(void)
{
    float dt,dx,dy,r1,r2,alpha;
    int   nlast,i,j,n,nop;

    dt=0.0005;
    nlast=400;  // ループ回数
    alpha=1.;  // 熱伝達汁


    dx=1./(float)(MX-1);
    dy=1./(float)(MY-1);
    r1=dt/(dx*dx);
    r2=dt/(dy*dy);
    nop = 5;  // 書き出すステップ間隔

    for(i=1;i<=MX;i++){
        for(j=1;j<=MY;j++){
            u[i][j]=0.;
            uu[i][j]=0.;
        }
    }

    // main loop
    for(n=0;n<=nlast;n++){
        
        // 初期条件を決める。今回は0で行くのでなし。
        // 書き込む場合はメッシュに相当する配列の部分に値を与えればいい。

        if( n%nop == 0){
            output_anime(u, MX, MY, n);
        }

        // 境界条件の設定
        for(i=1;i<=MX;i++){
            u[i][1]=0.5; u[i][MX]=0.;
        }
        for(j=1;j<=MY;j++){
            u[1][j]=1.0; u[MY][j]=0.;
        }
        for(i=2;i<=MX-1;i++){
            for(j=2;j<=MY-1;j++){
                uu[i][j]=u[i][j]+r1*(u[i][j+1]-2.*u[i][j]+u[i][j-1])
                        +r2*(u[i+1][j]-2.*u[i][j]+u[i-1][j]);
            }
        }
        for(i=2;i<=MX-1;i++){
            for(j=2;j<=MY-1;j++){
                u[i][j]=uu[i][j];
            }
        }
    }
    return 0;
}


void output_anime(const float uw[MX1][MY1], int mx, int my, int step)
{
    int i,j;
    char out_dir[256] = "./anime";
    char out_path[256];

    sprintf(out_path,"./anime/temp_%05d.dat",step);
    
    mkdir(out_dir, 0777); // sys.stat.hがないと機能しない。その場合は手動でanimeのフォルダを作ってコメントアウトする。
    
    FILE *file;
    file = fopen(out_path, "w");
    for(i=1; i<=mx; i++){
        for(j=1; j<=my; j++){
            fprintf(file, "\t%d\t%d\t%f\n",i,j,uw[i][j]);
        }
    }
    fclose(file);
}



void out(int i21, int j21)
{
    int i, ind[100], j;
    float qma, qmi;
    printf("  \n" );
    qma = u[1][1];
    qmi = u[1][1];
    for(j=1;j<=j21;j++){
        for(i=1;i<=i21;i++){
            if(qma<u[i][j]){
                qma = u[i][j];
            }
            if(qmi>u[i][j]){
                qmi = u[i][j];
            }
        }
    }

    for(j=j21;j>=1;j+=-1){
        for(i=1;i<=i21;i++){
            ind[i]=(int)((u[i][j]-qmi)/(qma-qmi)*9.9999);
        }
        for(i=1;i<=i21;i++){
            printf("%d ",ind[i]);
        }
        printf("\n");
    }
    return;
}