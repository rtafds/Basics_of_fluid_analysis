#include <ctype.h>
#include <math.h>
#include <stdio.h>

float a[100],b[100],c[100],d[100];

void thomas(int il,int iu);

int main(void)
{
    float h,x,u,err;
    int   n,i;
    FILE *file;
    file = fopen("output.dat", "w");
    
    //printf("How many mesh?\n");
    //scanf("%d",&n);
    n = 20;  // mesh number

    // 係数の初期値を入れる
    h=1./(float)n;
    for(i=2;i<=n;i++){
        a[i]=1.0; b[i]=h*h-2.0;
        c[i]=1.0; d[i]=-(float)(i-1)*h*h*h;
    }

    // thomasで、dに答えが入る感じになる。
    thomas(2,n);

    // 解析解を作って、差分近似で計算した解と照らし合わせる
    for(i=2;i<=n;i++){
        x=h*(float)(i-1);
        u=sin(x)/sin(1.0)-x;
        err=(d[i]-u)/u*100.0;
        fprintf(file, "\t%e\t%e\t%e\t%e\n",x,d[i],u,err);
    }
    fclose(file);

    return 0;
}
/* */
void thomas(int il,int iu)
{
    int   ip,i,j;
    float r;
    ip=il+1;  // 最初はg=b, s=dなのでそのまま使用。
    // 係数を算出していく
    for(i=ip;i<=iu;i++){
    r=c[i]/b[i-1]; 
        b[i]=b[i]-r*a[i-1];  // gを算出。最初がg=bなのでbをそのままg
        d[i]=d[i]-r*d[i-1];  // sを算出
    }

    // 最後の係数を算出
    d[iu]=d[iu]/b[iu];
    // 解を求める。
    for(i=ip;i<=iu;i++){
        j=iu-i+il;
        d[j]=(d[j]-a[j]*d[j+1])/b[j];
    }
    return;
}