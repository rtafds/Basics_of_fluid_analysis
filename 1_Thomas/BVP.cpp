#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

// ベクトルの参照渡し
void thomas(int il,int iu,vector<double>& a,vector<double>& b,vector<double>& c,vector<double>& d);

int main(void){
    
    float h,x,u,err;
    int   n,i;
    
    string filename = "output.dat";
    ofstream writing_file;
    writing_file.open(filename, ios::out);

    n = 20;  // mesh number

    // 要素数nで、0.0に初期化した配列を作成するが、なぜか2から始めるため一個多く確保する。
    vector<double> a(n+1,0.0),b(n+1,0.0),c(n+1,0.0),d(n+1,0.0);
    
    h=(double)1./n;
    for(i=2;i<=n;i++){
        a[i]=1.0; 
        b[i]=h*h-2.0;
        c[i]=1.0; 
        d[i]=(double)-(i-1)*h*h*h;
    }
    
    thomas(2,n, a,b,c,d);
    
    for(i=2;i<=n;i++){
        x=(double)h*(i-1);
        u=sin(x)/sin(1.0)-x;
        err=(d[i]-u)/u*100.0;
        // 出力
        writing_file
        <<"\t"<<scientific<<x
        <<"\t"<<scientific<<d[i]
        <<"\t"<<scientific<<u
        <<"\t"<<scientific<<err<<endl;
    }
    writing_file.close();
    return 0;
}

void thomas(int il,int iu, vector<double>& a,vector<double>& b,vector<double>& c,vector<double>& d)
{
    int   ip,i,j;
    float r;
    ip=il+1;  // 最初はg=b, s=dなのでそのまま使用。
    for(i=ip;i<=iu;i++){
        r=c[i]/b[i-1]; 
        b[i]=b[i]-r*a[i-1];  // gを算出。最初がg=bなのでbをそのままg
        d[i]=d[i]-r*d[i-1];  // sを算出
    }
    d[iu]=d[iu]/b[iu];
    for(i=ip;i<=iu;i++){
        j=iu-i+il;
        d[j]=(d[j]-a[j]*d[j+1])/b[j];
    }
    return;
}