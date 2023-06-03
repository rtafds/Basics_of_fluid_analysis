#include <stdio.h>

#define MX  21
#define MY  21
#define MX1  MX+1
#define MY1  MY+1
float u[MX1][MY1],uu[MX1][MY1];
void out(int, int);
int main(void)
{
    float dt,dx,dy,r1,r2;
    int   nlast,j,k,n;
    //printf(" DT?,NLAST? (0.0005,400) \n");
    //scanf("%f,%d",&dt,&nlast);
    dt=0.0005;
    nlast=400;
    dx=1./(float)(MX-1);
    dy=1./(float)(MY-1);
    r1=dt/(dx*dx);
    r2=dt/(dy*dy);
    for(k=1;k<=MY;k++){
        for(j=1;j<=MX;j++){
            u[j][k]=0.;
            uu[j][k]=0.;
        }
    }
    // main loop
    for(n=1;n<=nlast;n++){
        for(k=1;k<=MY;k++){
            u[1][k]=0.5; u[MX][k]=0.;
        }
        for(j=1;j<=MX;j++){
            u[j][1]=1.0; u[j][MY]=0.;
        }
        for(k=2;k<=MY-1;k++){
            for(j=2;j<=MX-1;j++){
                uu[j][k]=u[j][k]+r1*(u[j+1][k]-2.*u[j][k]+u[j-1][k])
                        +r2*(u[j][k+1]-2.*u[j][k]+u[j][k-1]);
            }
        }
        for(k=2;k<=MY-1;k++){
            for(j=2;j<=MX-1;j++){
                u[j][k]=uu[j][k];
            }
        }
        if(n-(int)(n/100)*100==0){
            out(MX,MY);
        }
    }
    return 0;
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