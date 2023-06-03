#include <ctype.h>
#include <math.h>
/***** 2D HEAT EQUATION -- ADI METHOD -- */
#define	MX	21
#define	MY	21
float a[100],b[100],c[100],d[100],u[MX][MY],uu[MX][MY];
main( )
{
 int j,k,n,nlast,imx,imy,icnst;
 float dt,dx,dy,r1,r2;
 imx=MX; imy=MY; icnst=2;

 printf("DT?,NLAST? \n");
 scanf("%f %d", &dt, &nlast );
 dx=1./(MX-1);
 dy=1./(MY-1);
 r1=.5*dt/(dx*dx);
 r2=.5*dt/(dy*dy);
 for(k=1;k<=MY;k++){
  for(j=1;j<=MX;j++){
   u[j][k]=0.;
   uu[j][k]=0.;
  }
 }

/* */
 for(n=1;n<=nlast;n++){
/* */
  for(k=1;k<=MY;k++){
   u[1][k]=0.;
   u[MX][k]=0.5;
   uu[1][k]=0.;
   uu[MX][k]=0.5;
  }
  for(j=1;j<=MX;j++){
   u[j][1]=0.;
   u[j][MY]=1.;
   uu[j][1]=0.;
   uu[j][MY]=1.;
  }
/* */
 for(k=2;k<=(MY-1);k++ ){
  for(j=2;j<=(MX-1);j++){
   a[j]=-r1;
   b[j]=2.*r1+1.;
   c[j]=-r1;
   d[j]=u[j][k]+r2*(u[j][k+1]-2.*u[j][k]+u[j][k-1]);
  }
   d[2]=d[2]+r1*u[1][k];
   d[MX-1]=d[MX-1]+r1*u[MX][k];
/* */
   thomas(icnst,imx-1);
/* */
   for(j=2;j<=(MX-1);j++){
    uu[j][k]=d[j];
   }
 }
/* */
 for(j=2;j<=(MX-1);j++){
  for(k=2;k<=(MY-1);k++){
   a[k]=-r2;
   b[k]=2.*r2+1.;
   c[k]=-r2;
   d[k]=uu[j][k]+r1*(uu[j+1][k]-2.*uu[j][k]+uu[j-1][k]);
  }
   d[2]=d[2]+r2*uu[j][1];
   d[MY-1]=d[MY-1]+r2*uu[j][MY];
/* */
   thomas(icnst,imy-1);
/* */
  for(k=2;k<=(MY-1);k++){
   u[j][k]=d[k];
  }
 }
  if((float)(n-(int)(n/50)*50)==0.0){
   out(imx,imy);
 }
}
}
/* */
thomas(il,iu)
int il,iu;
 {
  int   i,ip,j;
  float r;
  ip=il+1;
  for(i=ip;i<=iu;i++){
   r=c[i]/b[i-1];
   b[i]=b[i]-r*a[i-1];
   d[i]=d[i]-r*d[i-1];
  }

  d[iu]=d[iu]/b[iu];
  for(i=ip;i<=iu;i++ ){
   j=iu-i+il;
   d[j]=(d[j]-a[j]*d[j+1])/b[j];
  }
 }
 /* */
out(mx, my)
int mx, my;
{
int i, index[41][41], j;
float umax, umin;
	i = 1;
	j = 1;
	umin = u[i][j];
	umax = u[i][j];
	for( j = 1; j <= my; j++ ){
		for( i = 1; i <= mx; i++ ){
			if( umax < u[i][j] )
				umax = u[i][j];
			if( umin > u[i][j] )
				umin = u[i][j];
			}
		}
	for( j = 1; j <= my; j++ ){
		for( i = 1; i <= mx; i++ ){
			index[i][j]=(int)((u[i][j]-umin)/(umax-umin)*9.9999)*11;
			}
		}
	for( j = my; j >= 1; j-- ){
		for( i = 1; i <= 21; i++ ){
			printf("%2d", index[i][j] );
			}
		printf("\n" );
		}

	return;
}
