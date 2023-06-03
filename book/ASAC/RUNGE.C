
/* FOURTH ORDER RUNGE-KUTTA METHOD */
#include <stdio.h>
float f(float, float);
int main(void)
{
    float h,xlast,x,y;
    float yexact,err;
    float del1,del2,del3,del4;
    printf(" INPUT H & XLAST ?\n");
    //scanf("%f,%f",&h,&xlast);
    h=0.01;
    xlast=10.;
    x=0.; y=1.;
    do {
        yexact=4./(4.-2.*x-x*x);
        err=(yexact-y)/yexact*100.;
        printf(" %7.2f  %14.7f  %14.7f  %14.7f\n",x,y,yexact,err);
        del1=h*f(x,y); del2=h*f(x+h/2.,y+del1/2.);
        del3=h*f(x+h/2.,y+del2/2.); del4=h*f(x+h,y+del3);
        y=y+(del1+2.*del2+2.*del3+del4)/6.;
        x=x+h;
    } while ( x<=xlast);
    return 0;

}

float f(float a,float b)
{
    float c;
    c=.5*(1.+a)*b*b;
    return(c);
}