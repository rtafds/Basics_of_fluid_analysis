#include <stdio.h>
int main(void)
{
  float h,xlast,x,y;
  float yexact,err;
  float del1,del2,del3,del4;
  x=1.; y=2.;
  yexact=4./(4.-2.*x-x*x);
  printf("%f\n",yexact);
  return 0;
}