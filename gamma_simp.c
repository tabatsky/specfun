#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Pi 3.141592653589793238

int N;

double gamma(double x){
  int i;
  double yp,y,yn;
  double h=40.0/N;
  double gamma=0;
  if(x==1) yp=1;
    else yp=0;
  y=exp(-h+log(h)*(x-1));
  yn=exp(-2*h+log(2*h)*(x-1));
  gamma+=(yp+4*y+yn)*h/3;
  for(i=3; i<N; i+=2){
    yp=yn;
    y=exp(-i*h+log(i*h)*(x-1));
    yn=exp(-(i+1)*h+log((i+1)*h)*(x-1));
    gamma+=(yp+4*y+yn)*h/3;
  }
  return gamma;
}

int main(int argc, char *argv[]){
  double g1,g2,g3;
  if(argc<2) N=40000;
    else N=atoi(argv[1]);
  g1=gamma(1.0);
  g2=gamma(1.5);
  g3=gamma(2.0);
  printf("gamma(1.0) = %.17lf\n",g1);
  printf("mistake = %le\n",g1-1);
  printf("gamma(1.5) = %.17lf\n",g2);
  printf("mistake = %le\n",g2-sqrt(Pi)/2);
  printf("gamma(2.0) = %.17lf\n",g3);
  printf("mistake = %le\n",g3-1);
  return 0;
}