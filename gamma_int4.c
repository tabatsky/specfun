#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Pi 3.141592653589793228

int N;

double gamma(double x){
  int i;
  double ypp,yp,y,yn,ynn;
  double h=40.0/N;
  double gamma=0;
  if(x==1) ypp=1;
    else ypp=0;
  yp=exp(-h+log(h)*(x-1));
  y=exp(-2*h+log(2*h)*(x-1));
  yn=exp(-3*h+log(3*h)*(x-1));
  ynn=exp(-4*h+log(4*h)*(x-1));
  gamma+=(14*ypp+64*yp+24*y+64*yn+14*ynn)*h/45;
  for(i=6; i<N; i+=4){
    ypp=ynn;
    yp=exp(-(i-1)*h+log((i-1)*h)*(x-1));
    y=exp(-i*h+log(i*h)*(x-1));
    yn=exp(-(i+1)*h+log((i+1)*h)*(x-1));
    ynn=exp(-(i+2)*h+log((i+2)*h)*(x-1));
    gamma+=(14*ypp+64*yp+24*y+64*yn+14*ynn)*h/45;
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