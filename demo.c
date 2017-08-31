#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "specfun.h"

int main(int argc,char *argv[]){
  if(argc<2){
    printf("Usage:\ngamma x(double)\n");
    printf("beta x(double) y(double)\n");
    printf("besselJ nu(double) x(double)\n");
    printf("besselY nu(double) x(double)\n");
    printf("besselI nu(double) x(double)\n");
    printf("besselK nu(double) x(double)\n");
    return 0;
  }
  if(strcmp(argv[1],"gamma")==0){
    double x=atof(argv[2]);
    printf("gamma(%s) = %.17lf\n",argv[2],eulergamma(x));
  }
  if(strcmp(argv[1],"besselJ")==0){
    double nu=atof(argv[2]);
    double x=atof(argv[3]);
    printf("besselJ(%s,%s) = %.17lf\n",argv[2],argv[3],besselJ(nu,x));
  }
  return 0;
}
