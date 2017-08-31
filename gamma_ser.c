#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double pi = 3.141592653589793238;
double euler = 0.577215664901532860;
double zeta[] = {
1.644934066848226436, //2
1.202056903159594285, //3
1.082323233711138191, //4
1.036927755143369926, //5
1.017343061984449139, //6
1.008349277381922826, //7
1.004077356197944339, //8
1.002008392826082214, //9
1.000994575127818085, //10
1.000494188604119464, //11
1.000246086553308048, //12
1.000122713347578489, //13
1.000061248135058704, //14
1.000030588236307020, //15
1.000015282259408651, //16
1.000007637197637899, //17
1.000003817293264999, //18
1.000001908212716553, //19
1.000000953962033872, //20
1.000000476932986787, //21
1.000000238450502727, //22
1.000000119219925965, //23
1.000000059608189051, //24
1.000000029803503514, //25
1.000000014901554828, //26
1.000000007450711789, //27
1.000000003725334024, //28
1.000000001862659723, //29
1.000000000931327432, //30
};

double lngamma(double x){
  double t=x-1;
  double sum=-log(x)+t*(1-euler);
  double pow=t*t;
  int n;
  for(n=2; n<=30; n++){
    sum+=(zeta[n-2]-1)*pow/n;
    pow*=(-t);
  }  
  return sum;
}

double gamma(double x){
  if((x<=0.0)&&(round(x)==x)){
    printf("Error: argument cannot be zero or negative integer\n");
    exit(0);
  }
  if(x>1.5) return (x-1)*gamma(x-1);
  if(x<0.5) return gamma(x+1)/x;
  return exp(lngamma(x));
}

int main(int argc,char *argv[]){
  if(argc<2){
    double g1=gamma(0.5);
    double g2=gamma(1.0);
    double g3=gamma(1.5);
    printf("gamma(0.5) = %.17lf mistake = %le\n",g1,g1-sqrt(pi));
    printf("gamma(1.0) = %.17lf mistake = %le\n",g2,g2-1);
    printf("gamma(1.5) = %.17lf mistake = %le\n",g3,g3-sqrt(pi)/2);
  }else{
    double x=atof(argv[1]);
    printf("gamma(%lf) = %.17lf\n",x,gamma(x));
  }
  return 0;
}
