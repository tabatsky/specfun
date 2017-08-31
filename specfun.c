#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

double eulergamma(double x){
  if((x<=0.0)&&(round(x)==x)){
    printf("Error: argument of gamma cannot be zero or negative integer\n");
    exit(0);
  }
  if(x>1.5) return (x-1)*eulergamma(x-1);
  if(x<0.5) return eulergamma(x+1)/x;
  return exp(lngamma(x));
}

double eulerbeta(double x,double y){
  return eulergamma(x)*eulergamma(y)/eulergamma(x+y);
}

long long int fact(int n){
  int i;
  long long int fact=1;
  for(i=2; i<=n; i++) fact*=2;
  return fact;
}

double besselJ(double nu,double x){
  double next,sum,next1,next2,sum1,sum2,int1,int2,yp,y,yn,h;
  int m,N;
  if((nu<0)&&(round(nu)==nu))
    if(((int)round(nu))%2==1) return -besselJ(-nu,x);
      else return besselJ(-nu,x);
  if((nu==0.0)&&(x==0.0)) return 1.0;
  if(x==0.0) return 0.0;
 // if(nu>1) return 2*(nu-1)/x*besselJ(nu-1,x)-besselJ(nu-2,x);
 // if(nu<-1) return 2*(nu+1)/x*besselJ(nu+1,x)-besselJ(nu+2,x);
 // power series
  if(x<5){ 
    next=exp(log(x/2)*nu)/eulergamma(nu+1);
    sum=next;
    for(m=1; m<200; m++){
      next*=(-1)*(x*x/4)/(m*(m+nu));
      sum+=next;
    }
    return sum;
  }
  // asymptotic expansion
  if((x>400)&&(x>(40*(nu*nu-0.25)))){
    next1=1;
    sum1=next1;
    next2=(4*nu*nu-1)/(8*x);
    sum2=next2;
    for(m=1; m<=200; m++){
      next1*=-(4*nu*nu-(4*m-3)*(4*m-3))*(4*nu*nu-(4*m-1)*(4*m-1))/((2*m-1)*2*m*8*x*8*x);
      sum1+=next1;
      next2*=-(4*nu*nu-(4*m-1)*(4*m-1))*(4*nu*nu-(4*m+1)*(4*m+1))/(2*m*(2*m+1)*8*x*8*x);
      sum2+=next2;        
    }
    return sqrt(2/(pi*x))*(cos(x-pi*nu/2-pi/4)*sum1-sin(x-pi*nu/2-pi/4)*sum2);
  }
  //integral
  N=(int)400*round(x);
  h=pi/N;
  yp=1;
  y=cos(nu*h-x*sin(h));
  yn=cos(nu*2*h-x*sin(2*h));
  int1=h/3*(yp+4*y+yn);
  for(m=3; m<N; m+=2){
    yp=yn;
    y=cos(nu*m*h-x*sin(m*h));
    yn=cos(nu*(m+1)*h-x*sin((m+1)*h));
    int1+=h/3*(yp+4*y+yn);
  }
  if(nu==round(nu)) return int1/pi;
  h=0.000004;
  yp=1;
  y=exp(-x*sinh(h)-nu*h);
  yn=exp(-x*sinh(2*h)-nu*h);
  int2=h/3*(yp+4*y+yn);
  for(m=3; m<500000; m+=2){
    yp=yn;
    y=exp(-x*sinh(m*h)-nu*m*h);
    yn=exp(-x*sinh((m+1)*h)-nu*(m+1)*h);
    int2+=h/3*(yp+4*y+yn);
  }
  return int1/pi-sin(nu*pi)*int2/pi;
}

double besselY(double nu,double x){
  double next1,next2,sum1,sum2,int1,int2,yp,y,yn,h;
  int m,N,sign;
  if((nu<0)&&(round(nu)==nu))
    if(((int)round(nu))%2==1) return -besselY(-nu,x);
      else return besselY(-nu,x);
  if(round(nu)!=nu) return (besselJ(nu,x)*cos(nu*pi)-besselJ(-nu,x))/sin(nu*pi);
  // asymptotic expansion
  if((x>400)&&(x>(40*(nu*nu-0.25)))){
    next1=1;
    sum1=next1;
    next2=(4*nu*nu-1)/(8*x);
    sum2=next2;
    for(m=1; m<=200; m++){
      next1*=-(4*nu*nu-(4*m-3)*(4*m-3))*(4*nu*nu-(4*m-1)*(4*m-1))/((2*m-1)*2*m*8*x*8*x);
      sum1+=next1;
      next2*=-(4*nu*nu-(4*m-1)*(4*m-1))*(4*nu*nu-(4*m+1)*(4*m+1))/(2*m*(2*m+1)*8*x*8*x);
      sum2+=next2;        
    }
    return sqrt(2/(pi*x))*(sin(x-pi*nu/2-pi/4)*sum1+cos(x-pi*nu/2-pi/4)*sum2);
  }
  //integral
  N=(int)400*round(x);
  h=pi/N;
  yp=0;
  y=sin(-nu*h+x*sin(h));
  yn=sin(-nu*2*h+x*sin(2*h));
  int1=h/3*(yp+4*y+yn);
  for(m=3; m<N; m+=2){
    yp=yn;
    y=sin(-nu*m*h+x*sin(m*h));
    yn=sin(-nu*(m+1)*h+x*sin((m+1)*h));
    int1+=h/3*(yp+4*y+yn);
  }
  h=0.000004;
  sign=(((int)nu)%2==0?1:-1);
  yp=1+sign;
  y=exp(-x*sinh(h))*(exp(nu*h)+sign*exp(-nu*h));
  yn=exp(-x*sinh(2*h))*(exp(nu*2*h)+sign*exp(-nu*2*h));
  int2=h/3*(yp+4*y+yn);
  for(m=3; m<800000; m+=2){
    yp=yn;
    y=exp(-x*sinh(m*h))*(exp(nu*m*h)+sign*exp(-nu*m*h));
    yn=exp(-x*sinh((m+1)*h))*(exp(nu*(m+1)*h)+sign*exp(-nu*(m+1)*h));
    int2+=h/3*(yp+4*y+yn);
  }
  return int1/pi-int2/pi;  
}

double besselI(double nu,double x){
  double next,sum,next1,next2,sum1,sum2,int1,int2,yp,y,yn,h;
  int m,N;
  if((nu<0)&&(round(nu)==nu))
    if(((int)round(nu))%2==1) return -besselI(-nu,x);
      else return besselJ(-nu,x);
  if((nu==0.0)&&(x==0.0)) return 1.0;
  if(x==0.0) return 0.0;
 // if(nu>1) return 2*(nu-1)/x*besselI(nu-1,x)-besselI(nu-2,x);
 // if(nu<-1) return 2*(nu+1)/x*besselI(nu+1,x)-besselI(nu+2,x);
 // power series
  if(x<5){ 
    next=exp(log(x/2)*nu)/eulergamma(nu+1);
    sum=next;
    for(m=1; m<200; m++){
      next*=(x*x/4)/(m*(m+nu));
      sum+=next;
    }
    return sum;
  }
  // asymptotic expansion
  if((x>400)&&(x>(40*(nu*nu-0.25)))){
    next1=1;
    sum1=next1;
    //next2=(4*nu*nu-1)/(8*x);
    //sum2=next2;
    for(m=1; m<=200; m++){
      next1*=-(4*nu*nu-(2*m-1)*(2*m-1))/(m*8*x);
      sum1+=next1;
      //next2*=-(4*nu*nu-(4*m-1)*(4*m-1))*(4*nu*nu-(4*m+1)*(4*m+1))/(2*m*(2*m+1)*8*x*8*x);
      //sum2+=next2;        
    }
    return exp(x)/sqrt(2*pi*x)*sum1;
  }
  //integral
  N=(int)400*round(x);
  h=pi/N;
  yp=exp(x);
  y=exp(x*cos(h))*cos(nu*h);
  yn=exp(x*cos(2*h))*cos(nu*2*h);
  int1=h/3*(yp+4*y+yn);
  for(m=3; m<N; m+=2){
    yp=yn;
    y=exp(x*cos(m*h))*cos(nu*m*h);
    yn=exp(x*cos((m+1)*h))*cos(nu*(m+1)*h);
    int1+=h/3*(yp+4*y+yn);
  }
  if(nu==round(nu)) return int1/pi;
  h=0.000004;
  yp=1;
  y=exp(-x*cosh(h)-nu*h);
  yn=exp(-x*cosh(2*h)-nu*h);
  int2=h/3*(yp+4*y+yn);
  for(m=3; m<500000; m+=2){
    yp=yn;
    y=exp(-x*cosh(m*h)-nu*m*h);
    yn=exp(-x*cosh((m+1)*h)-nu*(m+1)*h);
    int2+=h/3*(yp+4*y+yn);
  }
  return int1/pi-sin(nu*pi)*int2/pi;
}

double besselK(double nu,double x){
  double next1,next2,sum1,sum2,int1,int2,yp,y,yn,h,max;
  int m,N,sign;
  //if(round(nu)!=nu) return (besselI(nu,x)-besselI(-nu,x))*pi/(2*sin(nu*pi));
  // asymptotic expansion
  if((x>400)&&(x>(40*(nu*nu-0.25)))){
    next1=1;
    sum1=next1;
    //next2=(4*nu*nu-1)/(8*x);
    //sum2=next2;
    for(m=1; m<=200; m++){
      next1*=(4*nu*nu-(2*m-1)*(2*m-1))/(m*8*x);
      sum1+=next1;
      //next2*=-(4*nu*nu-(4*m-1)*(4*m-1))*(4*nu*nu-(4*m+1)*(4*m+1))/(2*m*(2*m+1)*8*x*8*x);
      //sum2+=next2;        
    }
    return sqrt(pi/(2*x))*exp(x)*sum1;
  }
  //integral
  h=0.0004;
  yp=exp(-x);
  y=exp(-x*cosh(h))*cosh(nu*h);
  yn=exp(-x*cosh(2*h))*cosh(nu*2*h);
  int2=h/3*(yp+4*y+yn);
  max=(100000<2500*log(40/x)?2500*log(40/x):100000);
  for(m=3; m<max; m+=2){
    yp=yn;
    y=exp(-x*cosh(m*h))*cosh(nu*m*h);
    yn=exp(-x*cosh((m+1)*h))*cosh(nu*(m+1)*h);
    int2+=h/3*(yp+4*y+yn);
  }
  return int2;  
}

double airyAi(double x){
  if(x==0) return 0.35502805388781719;
  if(x>0){
    return 1/pi*sqrt(x/3)*besselK(1.0/3,2.0/3*x*sqrt(x));
  }else{
    return 1.0/3*sqrt(-x)*(besselJ(1.0/3,2.0/3*(-x)*sqrt(-x))+besselJ(-1.0/3,2.0/3*(-x)*sqrt(-x)));
  }
}

double airyBi(double x){
  if(x==0) return 0.61492662744600065;
  if(x>0){
    return sqrt(x/3)*(besselI(1.0/3,2.0/3*x*sqrt(x))+besselI(-1.0/3,2.0/3*x*sqrt(x)));
  }else{
    return sqrt(-x/3)*(besselJ(-1.0/3,2.0/3*(-x)*sqrt(-x))-besselJ(1.0/3,2.0/3*(-x)*sqrt(-x)));
  }
}

int main(int argc,char *argv[]){
  if(argc<2){
    printf("Usage:\ngamma x(double)\n");
    printf("beta x(double) y(double)\n");
    printf("besselJ nu(double) x(double)\n");
    printf("besselY nu(double) x(double)\n");
    printf("besselI nu(double) x(double)\n");
    printf("besselK nu(double) x(double)\n");
    printf("airyAi x(double)\n");
    printf("airyBi x(double)\n");
    return 0;
  }
  if(strcmp(argv[1],"gamma")==0){
    double x=atof(argv[2]);
    printf("gamma(%s) = %.17lf\n",argv[2],eulergamma(x));
  }
  if(strcmp(argv[1],"beta")==0){
    double x=atof(argv[2]);
    double y=atof(argv[3]);
    printf("beta(%s,%s) = %.17lf\n",argv[2],argv[3],eulerbeta(x,y));
  }
  if(strcmp(argv[1],"besselJ")==0){
    double nu=atof(argv[2]);
    double x=atof(argv[3]);
    printf("besselJ(%s,%s) = %.17lf\n",argv[2],argv[3],besselJ(nu,x));
  }
  if(strcmp(argv[1],"besselY")==0){
    double nu=atof(argv[2]);
    double x=atof(argv[3]);
    printf("besselY(%s,%s) = %.17lf\n",argv[2],argv[3],besselY(nu,x));
  }
  if(strcmp(argv[1],"besselI")==0){
    double nu=atof(argv[2]);
    double x=atof(argv[3]);
    printf("besselI(%s,%s) = %.17lf\n",argv[2],argv[3],besselI(nu,x));
  }
  if(strcmp(argv[1],"besselK")==0){
    double nu=atof(argv[2]);
    double x=atof(argv[3]);
    printf("besselK(%s,%s) = %.17lf\n",argv[2],argv[3],besselK(nu,x));
  }
  if(strcmp(argv[1],"airyAi")==0){
    double x=atof(argv[2]);
    printf("Ai(%s) = %.17lf\n",argv[2],airyAi(x));
  }
  if(strcmp(argv[1],"airyBi")==0){
    double x=atof(argv[2]);
    printf("Bi(%s) = %.17lf\n",argv[2],airyBi(x));
  }
  return 0;
}






