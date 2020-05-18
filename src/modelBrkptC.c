#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

static double min(double a, double b){
  if(a<b) return a;
  else return b;
}


double pmicTrue2(double m,double M1Test,double M2Test,double M1True,double M2True,double e){
  double prob=0;
  
  if (m<=M1True)  prob=pnorm(M1Test,m,e,TRUE,FALSE);
  if (m>M1True && m<M2True)  prob=pnorm(M2Test-1,m,e,TRUE,FALSE)-pnorm(M1Test,m,e,TRUE,FALSE);
  if (m>=M2True)  prob=1-pnorm(M2Test-1,m,e,TRUE,FALSE);
  return (prob);
}

double pdia2(double d,double m,double D1,double D2,double M1,double M2,double t){
  double prob=0;

  if (m<=M1)  prob=1-pnorm(D2-.5,d,t,TRUE,FALSE);
  if (m>M1 && m<M2)  prob=pnorm(D2-.5,d,t,TRUE,FALSE)-pnorm(D1+.5,d,t,TRUE,FALSE);
  if (m>=M2)  prob=pnorm(D1+.5,d,t,TRUE,FALSE);
  
  return(prob);
}

double calcLossTrueMin2(double D1,double D2,double *gridx,double *weights,double *f,double M1Test,double M2Test,
  double M1True, double M2True,double e,double t,int lgrid){

  double temp,sum=0;
  int i;
  
  for(i=0; i<lgrid; i++){
    temp=weights[i]*pow(min(0,pdia2(f[i],gridx[i],D1,D2,M1True,M2True,t)-pmicTrue2(gridx[i],M1Test,M2Test,M1True,M2True,e)),2);
    sum+=temp;
  }
  return(sum);
}


void findDIATrue(double *gridx,double *weights,double *fit,double *M1Test, double *M2Test,double *M1True,double *M2True,
  double *xsig,double *ysig,double *minDIA, double *maxDIA,int *lgrid,double *D1,double *D2, double *index, int *minWidth, int *maxWidth){
	
  double test,max=9999;
  int x,y;
  
  for(x=*minDIA; x<*maxDIA; x++){
    for(y=x+1; y<*maxDIA; y++){
      if((y-x)>=*minWidth & y-x<=*maxWidth){
        test=calcLossTrueMin2(x,y,gridx,weights,fit,*M1Test,*M2Test,*M1True,*M2True,*xsig,*ysig,*lgrid);
        if(test<max){max=test;*D1=x;*D2=y;*index=test;};
      }
    }
  }
}



