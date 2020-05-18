#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

static double min(double a, double b){
  if(a<b) return a;
  else return b;
}


double pmicTrue1(double m,double MTest,double MTrue,double e){
  double prob=0;
  
  if (m<=MTrue)  prob=pnorm(MTest,m,e,TRUE,FALSE);
  if (m>MTrue)  prob=1-pnorm(MTest,m,e,TRUE,FALSE);
  return (prob);
}

double pdia1(double d,double m,double DIABrkpt,double MTrue,double t){
  double prob=0;

  if (m<=MTrue)  prob=1-pnorm(DIABrkpt-.5,d,t,TRUE,FALSE);
  if (m>MTrue)  prob=pnorm(DIABrkpt-.5,d,t,TRUE,FALSE);
  
  return(prob);
}

double calcLossTrueMin1(double DIABrkpt,double *gridx,double *weights,double *f,double MTest,double MTrue,double e,double t,int lgrid){

  double temp,sum=0;
  int i;
  
  for(i=0; i<lgrid; i++){
    temp=weights[i]*pow(min(0,pdia1(f[i],gridx[i],DIABrkpt,MTrue,t)-pmicTrue1(gridx[i],MTest,MTrue,e)),2);
    sum+=temp;
  }
  return(sum);
}


void findDIATrueOne(double *gridx,double *weights,double *fit,double *MTest, double *MTrue,
  double *xsig,double *ysig,double *minDIA, double *maxDIA,int *lgrid,double *DIABrkpt, double *index){
	
  double x,test,max=9999;
  
  for(x=*minDIA; x<=*maxDIA; x++){
    test=calcLossTrueMin1(x,gridx,weights,fit,*MTest,*MTrue,*xsig,*ysig,*lgrid);
    if(test<=max){max=test;*DIABrkpt=x;*index=test;};
  }
 }



