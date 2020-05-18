#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

double findERBBrkptsOne(double DIABrkpt,double *MIC, double *DIA,double MICBrkpt,int N,double VM,double M){

  int i,SS=0,numVM=0,RR=0,numM=0;
  double index;

  for(i=0; i<N; i++){
    if(MIC[i]<=MICBrkpt & DIA[i]>=DIABrkpt) SS=SS+1;
    else if(MIC[i]>=MICBrkpt & DIA[i]<=DIABrkpt) RR=RR+1;
    else if(MIC[i]>=MICBrkpt & DIA[i]>=DIABrkpt) numVM=numVM+1;
    else if(MIC[i]<=MICBrkpt & DIA[i]<=DIABrkpt) numM=numM+1;
  }


  index=VM*numVM/N+M*numM/N;

  return(index);
}

void ERBOne(double *MIC, double *DIA, double *MICBrkpt, int *N, double *VM, double *M, double *DIABrkpt,
	double *minDIA, double *maxDIA, double *index){

  double x,test,max=9999;

  for(x=*minDIA; x<=*maxDIA; x++){
	test=findERBBrkptsOne(x,MIC,DIA,*MICBrkpt,*N,*VM,*M);
	if(test<=max){max=test;*DIABrkpt=x;*index=test;};
  }

}

