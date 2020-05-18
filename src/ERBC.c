#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

double findERBBrkpts(int D1,int D2,double *MIC, double *DIA,double MICBrkpt1,double MICBrkpt2, int N,
	double VM1,double M1,double m1,double VM2,double M2,double m2){

  int i,II=0,SS=0,VM=0,RR=0,M=0,m=0;
  double index;

  for(i=0; i<N; i++){
    if(MIC[i]<=MICBrkpt1 & DIA[i]>=D2)SS=SS+1;
    else if ((MIC[i]>MICBrkpt1 & MIC[i]<MICBrkpt2) & (DIA[i]>D1 & DIA[i]<D2)) II=II+1;
    else if(MIC[i]>=MICBrkpt2 & DIA[i]<=D1)RR=RR+1;
    else if((MIC[i]>MICBrkpt1 & MIC[i]<MICBrkpt2) & DIA[i]>=D2)m=m+1;
    else if(MIC[i]>=MICBrkpt2 & DIA[i]>=D2)VM=VM+1;
    else if(MIC[i]<=MICBrkpt1 & (DIA[i]>D1 & DIA[i]<D2))m=m+1;
    else if(MIC[i]>=MICBrkpt2 & (DIA[i]>D1 & DIA[i]<D2))m=m+1;
    else if(MIC[i]<=MICBrkpt1 & DIA[i]<=D1)M=M+1;
    else if((MIC[i]>MICBrkpt1 & MIC[i]<MICBrkpt2) & DIA[i]<=D1)m=m+1;
  }


  index=VM1*VM/N+M1*M/N+m1*m/N;

  return(index);
}

void ERB(int *minDIA, int *maxDIA, double *MICWithinTwo, double *MICOutsideTwo, double *DIAWithinTwo, double *DIAOutsideTwo, 
	double *MICBrkpt1, double *MICBrkpt2, int *N1, int *N2, double *VM1, double *M1, double *m1, 
	double *VM2, double *M2, double *m2,int *minWidth, int *maxWidth, int *D1, int *D2){

  int x, y;
  double minVal=9999;
  double idx1, idx2;
  *D1=0; *D2=1; 
  
  for(x=*minDIA; x<(*maxDIA-2); x++){
    for(y=(x+1); y<*maxDIA; y++){
      if((y-x)>=*minWidth & (y-x)<=*maxWidth){
        if(*N1>0)
          idx1=findERBBrkpts(x,y,MICWithinTwo,DIAWithinTwo,*MICBrkpt1,*MICBrkpt2,*N1,*VM1,*M1,*m1,*VM2,*M2,*m2);
        else idx1=0;
        if(*N2>0)
	  idx2=findERBBrkpts(x,y,MICOutsideTwo,DIAOutsideTwo,*MICBrkpt1,*MICBrkpt2,*N2,*VM1,*M1,*m1,*VM2,*M2,*m2);
        else idx2=0;
//        Rprintf("%d %d %f %f %d %d \n", x,y,idx1+idx2,minVal,*minDIA,*maxDIA); 
        if(idx1+idx2<=minVal)
          if(idx1+idx2==minVal){
            if(y-x>=*D2-*D1) {minVal=idx1+idx2; *D1=x; *D2=y;}
          }
          else {minVal=idx1+idx2; *D1=x; *D2=y;}
      } 
    }
  }

}

