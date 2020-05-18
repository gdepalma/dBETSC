#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <time.h>
#include <Rmath.h>


double max_array(double *a, int num_elements){

   int i;
   double max=-9999999;
   for (i=0; i<num_elements; i++)
     if (a[i]>max)
       max=a[i];

   return(max);
}

double max_arrayInt(int *a, int num_elements){

   int i;
   double max=-9999999;
   for (i=0; i<num_elements; i++)
     if (a[i]>max)
       max=a[i];

   return(max);
}

//sample one value from given probabilities
int sampleOne(int *a, double *probs, int n){
  
  int i;
  double cumProbs[n],ran;
  
  for(i=0; i<n; i++)
    if(i==0) cumProbs[i]=probs[i];
    else cumProbs[i]=probs[i]+cumProbs[i-1];
   
  //rescale probs
  for(i=0; i<n; i++){
    cumProbs[i]/=cumProbs[n-1];
//    Rprintf("%f \n", cumProbs[i]); 
  }

  ran = (double)rand()/(double)RAND_MAX;
  
  for(i=0; i<n; i++)
    if(ran<cumProbs[i]) return(a[i]);
  
  return(0);
}

double t_logpdf(double x, double mu, double v, double df){

  double temp, temp1, pi;
  pi = 4.0*atan(1.0);
  temp=lgamma(.5 * df+.5)-lgamma(.5 * df)-log(sqrt(df * pi * v));
  temp1=temp-.5 * (df+1) * log(1+(1/df)*pow((x-mu),2)/v);

  return(temp1);
}
  

void normalize_logprob(double *log_prob, int active, double *probs){

  int i;
  double max,temp=0,g;

  max=max_array(log_prob,active);
  
  for(i=0; i<active; i++)
    temp=temp+exp(log_prob[i]-max);
  
  g=log(temp)+max;
  for(i=0; i<active; i++)
    probs[i]=exp(log_prob[i]-g);
}

double loglik(double xn, double *x, double tau0, double beta0, double mu0, double kappa0,int N){

  double kappaN, tauN, betaN, muN,temp=0,xm,ssx=0;
  int i;

  kappaN=kappa0+N;
  tauN=tau0+N/2;
  if(N>0){
    for(i=0; i<N; i++)
      temp=temp+x[i];
    xm=temp/N;
    for(i=0; i<N; i++)
      ssx=ssx+pow(x[i]-xm,2);
    betaN = beta0 + 0.5 * ssx + (kappa0 * N *pow(xm - mu0,2)) / (2 * kappaN);
    muN = (kappa0 * mu0 + N * xm) / kappaN;
  }else{
    betaN=beta0;
    muN=mu0;
  }
 
  return(t_logpdf(xn, muN, betaN * (kappaN + 1) / (tauN * kappaN), 2 * kappaN));
}

void normPosterior(double *xgrid, double *x, double tau0, double beta0, double mu0, double kappa0, int N,double *tempNorm){
  
  double kappaN, tauN, betaN, muN,mu1,sig1,temp=0,xm,ssx=0;
  int i;

  kappaN=kappa0+N;
  tauN=tau0+N/2;

  if(N>0){
    for(i=0; i<N; i++){
      temp=temp+x[i];
    }
    xm=temp/N;
    for(i=0; i<N; i++)
      ssx=ssx+pow(x[i]-xm,2);
    betaN = beta0 + 0.5 * ssx + (kappa0 * N *pow(xm - mu0,2)) / (2 * kappaN);
    muN = (kappa0 * mu0 + N * xm) / kappaN;
  }else{
    betaN=beta0;
    muN=mu0;
  }

  mu1=rnorm(muN,1/kappaN);
  sig1=sqrt(1/rgamma(tauN,1/betaN));

  for(i=0; i<1200; i++){
   tempNorm[i]=dnorm(xgrid[i],mu1,sig1,FALSE);
  }
}


void DP(double *x, double *alpha, double *tau0, double *beta0, double *mu0, double *kappa0,
  int *C, int *N, int *m,int *nIter, int *NM, double *xgrid, double *densPosterior, int *burnin){


  int i,j,k,k_active,numComp,idx,idx2,idxPost=0,n=0;
  double temp;
  double cn[*N-1];
  double tempNorm[1200];
  double tempPost[1200];

  //loop over customers
  for(n=0; n < *N; n++){

	  Rprintf("\n\nCustomer: %d \n",n);
  
	  //remove customer, cn new customer array
	  for(i=0; i < *N-1; i++)
	    if(i < n)
	      cn[i]=C[i];
	    else
	      cn[i]=C[i+1];

	  //remove customer from membership
	  m[C[n]]=m[C[n]]-1;
	  

	  //active dishes + 1 new dish
	  k_active=0;
	  double kset[*NM];
	  for(i=0;i<*NM;i++)
	    if(m[i]>0){
	      kset[k_active]=i;
	      k_active++;
	    }
	  kset[k_active+1]=k_active+1;
	  k_active=k_active+1;
	  //reduce kset vector
	  int ksetRed[k_active];
	  for(i=0;i<k_active;i++){
	    ksetRed[i]=kset[i];
	  }
	  ksetRed[k_active-1]=max_arrayInt(ksetRed,k_active)+1;
	  Rprintf("Number active + 1: %d \n",k_active);
	  for(i=0;i<k_active;i++)
	    Rprintf("Active plus new active: %d \n", ksetRed[i]);


	  double prior[k_active];
	  double log_lik[k_active];
	  double post[k_active];
	  idx=0;
	  do{ 
	    //Rprintf("ksetRed %d \n",ksetRed[idx]);
	    //prior
	    if (m[ksetRed[idx]] > 0){
	      //prior for old dish
	      prior[idx] = log(m[ksetRed[idx]]);
	    }else{
	      //prior for new dish
	      prior[idx] = log(*alpha);
	    }
	//    Rprintf("%f \n",prior[idx]);
	    //log lik
	    //get obs in clusters - xtemp
	    double xtemp[m[ksetRed[idx]]];
	    idx2=0;
	    for(i=0; i<*N; i++)
	      if(C[i]==ksetRed[idx]){
		xtemp[idx2]=x[i];
		idx2++;
	      }
	//     Rprintf("%d \n",idx2);
	     log_lik[idx] = loglik(x[n],xtemp, *tau0, *beta0, *mu0, *kappa0,idx2);
	//     Rprintf("%f \n",log_lik[idx]);

	     //posterior
	     post[idx]=log_lik[idx]+prior[idx];
//	     Rprintf("%f \n",post[idx]);
	     idx=idx+1; 

	  }while(idx<k_active);

	  //normalize
	  double normPost[k_active];
	  normalize_logprob(post, k_active, normPost);

	  //update cluster assignment
	  C[n]=sampleOne(ksetRed,normPost, k_active);
	  Rprintf("%d %d \n",n,C[n]);

	  //insert customer into cluster
	  m[C[n]]=m[C[n]]+1;

	  for(i=0;i<*N;i++)
	    Rprintf("%d %d ", i,m[i]);

	  Rprintf("\n");

//	  for(i=0;i<*N;i++)
//	    Rprintf("%d %d ", i,C[i]);

	  Rprintf("\n");

	  for(i=0;i<k_active;i++){
	    Rprintf("Active plus new active: %d \n", ksetRed[i]);
          }
          Rprintf("\n");
	  for(i=0;i<k_active-1;i++){
	    Rprintf("M %d \n", m[ksetRed[i]]);
          }   

  } //end  loop over customers
 

  //calculate posterior
  k_active=0;
  double kset[*NM];
  for(i=0;i<*NM;i++)
    if(m[i]>0){
      kset[k_active]=i;
      k_active++;
    }
  kset[k_active]=k_active;
  //reduce kset vector
  int ksetRed[k_active];
  for(i=0;i<k_active;i++){
    ksetRed[i]=kset[i];
  }
  ksetRed[k_active-1]=max_arrayInt(ksetRed,k_active);
  numComp=0; i=0;
  do{ 
	  if (m[i] >0){
	    Rprintf("\n i %d \n", i);
	    double xtemp[m[i]];
	    idx2=0;
	    for(j=0; j<*N; j++)
	      if(C[j]==i){
		      xtemp[idx2]=x[j];
    		idx2++;
	    }
	    //call norm prob
	    Rprintf("m[i] %d \n", m[i]);
	    Rprintf("index idx2 %d %d \n", i,idx2);
	    Rprintf("N %d \n", *N);
	    Rprintf("prob %f \n", (float)idx2/ (float)*N);
	    normPosterior(xgrid,xtemp,*tau0,*beta0,*mu0,*kappa0,idx2,tempNorm);
	    for(j=0; j<1200; j++){
	      tempPost[j]=((float)idx2/ (float)*N)*tempNorm[j];
	      //Rprintf("%d %f \n", j,densPosterior[j]);
	    }
	    numComp++;
	  }
    i++;
  }while(numComp<k_active);
  
  for(i=0; i<1200; i++){
    densPosterior[idxPost]=tempPost[i];
    idxPost++;
  }
  
  Rprintf("idxPost %d \n", idxPost);

} //end function

void runDP(double *x, double *alpha, double *tau0, double *beta0, double *mu0, double *kappa0,
  int *C, int *N, int *m, int *nIter,int *NM, double *xgrid, double *densPosterior, int *burnin){
  
  int i,j,k,k_active,numComp,idx,idx2,idxPost=0,n=0;
  double temp;
  double cn[*N-1];
  double tempNorm[1200];
  double tempPost[1200];
  
  for(int niter=0; niter<*nIter; niter++){

    //loop over customers
    for(n=0; n < *N; n++){
      
      for(i=0; i<1200; i++)
        tempPost[i]=0;
    
  	  //remove customer, cn new customer array
  	  for(i=0; i < *N-1; i++)
  	    if(i < n)
  	      cn[i]=C[i];
  	    else
  	      cn[i]=C[i+1];
  
  	  //remove customer from membership
  	  m[C[n]]=m[C[n]]-1;
  	  
  
  	  //active dishes + 1 new dish
  	  k_active=0;
  	  double kset[*NM];
  	  for(i=0;i<*NM;i++)
  	    if(m[i]>0){
  	      kset[k_active]=i;
  	      k_active++;
  	    }
  	  kset[k_active+1]=k_active+1;
  	  k_active=k_active+1;
  	  //reduce kset vector
  	  int ksetRed[k_active];
  	  for(i=0;i<k_active;i++){
  	    ksetRed[i]=kset[i];
  	  }
  	  ksetRed[k_active-1]=max_arrayInt(ksetRed,k_active)+1;
  
  	  double prior[k_active];
  	  double log_lik[k_active];
  	  double post[k_active];
  	  idx=0;
  	  do{ 
  	    //prior
  	    if (m[ksetRed[idx]] > 0){
  	      //prior for old dish
  	      prior[idx] = log(m[ksetRed[idx]]);
  	    }else{
  	      //prior for new dish
  	      prior[idx] = log(*alpha);
  	    }
  	    //log lik
  	    //get obs in clusters - xtemp
  	    double xtemp[m[ksetRed[idx]]];
  	    idx2=0;
  	    for(i=0; i<*N; i++)
          if(C[i]==ksetRed[idx]){
            xtemp[idx2]=x[i];
            idx2++;
          }
  	     log_lik[idx] = loglik(x[n],xtemp, *tau0, *beta0, *mu0, *kappa0,idx2);
  
  	     //posterior
  	     post[idx]=log_lik[idx]+prior[idx];
  	     idx=idx+1; 
  
  	  }while(idx<k_active);
  
  	  //normalize
  	  double normPost[k_active];
  	  normalize_logprob(post, k_active, normPost);
  
  	  //update cluster assignment
  	  C[n]=sampleOne(ksetRed,normPost, k_active);
  
  	  //insert customer into cluster
  	  m[C[n]]=m[C[n]]+1;
  
    } //end  loop over customers
   
   
    if(niter>=*burnin){
//      Rprintf("iteration %d \n", niter);
  
      //calculate posterior after burnin
      k_active=0;
      double kset[*NM];
      for(i=0;i<*NM;i++)
        if(m[i]>0){
          kset[k_active]=i;
          k_active++;
        }
      kset[k_active]=k_active;
      //reduce kset vector
      int ksetRed[k_active];
      for(i=0;i<k_active;i++){
        ksetRed[i]=kset[i];
      }
      ksetRed[k_active-1]=max_arrayInt(ksetRed,k_active);
      numComp=0; i=0;
      do{ 
    	  if (m[i] >0){
    	    double xtemp[m[i]];
    	    idx2=0;
    	    for(j=0; j<*N; j++)
    	      if(C[j]==i){
    		      xtemp[idx2]=x[j];
        		idx2++;
    	    }
    	    //call norm prob
    	    normPosterior(xgrid,xtemp,*tau0,*beta0,*mu0,*kappa0,idx2,tempNorm);
    	    for(j=0; j<1200; j++){
    	      tempPost[j]=tempPost[j]+((float)idx2/ (float)*N)*tempNorm[j];
    	    }
    	    numComp++;
    	  }
        i++;
      }while(numComp<k_active);
      
      for(i=0; i<1200; i++){
        densPosterior[idxPost]=tempPost[i];
        idxPost++;
      }
      
//      Rprintf("idxPost %d \n", idxPost);
      
    } // end burnin
    
  } //end iterations

} //end function
  
