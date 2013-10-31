#include <math.h>
#include <R.h>
#include "headers.h"
#define MATRIX(x,j,k) x[*D*(j)+(k)]

void fruchterman_reingold(int *directed, int *N, int *D, int *steps, double *Y, double *X, 
                          double *repulserad, double *m, double *volume){
  double *dxdy=calloc(*N* *D, sizeof(double)); // position change vector
  double xd[*D], ded, af, rf, f;
  int i, j, k, d;
  double coolexp=1.5, t;
  double frk=pow(*volume/(double)*N, 1.0/(double)*D);
  for(i=*steps;i>0;i--) {
    t=*m*pow(i/(double)*steps,coolexp);
    /* calculate deltas for each undirected pair */
    for(j=0;j<*N;j++) {
    R_CheckUserInterrupt();
    for(k=0;k<*N;k++) 
      if ((k!=j) & (!isnan(Y[j* *N+k])))
	{
        /* Obtain difference vector */
	ded = 0.0;
        for (d=0;d<*D;d++)
          {
          xd[d]=MATRIX(X, j, d)-MATRIX(X, k, d);
          ded+=xd[d]*xd[d];
          }
        ded=sqrt(ded);
	for (d=0;d<*D;d++)
          xd[d]/=ded;                /* Rescale differences to length 1 */
	/*Calculate repulsive "force"*/
        rf=0.5*frk*frk*(1.0/ded-ded*ded/ *repulserad);
        af=Y[j* *N+k]*ded*ded/frk; // more attraction possible than standard FR routines
	f=af-rf;
  	for (d=0;d<*D;d++)
          {
          MATRIX(dxdy, j, d)-=xd[d]*f; /* Add to the position change vector */
          MATRIX(dxdy, k, d)+=xd[d]*f;
          }
        }
    }
    /* Dampen motion, if needed, and move the points */   
    for(j=0;j<*N;j++){
      ded = 0.0;
        for (d=0;d<*D;d++)
          ded += MATRIX(dxdy, j, d)*MATRIX(dxdy, j, d);
        ded=sqrt(ded);
        if(ded>t)
          {
          ded=t/ded;
          for (d=0;d<*D;d++)
            MATRIX(dxdy, j, d)*=ded;
          }
        for (d=0;d<*D;d++)
	//  if (!isnan(MATRIX(dxdy, j, d)))
	    MATRIX(X, j, d)+=MATRIX(dxdy, j, d); // Update positions 
    } // close nodes iterations
  } // close iterations loop
  free(dxdy);
  return;
}

