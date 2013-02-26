#include <math.h>
#include <R.h>
#include "headers.h"
#define MATRIX(x,j,k) x[*D*(j)+(k)]
#define damp 1.0e0
#define LIM 2.0e1

void log_like_forces(int *directed, int *N, int *D, int *steps, double *Y, double *X, double *B, double *m){
  double *dxdy=calloc(*N* *D, sizeof(double)); // position change vector
  double xd[*D], ded, f;
  int i, j, k, d, s;
  double coolexp = 1.5;
  double t;
  if (*steps > 0) // so that we can skip it if we like
  for(i=*steps;i>0;i--) {
    t=*m*pow(i/(double)*steps,coolexp);
    /* calculate deltas for each undirected pair */
    for(j=0;j<*N;j++) {
    R_CheckUserInterrupt();
    s=(*directed ? 0 : j+1);
    for(k=s;k<*N;k++) 
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
        f=-Y[j* *N+k]*(*B-ded)-log(1.0+exp(*B-ded));
        //f=-Y[j* *N+k]*(*B+1.0/ded)-log(1.0+exp(*B+1.0/ded)); // NEW dists
  	for (d=0;d<*D;d++)
          {
          MATRIX(dxdy, j, d)-=xd[d]*damp*(f); /* Add to the position change vector */
          MATRIX(dxdy, k, d)+=xd[d]*damp*(f);
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
	/* only needed when including isolates */
	/*
	ded=0.0;
        for (d=0;d<*D;d++)
	  ded += MATRIX(X, j, d)*MATRIX(X, j, d);
	ded = sqrt(ded);
	if (ded > LIM)  // don't let them move too far (isolated nodes do this)
          for (d=0;d<*D;d++)
	    MATRIX(X, j, d) = MATRIX(X, j, d)*LIM/ded;
	*/
    } // close nodes iterations
  } // close iterations loop
  free(dxdy);
  return;
}

