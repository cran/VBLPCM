#include <math.h>
#include <R.h>
#include "headers.h"  

void Y_to_E (int *N, int *directed, double *Y, int *E)
  {
  int j1, j2, i=0, s;
  for (j1=0;j1<*N;j1++)
    {
    s=(*directed==1) ? 0 : j1+1; 
    for (j2=s;j2<*N;j2++)
      if (Y[j1* *N+j2]>0.0)
        {
	E[i*2]=j1+1;
	E[i*2+1]=j2+1;
	i++;
	}
    } 
  return;
  }

void Y_to_nonE (int *N, int *directed, double *Y, int *nonE)
  {
  int j1, j2, i=0, s;
  for (j1=0;j1<*N;j1++)
    {
    s=(*directed==1) ? 0 : j1;
    for (j2=s;j2<*N;j2++)
      if (Y[j1* *N+j2]==0.0)
        {
	nonE[i*2]=j1+1;
	nonE[i*2+1]=j2+1;
	i++;
	}
    }
  return;
  }

void Y_to_M (int *N, int *directed, double *Y, int *M)
  {
  int j1, j2, i=0, s;
  for (j1=0;j1<*N;j1++)
    {
    s=(*directed==1) ? 0 : j1+1;
    for (j2=s;j2<*N;j2++)
      if (isnan(Y[j1* *N+j2]))
        {
        M[i*2]=j1+1;
        M[i*2+1]=j2+1;
	i++;
        }
    }
  return;
  }
void E_to_Y (int *N, int *NE, int *directed, int *E, double *Y)
  {
  int i;
  for (i=0;i<*NE;i++)
    {
    Y[(E[i*2]-1) * *N + E[i*2+1]-1]=1.0;
    if (*directed==0)
      Y[(E[i*2+1]-1) * *N + E[i*2]-1]=1.0;
    }
  return;
  }
