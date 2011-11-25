#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sf.h>
#include "headers.h"
#include <R_ext/Utils.h>


#define JP *params->P 


double diff_max(double *vec1, double *vec2, int n) // compute largest absolute difference b/w two vectors
  {
  int i;
  double diff = 0.0;
  for (i=0; i<n; i++)
    diff = (diff>fabs(vec1[i]-vec2[i]) ? diff:  fabs(vec1[i]-vec2[i]));
  return diff;
  }

double diff_mean(double *vec1, double *vec2, int n) // compute largest absolute difference b/w two vectors
  {
  int i;
  double diff = 0.0;
  for (i=0; i<n; i++)
    diff += fabs(vec1[i]-vec2[i]);
  return diff/n;
  }

void bb(double *lim, double *tol); // in bb.c
void optim(); 

void Rf_VB_bbs(int *steps,
  int *max_iter,
  int *P,
  int *D,
  int *N, 
  int *NE, // #edges
  int *NnonE, // #non-edges
  int *NM, //#missing edges
  int *G,
  double *Y,
  int *E,
  int *nonE,
  int *M,
  int *numedges,
  int *EnonE,
  int *diam,
  int *hopslist,
  double *XX,
  double *V_xi,
  double *V_psi2,
  double *V_z,
  double *V_sigma2,
  double *V_eta,
  double *V_lambda,
  double *V_omega2,
  double *V_nu,
  double *V_alpha,
  double *xi,
  double *psi2,
  double *sigma2,
  double *omega2,
  double *nu,
  double *alpha,
  double *inv_sigma02,
  double *tol,
  double *STRAT,
  double *seed,
  double *d_vector,
  int *conv)
    {
    double lim[2]={-1.0e2,1.0e2};
    int i, j, d, g, p, l, g1;
    double tmp, tmpsum1, tmpsum2;
    double mu_nought = 0.0;
  params=calloc(1,sizeof(Rf_params));
  params->MAX_ITER=max_iter;
  params->P=P;
  params->p=&p;
  params->D=D;
  params->d=&d;
  params->N=N;
  params->i=&i;
  params->NE=NE;
  params->NM=NM;
  params->NnonE=NnonE;
  params->G=G;
  params->Y=Y;
  params->E=E; // edges matrix. NE X 2
  params->nonE=nonE; // non-edges matrix. NM X 2
  params->M=M; // missing-edges matrix. NM X 2
  params->numedges=numedges;
  params->EnonE=EnonE;
  params->diam=diam; 
  params->hopslist=hopslist; 
  params->XX=XX; // design matrix for covariates. May also be used for sender / receiver effects, etc.
  params->V_xi=V_xi;
  params->V_psi2=V_psi2;
  params->V_z=V_z;
  params->V_sigma2=V_sigma2;
  params->V_eta=V_eta;
  params->V_lambda=V_lambda;
  params->V_omega2=V_omega2;
  params->V_nu=V_nu;
  params->V_alpha=V_alpha;
  params->xi=xi;
  params->psi2=psi2;
  params->sigma2=sigma2;
  params->omega2=omega2;
  params->nu=nu;
  params->alpha=alpha;
  params->inv_sigma02=inv_sigma02;
  params->dists=calloc(*N* *N, sizeof(double));
  params->STRAT=STRAT;
  params->seed=seed;
  params->conv=conv;
  int *samp_nodes = calloc(*N, sizeof(int));
  int *samp_groups= calloc(*G, sizeof(int));
  int *samp_coeffs= calloc(JP, sizeof(int));
  double *old_xi, *old_psi2, *old_z, *old_sigma2, *old_alpha,
         *old_nu, *old_eta, *old_lambda, *old_omega2;
  old_xi = calloc(*P, sizeof(double));
  old_psi2 = calloc(*P, sizeof(double));
  old_z = calloc(*N * *D, sizeof(double));
  old_sigma2 = calloc(*N, sizeof(double));
  old_alpha = calloc(*G, sizeof(double));
  old_nu = calloc(*G, sizeof(double));
  old_eta = calloc(*G * *D, sizeof(double));
  old_omega2 = calloc(*G, sizeof(double));
  old_lambda = calloc(*G * *N, sizeof(double));

  int flag;
  params->flag=&flag; 
  /*
  char **V_names = calloc(9, sizeof(char *));
  for (i=0; i<9; i++)
    V_names[i] = calloc(10, sizeof(char *));
  sprintf(V_names[0], "%s", "z");
  sprintf(V_names[1], "%s", "sigma2");
  sprintf(V_names[2], "%s", "lambda");
  sprintf(V_names[3], "%s", "eta");
  sprintf(V_names[4], "%s", "omega2");
  sprintf(V_names[5], "%s", "alpha");
  sprintf(V_names[6], "%s", "nu");
  sprintf(V_names[7], "%s", "xi");
  sprintf(V_names[8], "%s", "psi2");
  */

  tmp=0.0;
  
  for (l=0; l<*steps; l++) // number of cycles through the variational algorithm
  {
  R_CheckUserInterrupt();
  if (d_vector[0] > *tol)
    {
    flag=1;
    memcpy(old_z, V_z, *N * *D *sizeof(double));
    sample_permutation(*N, samp_nodes, params->seed);
    //#pragma omp parallel for
    for (i=0; i<*N; i++)
      { 
      params->i=&samp_nodes[i]; 
      optim();
      }
    d_vector[0]=diff_max(V_z, old_z, *N * *D);
    } else d_vector[0] = 0.0;
  R_CheckUserInterrupt();
  if (d_vector[3] > *tol)
    {
    memcpy(old_eta, V_eta, *G*sizeof(double));
    sample_permutation(*G, samp_groups, params->seed);
    for (g=0; g<*G; g++)
      {
      params->g=&samp_groups[g]; 
       for (d=0; d<*D; d++)
        {
        tmpsum1 = 0.0;
        tmpsum2 = 0.0;
        for (i=0; i<*N; i++)
          { 
          tmpsum1 += 0.5*V_lambda[samp_groups[g] * *N + i]* *inv_sigma02*V_alpha[samp_groups[g]]*V_z[i * *D + d];
          tmpsum2 += V_lambda[samp_groups[g] * *N + i]*0.5* *inv_sigma02*V_alpha[samp_groups[g]];
          }
        V_eta[samp_groups[g]* *D+d] = (tmpsum1 + 0.5*mu_nought/ *omega2)/(tmpsum2+0.5/ *omega2);
        }
      }
    d_vector[3]=diff_max(V_eta, old_eta, *G);
    } else d_vector[3] = 0.0;
  R_CheckUserInterrupt();
  if (d_vector[1] > *tol)
    {
    flag=2;
    memcpy(old_sigma2, V_sigma2, *N*sizeof(double));
    sample_permutation(*N, samp_nodes, params->seed);
    for (i=0; i<*N; i++)
      {
      params->i=&samp_nodes[i];
      lim[0]=1.0e-8; lim[1]=1.0e1;
      bb(lim, tol);
      }
    d_vector[1]=diff_max(V_sigma2, old_sigma2, *N);
    } else d_vector[1] = 0.0;
  // compute Euclidean distances + variances terms
  for (i=0;i<*N;i++)
    for (j=0;j<*N;j++) 
       {
       params->dists[i* *N + j] = 0.0;
       for (d = 0; d < *params->D; d++)
         params->dists[i* *N + j] += pow (V_z[i* *D + d] - V_z[j* *D + d], 2.0);
        params->dists[i* *N + j] = SQRT (params->dists[i* *N + j] +  *D*(params->V_sigma2[i] + params->V_sigma2[j]));
	}
  R_CheckUserInterrupt();
  if (d_vector[2] > *tol)
    {
    memcpy(old_lambda, V_lambda, *G * *N * sizeof(double));
    sample_permutation(*N, samp_nodes, params->seed);
    tmpsum1 = 0.0;
    for (g = 0; g < *G; g++)
      tmpsum1 += V_nu[g];
    for (i=0; i<*N; i++)
      {
      sample_permutation(*G, samp_groups, params->seed);
      for (g=0; g<*G; g++)
        {
        tmp = 0.0;
        for (d = 0; d < *D; d++)
          tmp += pow (V_z[samp_nodes[i] * *D + d] - V_eta[samp_groups[g] * *D + d], 2.0);
        //tmp = GSQRT (tmp);
        //tmp = GSQRT (tmp+*D*(V_sigma2[samp_nodes[i]]+V_omega2[samp_groups[g]]));
        tmp = GSQRT (tmp+(V_sigma2[samp_nodes[i]]+V_omega2[samp_groups[g]]));
        V_lambda[samp_groups[g]* *N + samp_nodes[i]] = 
                exp(-1.0-0.5* *inv_sigma02*V_alpha[samp_groups[g]]*(tmp)+
                       gsl_sf_psi(V_nu[samp_groups[g]]) - gsl_sf_psi(tmpsum1));
        }
      tmp = 0.0;
      for (g1=0; g1<*G; g1++)
        tmp += V_lambda[g1* *N + samp_nodes[i]];
      tmp = 1.0/tmp;
      for (g1=0; g1<*G; g1++)
        V_lambda[g1* *N + samp_nodes[i]] = V_lambda[g1* *N + samp_nodes[i]]*tmp;
      }
    d_vector[2]=diff_max(V_lambda, old_lambda, *N * *G);
    } else d_vector[2] = 0.0;
  
  R_CheckUserInterrupt();
  if (d_vector[4] > *tol)
    {
    memcpy(old_omega2, V_omega2, *G*sizeof(double));
    sample_permutation(*G, samp_groups, params->seed);
      for (g=0; g<*G; g++)
        {
        tmp = 0.0;
        for (i=0; i<*N; i++)
          tmp += V_lambda[samp_groups[g]* *N + i]* *inv_sigma02*V_alpha[samp_groups[g]];
	tmp = tmp/(*D/ *inv_sigma02) + 1.0/ *omega2;
        V_omega2[samp_groups[g]] = 1.0/tmp;
        }
    d_vector[4]=diff_max(V_omega2, old_omega2, *G);
      } else d_vector[4] = 0.0;
  R_CheckUserInterrupt();
  if (d_vector[5] > *tol)
    {
    flag=3;
    memcpy(old_alpha, V_alpha, *G*sizeof(double)); 
    sample_permutation(*G, samp_groups, params->seed);
    for (g=0; g<*G; g++)
      {        
      params->g=&samp_groups[g];
      lim[0]=1.0e-8; lim[1]=1.0e1;
      bb(lim, tol); 
      }
    d_vector[5]=diff_max(V_alpha, old_alpha, *G);
    } else d_vector[5] = 0.0;

  R_CheckUserInterrupt();
  if (d_vector[6] > *tol)
    {
    flag=4;
    memcpy(old_nu, V_nu, *G*sizeof(double)); 
    sample_permutation(*G, samp_groups, params->seed);
    for (g=0; g<*G; g++)
      {
      params->g=&samp_groups[g];
      lim[0]=1.0e-8; lim[1]=1.0e1;
      bb(lim, tol);
      }
    d_vector[6]=diff_max(V_nu, old_nu, *G);
    } else d_vector[6] = 0.0;
  R_CheckUserInterrupt();
  if (d_vector[7] > *tol)
    {
    flag=0;
    memcpy(old_xi, V_xi, *P*sizeof(double));
    sample_permutation(JP, samp_coeffs, params->seed);
    for (p=0;p<JP;p++) 
      {
      params->p=&samp_coeffs[p];
      //params->p=&p;
      lim[0]=-2.0e1; lim[1]=2.0e1;
      bb(lim, tol);
      }
    d_vector[7]=diff_max(V_xi, old_xi, *P);
    } else d_vector[7] = 0.0;
  R_CheckUserInterrupt();
  if (d_vector[8] > *tol)
    {
    flag=5;
    memcpy(old_psi2, V_psi2, *P*sizeof(double));
    sample_permutation(JP, samp_coeffs, params->seed);
    for (p=0;p<JP;p++) 
      {
      params->p=&samp_coeffs[p];
      //params->p=&p;
      lim[0]=1.0e-8; lim[1]=1.0e1;
      bb(lim, tol); 
      }
    d_vector[8]=diff_max(V_psi2, old_psi2, *P);
    } else d_vector[8] = 0.0;
  
  if ( (d_vector[0] < *tol) && (d_vector[1] < *tol) && (d_vector[2] < *tol) && (d_vector[3] < *tol) && 
       (d_vector[4] < *tol) && (d_vector[5] < *tol) && (d_vector[6] < *tol) && (d_vector[7] < *tol) && (d_vector[8] < *tol) )
    {
    *params->conv=1;
    printf("*************************************************\n");
    printf("All converged after %d iterations\n", l+1);
    printf("*************************************************\n");
    break;
    }
  /*
  printf("%d / %d iterations: KL = %lf\n", l+1, *steps, tmp);
  */
  printf("%d / %d iterations\n", l+1, *steps);
/*
  for (i=0; i<9; i++)
    {
    if (d_vector[i] > *tol)
      printf("%s change was %e\n", V_names[i], d_vector[i]);
    //else printf("%s converged after %d iterations\n", V_names[i], l+1);
    }
*/
  }
    //free(V_names);
    free(samp_nodes);
    free(samp_groups);
    free(samp_coeffs);
    free(old_xi);
    free(old_psi2);
    free(old_z);
    free(old_sigma2);
    free(old_alpha);
    free(old_nu);
    free(old_eta);
    free(old_omega2);
    free(old_lambda);
    free(params->dists);
    //free(params);
    return;
    }
