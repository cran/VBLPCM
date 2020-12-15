#include <math.h>
#include <gsl/gsl_sf.h>
#include "headers.h"

double
KL_V_xi_e (const gsl_vector *v_V_xi_e, Rf_params *params)
{
  params->V_xi_e[*params->p] = gsl_vector_get(v_V_xi_e, 0);
  double tmpsum = 0.0, tmp;
  double KL;
  tmpsum = loglikefunc(params);
  tmp = -pow (params->V_xi_e[*params->p] - *params->xi, 2.0);
  KL = fabs (tmpsum + 0.5*(tmp) / *params->psi2);
  //KL = fabs (tmpsum + (tmp) / *params->psi2);
  return KL;
}

double
KL_V_xi_n (const gsl_vector *v_V_xi_n, Rf_params *params)
{
  int P_n=*params->P_n;
  params->V_xi_n[*params->i* P_n+*params->p] = gsl_vector_get(v_V_xi_n, 0);
  double tmpsum = 0.0, tmp;
  double KL;
  tmpsum = loglikefunc(params);
  tmp = -pow (params->V_xi_n[*params->i + *params->N* *params->p] - *params->xi, 2.0);
  KL = fabs (tmpsum + 0.5*(tmp) / *params->psi2);
  //KL = fabs (tmpsum + (tmp) / *params->psi2);
  return KL;
}

void gr_KL_V_xi_e (const gsl_vector *v_V_xi_e, void *null, gsl_vector *df, Rf_params *params)
{
  int i, j, p=*params->p, pn;
  int P_e=*params->P_e, P_n=*params->P_n;
  params->V_xi_e[p] = gsl_vector_get(v_V_xi_e, 0);
  double tmpsum = 0.0, cov, cov2;
  int N = *params->N;
  double KL;
  int *sample_non_edges = calloc(*params->NnonE, sizeof(int));
  int NC1;
  for (i = 0; i < *params->NE; i++) // loop over all edges
    {
    cov = params->V_xi_e[p]*params->XX_e[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* P_e + p];
    cov2= params->V_psi2_e[p]*params->XX_e[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* P_e + p]; 
    if (*params->imodel==1)
      cov += params->V_xi_n[(params->E[i*2]-1)];
    if (*params->imodel==2)
      cov += params->V_xi_n[(params->E[i*2+1]-1)];
    if (*params->imodel==3)
      cov += params->V_xi_n[(params->E[i*2]-1)] + params->V_xi_n[N+(params->E[i*2+1]-1)];
    for (pn=0;pn<P_n;pn++)
      cov2+= params->V_psi2_n[pn];
    tmpsum += params->XX_e[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* P_e + p]*
              (1.0 - 1.0/(1.0 + exp (-cov + params->dists[((params->E[i*2]-1)*N + params->E[i*2+1]-1)] - 0.5 * cov2)));
    }
  sample_permutation(*params->NnonE, sample_non_edges, params->seed);
  NC1=MIN(*params->NnonE, (int)(*params->NC* *params->NE));
  for (j=0;j<NC1;j++) // loop over a sample of the non-edges
    {
    i=sample_non_edges[j];
    cov = params->V_xi_e[p]*params->XX_e[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)* P_e + p];
    cov2= params->V_psi2_e[p]*params->XX_e[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)* P_e + p]; 
    if (*params->imodel==1)
      cov += params->V_xi_n[(params->nonE[i*2]-1)];
    if (*params->imodel==2)
      cov += params->V_xi_n[(params->nonE[i*2+1]-1)];
    if (*params->imodel==3)
      cov += params->V_xi_n[(params->nonE[i*2]-1)] + params->V_xi_n[N+(params->nonE[i*2+1]-1)];
    for (pn=0;pn<P_n;pn++)
      cov2+= params->V_psi2_n[pn];
    tmpsum += *params->NnonE/NC1*(params->XX_e[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)* P_e + p]*
                (- 1.0/(1.0 + exp (-cov + params->dists[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)] - 0.5 * cov2))));
      }
  KL = tmpsum - 0.5*(params->V_xi_e[p] - *params->xi) / *params->psi2;
  free(sample_non_edges);
  gsl_vector_set(df, 0, -KL);
  return;
}
void gr_KL_V_xi_n (const gsl_vector *v_V_xi_n, void *null, gsl_vector *df, Rf_params *params)
{
  int i=*params->i, p=*params->p, P_e=*params->P_e, j, pe;
  int N = *params->N;
  params->V_xi_n[i+N*p] = gsl_vector_get(v_V_xi_n, 0);
  double tmpsum = 0.0, cov=0.0, cov2;
  int *sample_nodes, hsum, h, Nnon, k;
  double KL;
  int diam=*params->diam;
  int NC2;
  for (j = 2+diam; j < params->hopslist[i*(CONST+diam+N)+1]+2+diam; j++) // loop over all edges
    {
    // i->j part
    cov2 = 0.0;
    if (*params->imodel==1)
      cov = params->V_xi_n[i];
    if (*params->imodel==2)
       cov = params->V_xi_n[(params->hopslist[i*(CONST+diam+N)+j]-1)];
    if (*params->imodel==3)
      cov = params->V_xi_n[i] + params->V_xi_n[N+(params->hopslist[i*(CONST+diam+N)+j]-1)];
    for (pe=0;pe<P_e;pe++)
      {
      cov += params->V_xi_e[pe]*params->XX_e[(i*N + params->hopslist[i*(CONST+diam+N)+j]-1)* P_e + pe];  
      cov2+= params->V_psi2_e[pe]*params->XX_e[(i*N + params->hopslist[i*(CONST+diam+N)+j]-1)* P_e + pe];
      }
    cov2 += params->V_psi2_n[0];
    if (*params->imodel==3)
      cov2 += params->V_psi2_n[1];
    tmpsum += (params->Y[i*N+(params->hopslist[i*(CONST+diam+N)+j]-1)] - 
               1.0/(1.0 + exp (-cov + params->dists[i*N + (params->hopslist[i*(CONST+diam+N)+j]-1)] - 0.5 * cov2)));
    // j->i part
    cov2 = 0.0;
    if (*params->imodel==1)
      cov = params->V_xi_n[(params->hopslist[i*(CONST+diam+N)+j]-1)];
    if (*params->imodel==2)
      cov = params->V_xi_n[i];
    if (*params->imodel==3)
      cov = params->V_xi_n[(params->hopslist[i*(CONST+diam+N)+j]-1)] + params->V_xi_n[N+i];
    for (pe=0;pe<P_e;pe++)
      {
      cov += params->V_xi_e[pe]*params->XX_e[((params->hopslist[i*(CONST+diam+N)+j]-1)*N+i)* P_e + pe];  
      cov2+= params->V_psi2_e[pe]*params->XX_e[((params->hopslist[i*(CONST+diam+N)+j]-1)*N+i)* P_e + pe];
      }
    cov2 += params->V_psi2_n[0];
    if (*params->imodel==3)
      cov2 += params->V_psi2_n[1];
    tmpsum += (params->Y[(params->hopslist[i*(CONST+diam+N)+j]-1)*N+i] - 
               1.0/(1.0 + exp (-cov + params->dists[(params->hopslist[i*(CONST+diam+N)+j]-1)*N+i] - 0.5 * cov2)));
    }
  hsum = 0;
  for (h=2;h<1+diam;h++) // <1+diam because we don't sample unknown edges
    {
    Nnon=params->hopslist[i*(CONST+diam+N)+h];
    if (Nnon>0)
      {
      sample_nodes = calloc(Nnon, sizeof(int));
      sample_permutation(Nnon, sample_nodes, params->seed);
      NC2=MIN(Nnon, (int)(*params->NC*params->hopslist[i*(CONST+diam+N)+1]));
      for (k=0;k<NC2;k++)  // loop over some of the non-edges
        {
        j=params->hopslist[i*(CONST+diam+N)+2+diam+params->hopslist[i*(CONST+diam+N)+1]+hsum+sample_nodes[k]]-1;
	// this part for Yij = 0 
        if (*params->imodel==1)
          cov = params->V_xi_n[i]; 
        if (*params->imodel==2)
          cov = params->V_xi_n[j];
        if (*params->imodel==3)
          cov = params->V_xi_n[i] + params->V_xi_n[N+j];
    cov2 = params->V_psi2_n[0];
    if (*params->imodel==3)
      cov2 += params->V_psi2_n[1];
	for (pe=0;pe<P_e;pe++)
          {
          cov += params->V_xi_e[pe]*params->XX_e[(i*N + j)* P_e + pe];
          cov2+= params->V_psi2_e[pe]*params->XX_e[(i*N + j)* P_e + pe];
          }
        tmpsum += Nnon/NC2*(- 1.0/(1.0 + exp (-cov + params->dists[i*N + j] - 0.5 * cov2)));
	// this part for Yji = 0 
	if (*params->imodel==1)
          cov = params->V_xi_n[j];
        if (*params->imodel==2)
          cov = params->V_xi_n[i];
        if (*params->imodel==3)
          cov = params->V_xi_n[N+i] + params->V_xi_n[j];
    cov2 = params->V_psi2_n[0];
    if (*params->imodel==3)
      cov2 += params->V_psi2_n[1];
        for (pe=0;pe<P_e;pe++)
          {
          cov += params->V_xi_e[pe]*params->XX_e[(j*N + i)* P_e + pe];
          cov2+= params->V_psi2_e[pe]*params->XX_e[(j*N + i)* P_e + pe];
          }
        tmpsum += Nnon/NC2*(- 1.0/(1.0 + exp (-cov + params->dists[i*N + j] - 0.5 * cov2)));
        }
      hsum += Nnon;
      free(sample_nodes);
      }
    }
  KL = tmpsum - 0.5*(params->V_xi_n[i+N*p] - *params->xi) / *params->psi2;
  gsl_vector_set(df, 0, -KL);
  return;
}
 void
     xi_e_fdf (const gsl_vector *x, 
             void *null, double *f, gsl_vector *df, Rf_params *params)
     {
       *f = KL_V_xi_e(x, params);
       gr_KL_V_xi_e(x, NULL, df, params);
       return;
     }
void
     xi_n_fdf (const gsl_vector *x, 
             void *null, double *f, gsl_vector *df, Rf_params *params)
     {
       *f = KL_V_xi_n(x, params);
       gr_KL_V_xi_n(x, NULL, df, params);
       return;
     }

double KL_V_z_i (const gsl_vector *v_V_z_i, Rf_params *params)
{
  int i = *params->i, d = *params->d;
  int D = *params->D;	
  int g;
  double tmpsum = 0.0, tmp;
  int N = *params->N;
  double KL;
  for (d=0; d<D; d++)
    params->V_z[i*D+d] = gsl_vector_get(v_V_z_i, d);
  tmpsum = loglikefunc(params);
  KL = tmpsum;
  tmpsum = 0;
  for (g = 0; g < *params->G; g++)
    {
      tmp = 0.0;
        for (d=0; d<D; d++)
          tmp += pow(params->V_z[i*D+d] - params->V_eta[g * D + d], 2.0);
      tmp = GSQRT(tmp + params->V_sigma2[i] + params->V_omega2[g]);
      tmpsum +=
	params->V_lambda[g * N + i] * (D * gsl_sf_psi (0.5 * *params->inv_sigma02 *
				params->V_alpha[g]) -0.5 * *params->inv_sigma02 *
				     params->V_alpha[g] * tmp);
    }
  KL = fabs (KL + tmpsum);
  return KL;
}

void gr_KL_V_z_i (const gsl_vector *v_V_z_i, void *null, gsl_vector *df, Rf_params *params)
{
  int i = *params->i, j, k, d;
  int D = *params->D;     
  double tmp;
  for (d=0; d<D; d++)
    params->V_z[i * D + d] = gsl_vector_get(v_V_z_i, d);
  int g, p, G = *params->G;
  int N = *params->N, h, Nnon, hsum;
  int diam=*params->diam;
  int P_n=*params->P_n;
  double *KL = calloc(D, sizeof(double));
  double *tmpsum = calloc(D, sizeof(double));
  double cov, cov2;
  int *sample_nodes;
  int NC2;
  for (j = 2+diam; j < params->hopslist[i*(CONST+diam+N)+1]+2+diam; j++) // loop over all edges
    {
    tmp = 0.0;
    for (d=0; d<D; d++)
      tmp += pow(params->V_z[i * D + d] - params->V_z[(params->hopslist[i*(CONST+diam+N)+j]-1) * D + d], 2.0);
    tmp = SQRT(tmp + D *(params->V_sigma2[i] + params->V_sigma2[(params->hopslist[i*(CONST+diam+N)+j]-1)]));
    cov=0.0; cov2=0.0;
    if (params->Y[i*N+(params->hopslist[i*(CONST+diam+N)+j]-1)]>0) // i is sender and j is receiver
      {
      if (*params->imodel==1)
        cov += params->V_xi_n[i];
      if (*params->imodel==2)
        cov += params->V_xi_n[(params->hopslist[i*(CONST+diam+N)+j]-1)];
      if (*params->imodel==3)
        cov += params->V_xi_n[i] + params->V_xi_n[N+(params->hopslist[i*(CONST+diam+N)+j]-1)];
      for (p=0;p<*params->P_e;p++)
        {
        cov += params->V_xi_e[p]*params->XX_e[(i*N + (params->hopslist[i*(CONST+diam+N)+j]-1))* *params->P_e + p];
        cov2+= params->V_psi2_e[p]*params->XX_e[(i*N + (params->hopslist[i*(CONST+diam+N)+j]-1))* *params->P_e + p];
        }
      }
    if (params->Y[(params->hopslist[i*(CONST+diam+N)+j]-1)*N+i]>0) // j is sender and i is receiver
      {
      if (*params->imodel==1)
        cov += params->V_xi_n[(params->hopslist[i*(CONST+diam+N)+j]-1)];
      if (*params->imodel==2)
        cov += params->V_xi_n[i]; 
      if (*params->imodel==3)
        cov += params->V_xi_n[(params->hopslist[i*(CONST+diam+N)+j]-1)] + params->V_xi_n[N+i];
      for (p=0;p<*params->P_e;p++)
        {
        cov += params->V_xi_e[p]*params->XX_e[((params->hopslist[i*(CONST+diam+N)+j]-1)*N+i)* *params->P_e + p];
        cov2+= params->V_psi2_e[p]*params->XX_e[((params->hopslist[i*(CONST+diam+N)+j]-1)*N+i)* *params->P_e + p];
        }
      }
    for (p=0;p<P_n;p++)
      cov2+= params->V_psi2_n[p];
    for (d=0; d<D; d++)
        tmpsum[d] +=  (params->V_z[i * D + d] - params->V_z[(params->hopslist[i*(CONST+diam+N)+j]-1) * D + d]) *
                      (1.0 - 1.0/(1.0+exp(-cov + tmp - 0.5*cov2)));
    }
  hsum = 0;
  for (h=2;h<1+diam;h++) // <1+diam because we don't sample unknown edges
    {
    Nnon=params->hopslist[i*(CONST+diam+N)+h];
    if (Nnon>0)
      {
      sample_nodes = calloc(Nnon, sizeof(int));
      sample_permutation(Nnon, sample_nodes, params->seed);
      NC2=MIN(Nnon, (int)(*params->NC*params->hopslist[i*(CONST+diam+N)+1]));
      for (k=0;k<NC2;k++)  // loop over some of the non-edges
        {
        j=params->hopslist[i*(CONST+diam+N)+2+diam+params->hopslist[i*(CONST+diam+N)+1]+hsum+sample_nodes[k]]-1;
        tmp = 0.0;
        for (d=0; d<D; d++)
          tmp += pow(params->V_z[i * D + d] - params->V_z[j * D + d], 2.0);
        tmp = SQRT(tmp + D *(params->V_sigma2[i] + params->V_sigma2[j]));
        cov=0.0; cov2=0.0;
        if (*params->imodel==1)
          cov += params->V_xi_n[i];
        if (*params->imodel==2)
          cov += params->V_xi_n[j];
        if (*params->imodel==3)
          cov += params->V_xi_n[i] + params->V_xi_n[N+j];
        for (p=0;p<P_n;p++)
          cov2+= params->V_psi2_n[p];
        for (p=0;p<*params->P_e;p++)
          {
          cov += params->V_xi_e[p]*params->XX_e[(i*N + j)* *params->P_e + p];
          cov2+= params->V_psi2_e[p]*params->XX_e[(i*N + j)* *params->P_e + p];
          }
        for (d=0; d<D; d++)
            tmpsum[d] +=  Nnon/NC2*(
                          (params->V_z[i * D + d] - params->V_z[j * D + d]) *
                          (- 1.0/(1.0+exp(-cov + tmp - 0.5*cov2))));
        }
      hsum += Nnon;
      free(sample_nodes);
      }
    }
  for (d=0; d<D; d++)
     KL[d] = tmpsum[d];
  for (d=0; d<D; d++)
    tmpsum[d] = 0.0; // there's probably a function for this in C
  for (g = 0; g < G; g++)
    for (d=0; d<D; d++)
      tmpsum[d] = tmpsum[d] + params->V_lambda[g * N + i] *
      0.5* *params->inv_sigma02 * params->V_alpha[g] *
      fabs(params->V_eta[g * D + d] - params->V_z[i * D + d]);
  for (d=0; d<D; d++)
    KL[d] -= tmpsum[d];
    //KL[d] += tmpsum[d];
  for (d=0; d<D; d++)
    gsl_vector_set(df, d, -KL[d]);
  free(tmpsum);
  free(KL);
  return;
}


 void
     z_i_fdf (const gsl_vector *x, 
             void *null, double *f, gsl_vector *df, Rf_params *params)
     {
       *f = KL_V_z_i(x, params);
       gr_KL_V_z_i(x, NULL, df, params);
       return;
     }

double
KL_V_sigma2_i (const gsl_vector *v_V_sigma2_i, Rf_params *params)
{
  int i = *params->i, j;
  int D = *params->D;	
  int g, G = *params->G;
  double tmpsum = 0.0;
  double KL;
  int N = *params->N;
  params->V_sigma2[*params->i] = gsl_vector_get(v_V_sigma2_i, 0);
  tmpsum = loglikefunc(params);
  double tmpsum1 = 0.0;
  for (g = 0; g < G; g++)
    tmpsum1 =
      tmpsum1 -
      params->V_lambda[g * N + i] * 0.5 * *params->inv_sigma02 * params->V_alpha[g] * params->V_sigma2[i];
  
  double tmpsum2 = 0.0;
  for (j = 0; j < SUBSET; j++) if (j != i)
      tmpsum2 += log (params->V_sigma2[j]);
  tmpsum2 += log (params->V_sigma2[*params->i]);
  
  KL = fabs (tmpsum + tmpsum1  + 0.5 * D * tmpsum2);
  return KL;
}

void gr_KL_V_sigma2_i (const gsl_vector *v_V_sigma2_i, void *null, gsl_vector *df, Rf_params *params)
{
  int i = *params->i, j, k, d = *params->d;
  int D = *params->D;     
  int g, G = *params->G, hsum, h;
  double tmpsum = 0.0, tmp;
  int p, N = *params->N;
  int diam=*params->diam;
  int P_n=*params->P_n;
  double KL;
  double cov, cov2;
  double V_sigma2_i = gsl_vector_get(v_V_sigma2_i, 0);
  int Nnon=N-params->hopslist[i*(CONST+diam+N)+1]-params->hopslist[i*(CONST+diam+N)+1+diam];
  int *sample_nodes = 0;
  int NC2;
  for (j = 2+diam; j < params->hopslist[i*(CONST+diam+N)+1]+2+diam; j++) // loop over all edges
    {
    tmp = 0.0;
    for (d = 0; d < D; d++)
      tmp += pow (params->V_z[i * D + d] - params->V_z[(params->hopslist[i*(CONST+diam+N)+j]-1) * D + d], 2.0);
    cov=0.0; cov2=0.0;
    if (params->Y[i*N+(params->hopslist[i*(CONST+diam+N)+j]-1)]>0) // i is sender and j is receiver
      {
      if (*params->imodel==1)
        cov += params->V_xi_n[i];
      if (*params->imodel==2)
        cov += params->V_xi_n[(params->hopslist[i*(CONST+diam+N)+j]-1)];
      if (*params->imodel==3)
        cov += params->V_xi_n[i] + params->V_xi_n[N+(params->hopslist[i*(CONST+diam+N)+j]-1)];
      for (p=0;p<*params->P_e;p++)
        {
        cov += params->V_xi_e[p]*params->XX_e[(i*N + (params->hopslist[i*(CONST+diam+N)+j]-1))* *params->P_e + p];
        cov2+= params->V_psi2_e[p]*params->XX_e[(i*N + (params->hopslist[i*(CONST+diam+N)+j]-1))* *params->P_e + p];
        }
      }
    if (params->Y[(params->hopslist[i*(CONST+diam+N)+j]-1)*N+i]>0) // j is sender and i is receiver
      {
      if (*params->imodel==1)
        cov += params->V_xi_n[(params->hopslist[i*(CONST+diam+N)+j]-1)];
      if (*params->imodel==2)
        cov += params->V_xi_n[i];
      if (*params->imodel==3)
        cov += params->V_xi_n[(params->hopslist[i*(CONST+diam+N)+j]-1)] + params->V_xi_n[N+i];
      for (p=0;p<*params->P_e;p++)
        {
        cov += params->V_xi_e[p]*params->XX_e[((params->hopslist[i*(CONST+diam+N)+j]-1)*N+i)* *params->P_e + p];
        cov2+= params->V_psi2_e[p]*params->XX_e[((params->hopslist[i*(CONST+diam+N)+j]-1)*N+i)* *params->P_e + p];
        }
      }
    for (p=0;p<P_n;p++)
      cov2+= params->V_psi2_n[p];
    tmp = SQRT (tmp + D * (V_sigma2_i + params->V_sigma2[(params->hopslist[i*(CONST+diam+N)+j]-1)]));
    tmpsum = tmpsum - D + D / (V_sigma2_i * (1.0 + exp (-cov + tmp - 0.5 * cov2)));
    }
  hsum = 0;
  for (h=2;h<1+diam;h++) // <1+diam because we don't sample unknown edges
    {
    Nnon=params->hopslist[i*(CONST+diam+N)+h];
    if (Nnon>0)
      {
      sample_nodes = calloc(Nnon, sizeof(int));
      sample_permutation(Nnon, sample_nodes, params->seed);
      NC2=MIN(Nnon, (int)(*params->NC*params->hopslist[i*(CONST+diam+N)+1]));
      for (k=0;k<NC2;k++)  // loop over some of the non-edges
        {
        j=params->hopslist[i*(CONST+diam+N)+2+diam+params->hopslist[i*(CONST+diam+N)+1]+hsum+sample_nodes[k]]-1;
        tmp = 0.0;
        for (d = 0; d < D; d++)
          tmp += pow (params->V_z[i * D + d] - params->V_z[j * D + d], 2.0);
        cov=0.0; cov2=0.0;
        if (*params->imodel==1)
          cov += params->V_xi_n[i];
        if (*params->imodel==2)
          cov +=  params->V_xi_n[j];
        if (*params->imodel==3)
          cov += params->V_xi_n[i] + params->V_xi_n[N+j];
        for (p=0;p<P_n;p++)
          cov2+= params->V_psi2_n[p];
        for (p=0;p<*params->P_e;p++)
          {
          cov += params->V_xi_e[p]*params->XX_e[(i*N + j)* *params->P_e + p];
          cov2+= params->V_psi2_e[p]*params->XX_e[(i*N + j)* *params->P_e + p];
          }
          tmp = SQRT (tmp + D * (V_sigma2_i + params->V_sigma2[j]));
          tmpsum = tmpsum + Nnon/NC2*(D / (V_sigma2_i * (1.0 + exp (-cov + tmp - 0.5 * cov2))));
        }
      hsum += Nnon;
      free(sample_nodes);
      }
    }
  tmp = 0.0;
  for (g = 0; g < G; g++)
    tmp += params->V_lambda[g * N + i] * *params->inv_sigma02 * params->V_alpha[g];
  KL = tmpsum - 0.5 * tmp + 0.5 * D / V_sigma2_i;
  gsl_vector_set(df, 0, -KL);
  return;
}


 void
     sigma2_i_fdf (const gsl_vector *x, 
             void *null, double *f, gsl_vector *df, Rf_params *params)
     {
       *f = KL_V_sigma2_i(x, params);
       gr_KL_V_sigma2_i(x, NULL, df, params);
       return;
     }


double
KL_V_alpha_g (const gsl_vector *v_V_alpha_g, Rf_params *params)
{
  int i = *params->i;
  int d;
  int G = *params->G, g = *params->g;
  double tmpsum = 0.0, tmp;
  int N = *params->N;
  double KL;
  double V_alpha_g = gsl_vector_get(v_V_alpha_g, 0);
  for (i = 0; i < SUBSET; i++)
    {
  tmp = 0.0;
  for (d=0; d<*params->D; d++)
    tmp += pow(params->V_z[i * *params->D + d] - params->V_eta[g * *params->D + d], 2.0);
      tmpsum +=
	params->V_lambda[*params->g * N + i] * 
	              (*params->D*gsl_sf_psi(*params->inv_sigma02*V_alpha_g) -
                       0.5* *params->inv_sigma02*V_alpha_g*(params->V_sigma2[i]+params->V_omega2[g]+tmp));
    }
  tmp = 0.0;
  for (g=0; g<G; g++)
    if (g!=*params->g)
      tmp += params->V_alpha[g];
  KL =
    fabs (tmpsum*(tmp+V_alpha_g) + lgamma(0.5*V_alpha_g) +
   0.5*(params->alpha[*params->g]-V_alpha_g)*(gsl_sf_psi(0.5*V_alpha_g)));
   //-gsl_sf_psi(params->alpha[*params->g])));
  return KL;
}

void gr_KL_V_alpha_g (const gsl_vector *v_V_alpha_g, void *null, gsl_vector *df, Rf_params *params)
{
  int i = *params->i;
  int g = *params->g;
  double tmpsum = 0.0, tmp;
  int N = *params->N;
  int d;
  double KL;
  double V_alpha_g = gsl_vector_get(v_V_alpha_g, 0);
  for (i = 0; i < SUBSET; i++)
    {
  tmp = 0.0;
  for (d=0; d<*params->D; d++)
    tmp += pow(params->V_z[i * *params->D + d] - params->V_eta[g * *params->D + d], 2.0);
      tmpsum +=
	params->V_lambda[*params->g * N + i] * 
	(*params->D* *params->inv_sigma02* gsl_sf_psi_1(0.5* *params->inv_sigma02*V_alpha_g)-
                       0.5* *params->inv_sigma02*(params->V_sigma2[i]+params->V_omega2[g]+tmp));
    }
  KL = tmpsum + (0.5*(V_alpha_g-params->alpha[*params->g])* gsl_sf_psi_1(0.5*V_alpha_g));// - 
		     //gsl_sf_psi(0.5*V_alpha_g)) + gsl_sf_psi(0.5*V_alpha_g)/lgamma(0.5*V_alpha_g);
  gsl_vector_set(df, 0, -KL);
  return;
}

 void
     alpha_g_fdf (const gsl_vector *x, 
             void *null, double *f, gsl_vector *df, Rf_params *params)
     {
       *f = KL_V_alpha_g(x, params);
       gr_KL_V_alpha_g(x, NULL, df, params);
       return;
     }


double
KL_V_nu_g (const gsl_vector *v_V_nu_g, Rf_params *params)
{
  int i = *params->i;	
  int g, G = *params->G;
  double tmpsum = 0.0, tmp;
  int N = *params->N;
  double KL;
  double V_nu_g = gsl_vector_get(v_V_nu_g, 0);
  tmp = 0.0;
  for (g = 0; g < G; g++)
    if (g!=*params->g)
      tmp += params->V_nu[g];
  for (i = 0; i < SUBSET; i++)
    tmpsum +=
      params->V_lambda[*params->g * N +
			     i] * (gsl_sf_psi (V_nu_g) - gsl_sf_psi (tmp+V_nu_g));
  KL = tmpsum - lgamma (tmp+V_nu_g);
  tmp = 0.0;
  for (g = 0; g < G; g++)
    if (g!=*params->g)
      tmp += lgamma (params->V_nu[g]);
  KL = fabs (KL + tmp + lgamma (V_nu_g) -
	     (V_nu_g -
	      params->nu[*params->g]) *
	     (gsl_sf_psi (V_nu_g) -
	      gsl_sf_psi (params->nu[*params->g])));
  return KL;
}

void gr_KL_V_nu_g (const gsl_vector *v_V_nu_g, void *null, gsl_vector *df, Rf_params *params)
{
  int i = *params->i;	//, j, d = *params->d;
  int g, G = *params->G;
  double tmpsum = 0.0, tmp;
  int N = *params->N;
  double KL;
  double V_nu_g = gsl_vector_get(v_V_nu_g, 0);
  tmp = 0.0;
  for (g = 0; g < G; g++)
    if (g!=*params->g)
      tmp += params->V_nu[g];
  for (i = 0; i < SUBSET; i++)
    tmpsum =
      tmpsum + params->V_lambda[*params->g * N + i] * (gsl_sf_psi_1 (V_nu_g) -
						    gsl_sf_psi_1 (tmp+V_nu_g));
  KL =
    tmpsum - gsl_sf_psi (tmp+V_nu_g) - gsl_sf_psi (params->nu[*params->g]) -
    V_nu_g * gsl_sf_psi_1 (V_nu_g) +
    params->nu[*params->g] * gsl_sf_psi_1 (V_nu_g);
  gsl_vector_set(df, 0, -KL);
  return;
}

 void
     nu_g_fdf (const gsl_vector *x, 
             void *null, double *f, gsl_vector *df, Rf_params *params)
     {
       *f = KL_V_nu_g(x, params);
       gr_KL_V_nu_g(x, NULL, df, params);
       return;
     }


double
KL_V_psi2_e (const gsl_vector *v_V_psi2_e, Rf_params *params)
{
  int D = *params->D;
  params->V_psi2_e[*params->p] = gsl_vector_get(v_V_psi2_e, 0);
  double tmpsum = 0.0;
  tmpsum = loglikefunc(params);
  double KL = tmpsum +
	    0.5*( D * (log (params->V_psi2_e[*params->p])-log(*params->psi2)) -
	    D * params->V_psi2_e[*params->p] / *params->psi2);
  KL = fabs (KL);
  return KL;
}
double
KL_V_psi2_n (const gsl_vector *v_V_psi2_n, Rf_params *params)
{
  int D = *params->D;
  params->V_psi2_n[*params->p] = gsl_vector_get(v_V_psi2_n, 0);
  double tmpsum = 0.0;
  tmpsum = loglikefunc(params);
  double KL = tmpsum +
	    0.5*( D * (log (params->V_psi2_n[*params->p])-log(*params->psi2)) -
	    D * params->V_psi2_n[*params->p] / *params->psi2);
  KL = fabs (KL);
  return KL;
}

void gr_KL_V_psi2_e (const gsl_vector *v_V_psi2_e, void * null, gsl_vector *df, Rf_params *params)
{
  int j, i = *params->i, p=*params->p, pn;
  int D = *params->D, P_e = *params->P_e, P_n = *params->P_n;
  double tmp1=0.0, cov, cov2;
  int N = *params->N;
  double KL;
  int *sample_non_edges = calloc(*params->NnonE, sizeof(int));
  int NC1;
  params->V_psi2_e[*params->p] = gsl_vector_get(v_V_psi2_e, 0);
  KL = 0.0;
  for (i = 0; i < *params->NE; i++) // loop over all edges
      {
      cov = params->V_xi_e[p]*params->XX_e[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* P_e + p];
      cov2= params->V_psi2_e[p]*params->XX_e[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* P_e + p];
      if (*params->imodel==1)
        cov += params->V_xi_n[(params->E[i*2]-1)];
      if (*params->imodel==2)
        cov += params->V_xi_n[(params->E[i*2+1]-1)];
      if (*params->imodel==3)
        cov += params->V_xi_n[(params->E[i*2]-1)] + params->V_xi_n[N+(params->E[i*2+1]-1)];
      for (pn=0;pn<P_n;pn++)
        cov2+= params->V_psi2_n[pn];
      tmp1 += -0.5*params->XX_e[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* P_e + *params->p]/
               (1.0+exp(-cov + params->dists[((params->E[i*2]-1)*N + params->E[i*2+1]-1)] - 0.5*cov2));
      }
  sample_permutation(*params->NnonE, sample_non_edges, params->seed);
  NC1=MIN(*params->NnonE, (int)(*params->NC* *params->NE));
  for (j=0;j<NC1;j++) // loop over a sample of the non-edges
      {
      i=sample_non_edges[j];
      cov = params->V_xi_e[p]*params->XX_e[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)* P_e + p];
      cov2= params->V_psi2_e[p]*params->XX_e[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)* P_e + p];
      if (*params->imodel==1)
        cov += params->V_xi_n[(params->nonE[i*2]-1)];
      if (*params->imodel==2)
        cov += params->V_xi_n[(params->nonE[i*2+1]-1)];
      if (*params->imodel==3)
        cov += params->V_xi_n[(params->nonE[i*2]-1)] + params->V_xi_n[N+(params->nonE[i*2+1]-1)];
      for (pn=0;pn<P_n;pn++)
        cov2+= params->V_psi2_n[pn];
      tmp1 += *params->NnonE/NC1*(-0.5*params->XX_e[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)* P_e + *params->p]/
               (1.0+exp(-cov + params->dists[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)] - 0.5*cov2)));
      }
  KL = tmp1 + 0.5 * (D / params->V_psi2_e[*params->p] - D / *params->psi2);
  gsl_vector_set(df, 0, -KL);
  free(sample_non_edges);
  return;
}

void gr_KL_V_psi2_n (const gsl_vector *v_V_psi2_n, void * null, gsl_vector *df, Rf_params *params)
{
  int i, p=*params->p, j, pe;
  int D = *params->D;     
  double tmp1=0.0, cov=0.0, cov2;
  int N = *params->N;
  int P_e = *params->P_e;
  double KL;
  int *sample_non_edges = calloc(*params->NnonE, sizeof(int));
  int NC1;
  params->V_psi2_n[*params->p] = gsl_vector_get(v_V_psi2_n, 0);
  KL = 0.0;
  for (i = 0; i < *params->NE; i++) // loop over all edges
    {
    if (*params->imodel==1)
      cov = params->V_xi_n[(params->E[i*2]-1)];
    if (*params->imodel==2)
      cov = params->V_xi_n[(params->E[i*2+1]-1)];
    if (*params->imodel==3)
      cov = params->V_xi_n[(params->E[i*2]-1)] + params->V_xi_n[N+(params->E[i*2+1]-1)];
    cov2 = params->V_psi2_n[p];
    for (pe=0;pe<P_e;pe++)
      {
      cov += params->V_xi_e[pe]*params->XX_e[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* P_e + pe];
      cov2+= params->V_psi2_e[pe]*params->XX_e[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* P_e + pe];
      }
    tmp1 += -0.5/(1.0+exp(-cov + params->dists[((params->E[i*2]-1)*N + params->E[i*2+1]-1)] - 0.5*cov2));
    }
  sample_permutation(*params->NnonE, sample_non_edges, params->seed);
  NC1=MIN(*params->NnonE, (int)(*params->NC* *params->NE));
  for (j=0;j<NC1;j++) // loop over a sample of the non-edges
    {
    i=sample_non_edges[j];
    if (*params->imodel==1)
      cov = params->V_xi_n[(params->nonE[i*2]-1)];
    if (*params->imodel==2)
      cov = params->V_xi_n[(params->nonE[i*2+1]-1)];
    if (*params->imodel==3)
      cov = params->V_xi_n[(params->nonE[i*2]-1)] + params->V_xi_n[N+(params->nonE[i*2+1]-1)];
    cov2= params->V_psi2_n[p];
    for (pe=0;pe<P_e;pe++)
      {
      cov += params->V_xi_e[pe]*params->XX_e[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)* P_e + pe];
      cov2+= params->V_psi2_e[pe]*params->XX_e[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)* P_e + pe];
      }
    tmp1 += *params->NnonE/NC1*(-0.5/(1.0+exp(-cov + params->dists[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)] - 0.5*cov2)));
    }
  KL = tmp1 + 0.5 * (D / params->V_psi2_n[*params->p] - D / *params->psi2);
  gsl_vector_set(df, 0, -KL);
  free(sample_non_edges);
  return;
}

 void
     psi2_e_fdf (const gsl_vector *x, 
             void *null, double *f, gsl_vector *df, Rf_params *params)
     {
       *f = KL_V_psi2_e(x, params);
       gr_KL_V_psi2_e(x, NULL, df, params);
       return;
     }
 void
     psi2_n_fdf (const gsl_vector *x, 
             void *null, double *f, gsl_vector *df, Rf_params *params)
     {
       *f = KL_V_psi2_n(x, params);
       gr_KL_V_psi2_n(x, NULL, df, params);
       return;
     }
void KL_total (int *imodel,
  int *P_n,
  int *P_e,
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
  double *XX_e,
  double *V_xi_n,
  double *V_xi_e,
  double *V_psi2_n,
  double *V_psi2_e,
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
  double *seed,
  int *NC,
  double *KL)
{ 
  int p, i, g, d, flag=0;
  double tmp;
  Rf_params *tmp_params;
  tmp_params=calloc(1,sizeof(Rf_params));
  tmp_params->seed=seed;
  tmp_params->p=&p;
  tmp_params->i=&i;
  tmp_params->g=&g;
  tmp_params->d=&d;
  tmp_params->flag=&flag;
  tmp_params->P_n=P_n;
  tmp_params->P_e=P_e;
  tmp_params->D=D;
  tmp_params->N=N;
  tmp_params->NE=NE;
  tmp_params->NM=NM;
  tmp_params->NnonE=NnonE;
  tmp_params->G=G;
  tmp_params->Y=Y;
  tmp_params->E=E; // edges matrix. NE X 2
  tmp_params->nonE=nonE; // non-edges matrix. NM X 2
  tmp_params->M=M; // missing-edges matrix. NM X 2
  tmp_params->numedges=numedges;
  tmp_params->EnonE=EnonE;
  tmp_params->diam=diam; 
  tmp_params->hopslist=hopslist; 
  tmp_params->XX_e=XX_e; // design matrix for edge covariates. 
  tmp_params->V_xi_n=V_xi_n;
  tmp_params->V_xi_e=V_xi_e;
  tmp_params->V_psi2_n=V_psi2_n;
  tmp_params->V_psi2_e=V_psi2_e;
  tmp_params->V_z=V_z;
  tmp_params->V_sigma2=V_sigma2;
  tmp_params->V_eta=V_eta;
  tmp_params->V_lambda=V_lambda;
  tmp_params->V_omega2=V_omega2;
  tmp_params->V_nu=V_nu;
  tmp_params->V_alpha=V_alpha;
  tmp_params->xi=xi;
  tmp_params->psi2=psi2;
  tmp_params->sigma2=sigma2;
  tmp_params->omega2=omega2;
  tmp_params->nu=nu;
  tmp_params->alpha=alpha;
  tmp_params->inv_sigma02=inv_sigma02;
  tmp_params->NC=NC;
  tmp_params->imodel=imodel;
  flag=0;
  // p1
  *KL = loglikefunc(tmp_params); 
  // p2
  for (g=0;g<*G;g++)
    for (i=0;i<*N;i++)
      {
      tmp = 0.0;
      for (d=0; d<*D; d++)
        tmp += pow(V_z[i* *D+d] - V_eta[g * *D + d], 2.0);
      tmp = GSQRT(tmp + V_sigma2[i] + V_omega2[g]);
      *KL += V_lambda[g* *N+i]*(*D * gsl_sf_psi (0.5 * *inv_sigma02 *V_alpha[g])- 0.5 * *inv_sigma02 * tmp_params->V_alpha[g] * tmp);
      }
  // p3
  tmp=0;
  for (g=0;g<*G;g++)
    tmp+=V_nu[g];
  for (g=0;g<*G;g++)
    for (i=0;i<*N;i++)
      *KL += V_lambda[g* *N+i]*(gsl_sf_psi(V_nu[g])-gsl_sf_psi(tmp));
  // -q2
  tmp=0; 
  for (i=0;i<*N;i++)
    tmp+=log(V_sigma2[i]);
  *KL += 0.5* *N *(1.0-log(2.0*M_PI)+0.5* *D *tmp);
  // -q3
  for (g=0;g<*G;g++)
    for (i=0;i<*N;i++)
      *KL += V_lambda[g* *N+i]*log(V_lambda[g* *N+i]);
  // p4-q4
  tmp=0;
  for (g=0;g<*G;g++)
    tmp+=V_nu[g];
  *KL += -lgamma(tmp);
  tmp=0;
  for (g=0;g<*G;g++)
    tmp+=nu[g];
  *KL += lgamma(tmp);
  for (g=0;g<*G;g++)
    *KL += lgamma(V_nu[g]) - lgamma(nu[g]) - (V_nu[g]-nu[g])*(gsl_sf_psi(V_nu[g])-gsl_sf_psi(nu[g]));
  // p5-q5
  if (*P_n > 0)
    for (p=0;p<*P_n;p++)
      {
      *KL += 0.5*(*D*(log(V_psi2_n[p])-log(*psi2)) - *D*(V_psi2_n[p]/ *psi2));
      for (i=0;i<*N;i++)
        *KL += 0.5*(- pow(V_xi_n[i+ *N*p]-xi[p],2.0)/ *psi2 +*D);
      }
  for (p=0;p<*P_e;p++)
    *KL += 0.5*(*D*(log(V_psi2_e[p])-log(*psi2)) - *D*(V_psi2_e[p]/ *psi2) - pow(V_xi_e[p]-xi[p],2.0)/ *psi2 +*D);
  // p6-q6
  for (g=0;g<*G;g++)
    *KL += 0.5*(*D*(log(V_omega2[g])-log(*omega2)) - *D*(V_omega2[g]/ *omega2) - pow(V_eta[g]-0.0,2.0)/ *omega2 +*D);
  // p7-q7
  for (g=0;g<*G;g++)
    *KL += lgamma(0.5*V_alpha[g])-lgamma(0.5*alpha[g]) + 0.5*(alpha[g]-V_alpha[g])*(gsl_sf_psi(0.5*V_alpha[g]));
  //*KL = -*KL;
  return;
}


