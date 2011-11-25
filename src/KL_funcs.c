#include <math.h>
#include <gsl/gsl_sf.h>
#include "headers.h"

double
KL_V_xi (const gsl_vector *v_V_xi)
{
  params->V_xi[*params->p] = gsl_vector_get(v_V_xi, 0);
  double tmpsum = 0.0, tmp;
  double KL;
  tmpsum = loglikefunc();
  tmp = -pow (params->V_xi[*params->p] - *params->xi, 2.0);
  KL = fabs (tmpsum + 0.5*(tmp) / *params->psi2);
  //KL = fabs (tmpsum + (tmp) / *params->psi2);
  return KL;
}

void gr_KL_V_xi (const gsl_vector *v_V_xi, void *null, gsl_vector *df)
{
  params->V_xi[*params->p] = gsl_vector_get(v_V_xi, 0);
  int i, j, p;
  double tmpsum = 0.0, cov, cov2;
  int N = *params->N;
  double KL;
  int *sample_non_edges = calloc(*params->NnonE, sizeof(int));
  for (i = 0; i < *params->NE; i++) // loop over all edges
      {
      cov=0.0; cov2=0.0;
      for (p=0;p<*params->P;p++)
        {
        cov += params->V_xi[p]*params->XX[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* 
	                            *params->P + p];
        cov2+= params->V_psi2[p]*params->XX[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* 
	                              *params->P + p];
        }
      tmpsum += params->XX[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* 
                                   *params->P + *params->p]*
                (1.0 - 1.0/(1.0 + exp (-cov + 
		params->dists[((params->E[i*2]-1)*N + params->E[i*2+1]-1)] - 0.5 * cov2)));
      }
  sample_permutation(*params->NnonE, sample_non_edges, params->seed);
  for (j=0;j<STRATSUB1;j++) // loop over a sample of the non-edges
      {
      i=sample_non_edges[j];
      cov=0.0; cov2=0.0;
      for (p=0;p<*params->P;p++)
        {
        cov += params->V_xi[p]*params->XX[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)*
                                    *params->P + p];
        cov2+= params->V_psi2[p]*params->XX[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)*
                                      *params->P + p];
        }
      tmpsum += *params->NnonE/STRATSUB1*(params->XX[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)*
                                   *params->P + *params->p]*
                (- 1.0/(1.0 + exp (-cov +
                params->dists[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)] - 0.5 * cov2))));
      }
  KL = tmpsum - 0.5*fabs(params->V_xi[*params->p] - *params->xi) / *params->psi2;
  free(sample_non_edges);
  gsl_vector_set(df, 0, -KL);
  return;
}
 void
     xi_fdf (const gsl_vector *x, 
             void *null, double *f, gsl_vector *df)
     {
       *f = KL_V_xi(x);
       gr_KL_V_xi(x, NULL, df);
       return;
     }

double KL_V_z_i (const gsl_vector *v_V_z_i)
{
  int i = *params->i, d = *params->d;
  int D = *params->D;	//, P = *params->P;
  int g;
  double tmpsum = 0.0, tmp;
  int N = *params->N;
  double KL;
  for (d=0; d<D; d++)
    params->V_z[i*D+d] = gsl_vector_get(v_V_z_i, d);
  tmpsum = loglikefunc();
  /*
  if (isnan(tmpsum))
    {
    printf("%d=(%lf,%lf)\n", i, params->V_z[i*D], params->V_z[i*D+1]);
    getchar();
    }
  */
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
  //printf("%d=(%lf,%lf)\n", i, params->V_z[i*D], params->V_z[*i*D+1]);
  return KL;
}

void gr_KL_V_z_i (const gsl_vector *v_V_z_i, void *null, gsl_vector *df)
{
  int i = *params->i, j, k, d;//= *params->d;
  int D = *params->D;     //, P = *params->p;
  double tmp;
  for (d=0; d<D; d++)
    params->V_z[i * D + d] = gsl_vector_get(v_V_z_i, d);
  int g, p, G = *params->G;
  int N = *params->N, h, Nnon, hsum;
  int diam=*params->diam;
  double *KL = calloc(D, sizeof(double));
  double *tmpsum = calloc(D, sizeof(double));
  double cov, cov2;
  int *sample_nodes;
  for (j = 2+diam; j < params->hopslist[i*(CONST+diam+N)+1]+2+diam; j++) // loop over all edges
      {
      tmp = 0.0;
        for (d=0; d<D; d++)
          tmp += pow(params->V_z[i * D + d] - params->V_z[(params->hopslist[i*(CONST+diam+N)+j]-1) * D + d], 2.0);
      tmp = SQRT(tmp + D *(params->V_sigma2[i] + params->V_sigma2[(params->hopslist[i*(CONST+diam+N)+j]-1)]));
      cov=0.0; cov2=0.0;
      for (p=0;p<*params->P;p++)
        {
        cov += params->V_xi[p]*params->XX[(i*N + (params->hopslist[i*(CONST+diam+N)+j]-1))* *params->P + p];
        cov2+= params->V_psi2[p]*params->XX[(i*N + (params->hopslist[i*(CONST+diam+N)+j]-1))* *params->P + p];
        }
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
      for (k=0;k<STRATSUB2;k++)  // loop over some of the non-edges
        {
        j=params->hopslist[i*(CONST+diam+N)+2+diam+params->hopslist[i*(CONST+diam+N)+1]+hsum+sample_nodes[k]]-1;
        tmp = 0.0;
        for (d=0; d<D; d++)
          tmp += pow(params->V_z[i * D + d] - params->V_z[j * D + d], 2.0);
        tmp = SQRT(tmp + D *(params->V_sigma2[i] + params->V_sigma2[j]));
        cov=0.0; cov2=0.0;
        for (p=0;p<*params->P;p++)
          {
          cov += params->V_xi[p]*params->XX[(i*N + j)* *params->P + p];
          cov2+= params->V_psi2[p]*params->XX[(i*N + j)* *params->P + p];
          }
        for (d=0; d<D; d++)
            tmpsum[d] +=  Nnon/STRATSUB2*(
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
             void *null, double *f, gsl_vector *df)
     {
       *f = KL_V_z_i(x);
       gr_KL_V_z_i(x, NULL, df);
       return;
     }

double
KL_V_sigma2_i (const gsl_vector *v_V_sigma2_i)
{
  int i = *params->i, j;
  int D = *params->D;	//, P = *params->P;
  int g, G = *params->G;
  double tmpsum = 0.0;
  double KL;
  int N = *params->N;
  params->V_sigma2[*params->i] = gsl_vector_get(v_V_sigma2_i, 0);
  tmpsum = loglikefunc();
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

void gr_KL_V_sigma2_i (const gsl_vector *v_V_sigma2_i, void *null, gsl_vector *df)
{
  int i = *params->i, j, k, d = *params->d;
  int D = *params->D;     //, P = *params->p;
  int g, G = *params->G, hsum, h;
  double tmpsum = 0.0, tmp;
  int p, N = *params->N;
  int diam=*params->diam;
  double KL;
  double cov, cov2;
  double V_sigma2_i = gsl_vector_get(v_V_sigma2_i, 0);
  int Nnon=N-params->hopslist[i*(CONST+diam+N)+1]-params->hopslist[i*(CONST+diam+N)+1+diam];
  int *sample_nodes = 0;
  for (j = 2+diam; j < params->hopslist[i*(CONST+diam+N)+1]+2+diam; j++) // loop over all edges
      {
        tmp = 0.0;
        for (d = 0; d < D; d++)
          tmp += pow (params->V_z[i * D + d] - params->V_z[(params->hopslist[i*(CONST+diam+N)+j]-1) * D + d], 2.0);
      cov=0.0; cov2=0.0;
      for (p=0;p<*params->P;p++)
        {
        cov += params->V_xi[p]*params->XX[(i*N + (params->hopslist[i*(CONST+diam+N)+j]-1))* *params->P + p];
        cov2+= params->V_psi2[p]*params->XX[(i*N + (params->hopslist[i*(CONST+diam+N)+j]-1))* *params->P + p];
        }
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
      for (k=0;k<STRATSUB2;k++)  // loop over some of the non-edges
        {
        j=params->hopslist[i*(CONST+diam+N)+2+diam+params->hopslist[i*(CONST+diam+N)+1]+hsum+sample_nodes[k]]-1;
        tmp = 0.0;
        for (d = 0; d < D; d++)
          tmp += pow (params->V_z[i * D + d] - params->V_z[j * D + d], 2.0);
        cov=0.0; cov2=0.0;
        for (p=0;p<*params->P;p++)
          {
          cov += params->V_xi[p]*params->XX[(i*N + j)* *params->P + p];
          cov2+= params->V_psi2[p]*params->XX[(i*N + j)* *params->P + p];
          }
          tmp = SQRT (tmp + D * (V_sigma2_i + params->V_sigma2[j]));
          tmpsum = tmpsum + Nnon/STRATSUB2*(D / (V_sigma2_i * (1.0 + exp (-cov + tmp - 0.5 * cov2))));
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
             void *null, double *f, gsl_vector *df)
     {
       *f = KL_V_sigma2_i(x);
       gr_KL_V_sigma2_i(x, NULL, df);
       return;
     }


double
KL_V_alpha_g (const gsl_vector *v_V_alpha_g)
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
    fabs (tmpsum - log(2.0)*(tmp+V_alpha_g) + lgamma(V_alpha_g) +
   (params->alpha[*params->g]-V_alpha_g)*(gsl_sf_psi(V_alpha_g)-log(2.0)));
   //-gsl_sf_psi(params->alpha[*params->g])));
  return KL;
}

void gr_KL_V_alpha_g (const gsl_vector *v_V_alpha_g, void *null, gsl_vector *df)
{
  int i = *params->i;
  int g = *params->g, G=*params->G;
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
  KL = tmpsum + ((V_alpha_g-params->alpha[*params->g])* gsl_sf_psi_1(V_alpha_g) - 
		     gsl_sf_psi(V_alpha_g)) + gsl_sf_psi(V_alpha_g)/lgamma(V_alpha_g) + G*log(2.0);
  gsl_vector_set(df, 0, -KL);
  return;
}

 void
     alpha_g_fdf (const gsl_vector *x, 
             void *null, double *f, gsl_vector *df)
     {
       *f = KL_V_alpha_g(x);
       gr_KL_V_alpha_g(x, NULL, df);
       return;
     }


double
KL_V_nu_g (const gsl_vector *v_V_nu_g)
{
  int i = *params->i;	//, j, d = *params->d;
  //int D = *params->D, P = *params->P;
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

void gr_KL_V_nu_g (const gsl_vector *v_V_nu_g, void *null, gsl_vector *df)
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
             void *null, double *f, gsl_vector *df)
     {
       *f = KL_V_nu_g(x);
       gr_KL_V_nu_g(x, NULL, df);
       return;
     }


double
KL_V_psi2 (const gsl_vector *v_V_psi2)
{
  int D = *params->D;
  params->V_psi2[*params->p] = gsl_vector_get(v_V_psi2, 0);
  double tmpsum = 0.0;
  tmpsum = loglikefunc();
  double KL = tmpsum +
	    0.5*( D * (log (params->V_psi2[*params->p])-log(*params->psi2)) -
	    D * params->V_psi2[*params->p] / *params->psi2);
  KL = fabs (KL);
  return KL;
}

void gr_KL_V_psi2 (const gsl_vector *v_V_psi2, void * null, gsl_vector *df)
{
  int j, i = *params->i, p;
  int D = *params->D;     //, P = *params->p;
  //int g, G = *params->G;
  //double tmpsum = 0.0;
  double tmp1=0.0, cov, cov2;
  int N = *params->N;
  double KL;
  int *sample_non_edges = calloc(*params->NnonE, sizeof(int));
  params->V_psi2[*params->p] = gsl_vector_get(v_V_psi2, 0);
  KL = 0.0;
  for (i = 0; i < *params->NE; i++) // loop over all edges
      {
      cov=0.0; cov2=0.0;
      for (p=0;p<*params->P;p++)
        {
        cov += params->V_xi[p]*params->XX[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* *params->P + p];
        cov2+= params->V_psi2[p]*params->XX[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* *params->P + p];
        }
      tmp1 += -0.5*params->XX[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* *params->P + *params->p]/
               (1.0+exp(-cov + params->dists[((params->E[i*2]-1)*N + params->E[i*2+1]-1)] - 0.5*cov2));
      }
  sample_permutation(*params->NnonE, sample_non_edges, params->seed);
  for (j=0;j<STRATSUB1;j++) // loop over a sample of the non-edges
      {
      i=sample_non_edges[j];
       cov=0.0; cov2=0.0;
      for (p=0;p<*params->P;p++)
        {
        cov += params->V_xi[p]*params->XX[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)* *params->P + p];
        cov2+= params->V_psi2[p]*params->XX[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)* *params->P + p];
        }
      tmp1 += *params->NnonE/STRATSUB1*(-0.5*params->XX[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)* *params->P + *params->p]/
               (1.0+exp(-cov + params->dists[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)] - 0.5*cov2)));
      }
  KL = tmp1 + 0.5 * (D / params->V_psi2[*params->p] - D / *params->psi2);
  gsl_vector_set(df, 0, -KL);
  free(sample_non_edges);
  return;
}



 void
     psi2_fdf (const gsl_vector *x, 
             void *null, double *f, gsl_vector *df)
     {
       *f = KL_V_psi2(x);
       gr_KL_V_psi2(x, NULL, df);
       return;
     }

void KL_total (int *P,
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
  double *dists,
  double *STRAT,
  double *KL)
{ 
  int p, i, g, d, flag;
  double tmp;
  params->p = &p;
  params->i = &i;
  params->g = &g;
  params->d = &d;
  params->flag = &flag;
  params->P=P;
  params->D=D;
  params->N=N;
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
  params->dists=dists;
  params->STRAT=STRAT;
  // p1
  *KL = loglikefunc(); 
  // p2
  for (g=0;g<*G;g++)
    for (i=0;i<*N;i++)
      {
      tmp = 0.0;
      for (d=0; d<*D; d++)
        tmp += pow(V_z[i* *D+d] - V_eta[g * *D + d], 2.0);
      tmp = GSQRT(tmp + V_sigma2[i] + V_omega2[g]);
      *KL += V_lambda[g* *N+i]*(*D * gsl_sf_psi (0.5 * *inv_sigma02 *V_alpha[g])- 0.5 * *inv_sigma02 * params->V_alpha[g] * tmp);
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
  for (p=0;p<*P;p++)
    *KL += 0.5*(*D*(log(V_psi2[p])-log(*psi2)) - *D*(V_psi2[p]/ *psi2) - pow(V_xi[p]-xi[p],2.0)/ *psi2 +*D);
  // p6-q6
  for (g=0;g<*G;g++)
    *KL += 0.5*(*D*(log(V_omega2[g])-log(*omega2)) - *D*(V_omega2[g]/ *omega2) - pow(V_eta[g]-0.0,2.0)/ *omega2 +*D);
  // p7-q7
  for (g=0;g<*G;g++)
    *KL += (V_alpha[g]-alpha[g])*log(2.0) + lgamma(V_alpha[g])-lgamma(alpha[g]) + (alpha[g]-V_alpha[g])*(gsl_sf_psi(V_alpha[g])-log(2.0));
  *KL = -*KL;
  return;
}


