#include <stdio.h>
#include <R.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include "gsl/gsl_multimin.h"
#include "headers.h"

double KL_V_xi_n(const gsl_vector *V_xi_n, void *null);
double KL_V_xi_e(const gsl_vector *V_xi_e, void *null);
double KL_V_z_i(const gsl_vector *V_z_i, void *null);
double KL_V_sigma2_i(const gsl_vector *log_V_sigma2_i, void *null);
double KL_V_alpha_g(const gsl_vector *log_V_alpha_g, void *null);
double KL_V_nu_g(const gsl_vector *log_V_nu_g, void *null);
double KL_V_psi2_n(const gsl_vector *log_V_psi2_n, void *null);
double KL_V_psi2_e(const gsl_vector *log_V_psi2_e, void *null);
void gr_KL_V_xi_n(const gsl_vector *V_xi_n, void *null, gsl_vector *df);
void gr_KL_V_xi_e(const gsl_vector *V_xi_e, void *null, gsl_vector *df);
void gr_KL_V_z_i(const gsl_vector *V_z_i, void *null, gsl_vector *df);
void gr_KL_V_sigma2_i(const gsl_vector *V_sigma2_i, void *null, gsl_vector *df);
void gr_KL_V_alpha_g(const gsl_vector *V_alpha_g, void *null, gsl_vector *df);
void gr_KL_V_nu_g(const gsl_vector *V_nu_g, void *null, gsl_vector *df);
void gr_KL_V_psi2_n(const gsl_vector *V_psi2_n, void *null, gsl_vector *df);
void gr_KL_V_psi2_e(const gsl_vector *V_psi2_e, void *null, gsl_vector *df);
void xi_n_fdf(const gsl_vector *V_xi_n, void *null, double *f, gsl_vector *df);
void xi_e_fdf(const gsl_vector *V_xi_e, void *null, double *f, gsl_vector *df);
void z_i_fdf(const gsl_vector *V_z_i, void *null, double *f, gsl_vector *df);
void sigma2_i_fdf(const gsl_vector *log_V_sigma2_i, void *null, double *f, gsl_vector *df);
void alpha_g_fdf(const gsl_vector *log_V_alpha_g, void *null, double *f, gsl_vector *df);
void nu_g_fdf(const gsl_vector *log_V_nu_g, void *null, double *f, gsl_vector *df);
void psi2_n_fdf(const gsl_vector *log_V_psi2_n, void *null, double *f, gsl_vector *df);
void psi2_e_fdf(const gsl_vector *log_V_psi2_e, void *null, double *f, gsl_vector *df);


void optim()
{
  int dd, status;
  size_t iter = 0;
  int *max_iter; // ok to be since they're all unidimensional updates / minimisations
  gsl_multimin_fdfminimizer *s;
  
  gsl_multimin_function_fdf F;
  max_iter=params->MAX_ITER;
  gsl_vector *x;
  int flag = *params->flag, SIZE;
  double *tmp, sum_tmp;
  SIZE=(flag==1 ? *params->D : 1);
  x = gsl_vector_alloc (SIZE);
  tmp=calloc(SIZE, sizeof(double));
  F.n = SIZE;
  switch (flag)
    {
    case 0: gsl_vector_set(x, 0, params->V_xi_n[*params->i+*params->N* *params->p]);
      F.f = &KL_V_xi_n;
      F.df = &gr_KL_V_xi_n;
      F.fdf = &xi_n_fdf;
      break;
    case 6: gsl_vector_set(x, 0, params->V_xi_e[*params->p]);
      F.f = &KL_V_xi_e;
      F.df = &gr_KL_V_xi_e;
      F.fdf = &xi_e_fdf;
      break;
    case 1: 
      for (dd=0; dd<*params->D; dd++)
        gsl_vector_set(x, dd, params->V_z[*params->i* *params->D+dd]); 
      F.f = &KL_V_z_i;
      F.df = &gr_KL_V_z_i;
      F.fdf = &z_i_fdf;
      break;
    case 2: gsl_vector_set(x, 0, (params->V_sigma2[*params->i])); 
      F.f = &KL_V_sigma2_i;
      F.df = &gr_KL_V_sigma2_i;
      F.fdf = &sigma2_i_fdf;
      break;
    case 3: gsl_vector_set(x, 0, (params->V_alpha[*params->g])); 
      F.f = &KL_V_alpha_g;
      F.df = &gr_KL_V_alpha_g;
      F.fdf = &alpha_g_fdf;
      break;
    case 4: gsl_vector_set(x, 0, (params->V_nu[*params->g])); 
      F.f = &KL_V_nu_g;
      F.df = &gr_KL_V_nu_g;
      F.fdf = &nu_g_fdf;
      break;
    case 5: gsl_vector_set(x, 0, (params->V_psi2_n[*params->p])); 
      F.f = &KL_V_psi2_n;
      F.df = &gr_KL_V_psi2_n;
      F.fdf = &psi2_n_fdf;
      break;
    case 7: gsl_vector_set(x, 0, (params->V_psi2_e[*params->p])); 
      F.f = &KL_V_psi2_e;
      F.df = &gr_KL_V_psi2_e;
      F.fdf = &psi2_e_fdf;
      break;
    default: gsl_vector_set(x, 0, 6);
      break;
    }
  
  s = gsl_multimin_fdfminimizer_alloc (gsl_multimin_fdfminimizer_conjugate_fr, SIZE);
  gsl_multimin_fdfminimizer_set (s, &F, x, 1.0e-2, 1.0e-4);
  for (dd=0; dd<SIZE; dd++)
        tmp[dd] = gsl_vector_get(s->x,dd)-1.0;
  do
    {
      iter++;
      for (dd=0; dd<SIZE; dd++)
        tmp[dd] = gsl_vector_get(s->x,dd);
      status = gsl_multimin_fdfminimizer_iterate (s);
      sum_tmp=0.0;
      for (dd=0; dd<SIZE; dd++)
        sum_tmp += fabs(gsl_vector_get(s->x,dd)-tmp[dd]);
      if (sum_tmp<eps)
        break;
      if (status)
             break;
      status = gsl_multimin_test_gradient (s->gradient, 1.0e-4); 
  }
  while (status == GSL_CONTINUE && iter < *max_iter);
  switch (flag)
    {
    case 0: params->V_xi_n[*params->i+*params->N* *params->p] = gsl_vector_get (s->x, 0); 
      break;
    case 6: params->V_xi_e[*params->p] = gsl_vector_get (s->x, 0); 
      break;
    case 1: 
      for (dd=0; dd<*params->D; dd++) 
        {
        params->V_z[*params->i* *params->D+dd] = gsl_vector_get (s->x, dd); 
	}
      //getchar();
      break;
    case 2: params->V_sigma2[*params->i] = (gsl_vector_get (s->x, 0)); 
      break;
    case 3: params->V_alpha[*params->g] = (gsl_vector_get (s->x, 0)); 
      break;
    case 4: params->V_nu[*params->g] = (gsl_vector_get (s->x, 0)); 
      break;
    case 5: params->V_psi2_n[*params->p] = (gsl_vector_get (s->x, 0)); 
      break;
    case 7: params->V_psi2_e[*params->p] = (gsl_vector_get (s->x, 0)); 
      break;
    }
  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);
  free(tmp);
  return;
}
