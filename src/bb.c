#include <R.h>
#include <gsl/gsl_vector.h>
#include "headers.h"

void gr_KL_V_xi_n(const gsl_vector *V_xi_n, void *null, gsl_vector *df);
void gr_KL_V_xi_e(const gsl_vector *V_xi_e, void *null, gsl_vector *df);
void gr_KL_V_z_i(const gsl_vector *V_z_i, void *null, gsl_vector *df);
void gr_KL_V_sigma2_i(const gsl_vector *V_sigma2_i, void *null, gsl_vector *df);
void gr_KL_V_alpha_g(const gsl_vector *V_alpha_g, void *null, gsl_vector *df);
void gr_KL_V_nu_g(const gsl_vector *V_nu_g, void *null, gsl_vector *df);
void gr_KL_V_psi2_n(const gsl_vector *V_psi2_n, void *null, gsl_vector *df);
void gr_KL_V_psi2_e(const gsl_vector *V_psi2_e, void *null, gsl_vector *df);

int signum(double x)
{
  if (x<0) {
    return -1;
  }
  else if (x==0) {
    return 0;
  }
  else if (x>0) {
    return 1;
  }
  else {
    /* must be that x is Not-a-Number */
    return 2;
  }
}

void F (const gsl_vector *V, void *null, gsl_vector *df)
  {
  switch (*params->flag)
    {
    case 0: 
      gr_KL_V_xi_n(V,NULL,df);
      break;
    case 6: 
      gr_KL_V_xi_e(V,NULL,df);
      break;
    case 1: 
      gr_KL_V_z_i(V,NULL,df);
      break;
    case 2: 
      gr_KL_V_sigma2_i(V,NULL,df);
      break;
    case 3: 
      gr_KL_V_alpha_g(V,NULL,df);
      break;
    case 4: 
      gr_KL_V_nu_g(V,NULL,df);
      break;
    case 5: 
      gr_KL_V_psi2_n(V,NULL,df);
      break;
    case 7: 
      gr_KL_V_psi2_e(V,NULL,df);
      break;
    default: 
      break;
    }
  return;
  }

void bb(double *lim, double *tol)
{
  size_t iter = 0; 
  double val[2]={0.0,0.0};
  gsl_vector *a, *b, *ab, *tmp;
  a = gsl_vector_alloc (1);
  b = gsl_vector_alloc (1);
  ab = gsl_vector_alloc (1);
  tmp = gsl_vector_alloc (1);
  
    gsl_vector_set(a,0,lim[0]);
    gsl_vector_set(b,0,lim[1]);
    do
      {
      iter++;
      gsl_vector_set(ab, 0, 0.5*(gsl_vector_get(a,0)+gsl_vector_get(b,0)));
      F(a,NULL,tmp);
      val[0] = gsl_vector_get(tmp,0);
      F(ab,NULL,tmp);
      val[1] = gsl_vector_get(tmp,0);
      if (signum(val[0])==signum(val[1])) gsl_vector_set(a,0,gsl_vector_get(ab,0)); 
      else gsl_vector_set(b,0,gsl_vector_get(ab,0));
      gsl_vector_set(ab, 0, 0.5*(gsl_vector_get(a,0)+gsl_vector_get(b,0)));
      if (fabs(gsl_vector_get(a,0)-gsl_vector_get(b,0)) < *tol)
        break;
      }
      while (iter < *params->MAX_ITER);
  // assign the new value
  switch (*params->flag)
    {
    case 0: params->V_xi_n[*params->i+ *params->N* *params->p] = gsl_vector_get (ab, 0); 
      break;
    case 6: params->V_xi_e[*params->p] = gsl_vector_get (ab, 0); 
      break;
    case 1: params->V_z[*params->i* *params->D + *params->d] = gsl_vector_get (ab, 0); 
      break;
    case 2: params->V_sigma2[*params->i] = gsl_vector_get (ab, 0); 
      break;
    case 3: params->V_alpha[*params->g] = gsl_vector_get (ab, 0); 
      break;
    case 4: params->V_nu[*params->g] = gsl_vector_get (ab, 0); 
      break;
    case 5: params->V_psi2_n[*params->p] = gsl_vector_get (ab, 0); 
      break;
    case 7: params->V_psi2_e[*params->p] = gsl_vector_get (ab, 0); 
      break;
    }
  gsl_vector_free (a);
  gsl_vector_free (b);
  gsl_vector_free (ab);
  gsl_vector_free (tmp);
  return;
}
