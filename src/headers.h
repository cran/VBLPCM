#include <stdlib.h>
//#define SQRT 
#define SQRT sqrt
#define GSQRT 
#define CONST 1

/* plain old macros for general use */
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b)) 

#define SUBSET N // for KL_funcs
#define STRATSUB1 MAX((int)(*params->NnonE* *params->STRAT),1)
#define STRATSUB2 MAX((int)(Nnon* *params->STRAT),1)

typedef struct 
{
  int *MAX_ITER;
  int *p;
  int *P;
  int *d;
  int *D;
  int *N; 
  int *NE; // #edges
  int *NnonE; // #non-edges
  int *NM; //#missing edges
  int *g;
  int *G;
  int *i;
  int *numedges;
  int *diam;
  int *hopslist;
  double *Y;
  double *dists;
  int *E; // edges matrix. NE X 2
  int *nonE; // non-edges matrix. NM X 2
  int *M; // missing-edges matrix. NM X 2
  int *EnonE; //edges and non-edges matrix for nodes
  double *XX; // design matrix for covariates. May also be used for sender / receiver effects, etc. 
  double *V_xi;
  double *V_psi2;
  double *V_z;
  double *V_sigma2;
  double *V_eta;
  double *V_lambda;
  double *V_omega2;
  double *V_alpha;
  double *V_nu;
  double *xi;
  double *psi2;
  double *sigma2;
  double *omega2;
  double *alpha;
  double *nu;
  double *inv_sigma02;
  int *flag;
  int *model;
  double *STRAT;
  double *seed;
  int *conv;
  }Rf_params;  

int flag;
Rf_params *params;


double logistic_log_like();

#define eps 1.0e-6
#define loglikefunc logistic_log_like

void sample_permutation(int N, int *samp, double *seed);


