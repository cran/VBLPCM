#include <stdlib.h>
//#define SQRT 
#define SQRT sqrt
#define GSQRT 
#define CONST 1

/* plain old macros for general use */
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b)) 

#define SUBSET N // for KL_funcs

typedef struct 
{
  int *MAX_ITER;
  int *p;
  int *P_n;
  int *P_e;
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
  double *XX_n; // design matrix for node covariates. May also be used for sender / receiver effects, etc. 
  double *XX_e; // design matrix for edge covariates. 
  double *V_xi_n;
  double *V_xi_e;
  double *V_psi2_n;
  double *V_psi2_e;
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
  int *imodel;
  int *NC;
  double *seed;
  int *conv;
  }Rf_params;

#define eps 1.0e-6

void bb(double *lim, double *tol, Rf_params *params); // in bb.c
void sample_permutation(int N, int *samp, double *seed); // in funcs.c
double loglikefunc(Rf_params *params); // in likelihoods.c
void Y_to_E (int *N, int *directed, double *Y, int *E);
void Y_to_nonE (int *N, int *directed, double *Y, int *nonE);
void Y_to_M (int *N, int *directed, double *Y, int *M);
void E_to_Y (int *N, int *NE, int *directed, int *E, double *Y);
void fruchterman_reingold(int *directed, int *N, int *D, int *steps, double *Y, double *X, double *repulserad, double *m, double *volume);
void log_like_forces(int *directed, int *N, int *D, int *steps, double *Y, double *X, double *B, double *m, Rf_params *params);
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
  double *KL);

void Rf_VB_bbs(int *imodel,
  int *steps,
  int *max_iter,
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
  double *tol,
  int *NC,
  double *seed,
  double *d_vector,
  int *conv);
