#include "headers.h"
#include <math.h>

double loglikefunc(Rf_params *params)
  {
  double log_like = 0.0, tmp, cov, cov2;
  int i=*params->i, j, d, p, k;
  int N=*params->N;
  int P_n=*params->P_n;
  int diam=*params->diam, h, hsum, Nnon;
  int *sample_non_edges;
  int *sample_nodes;
  int NC1, NC2;
  if (*params->flag==1 || *params->flag==2)
    {
    for (j = 2+diam; j < params->hopslist[i*(CONST+diam+N)+1]+2+diam; j++) // loop over all edges
      {
      tmp = 0.0;
      for (d = 0; d < *params->D; d++)
        tmp += pow (params->V_z[i* *params->D + d] - params->V_z[(params->hopslist[i*(CONST+diam+N)+j]-1)* *params->D + d], 2.0);
      tmp = SQRT (tmp + *params->D*(params->V_sigma2[i] + params->V_sigma2[params->hopslist[i*(CONST+diam+N)+j]-1]));
      //tmp = -1.0/SQRT (tmp + *params->D*(params->V_sigma2[i] + params->V_sigma2[params->hopslist[i*(CONST+diam+N)+j]-1])); // NEW dists
      // i->j part
      cov=0.0; cov2=0.0;
      if (*params->imodel==1)
        cov += params->V_xi_n[i];
      if (*params->imodel==2)
        cov += params->V_xi_n[(params->hopslist[i*(CONST+diam+N)+j]-1)];
      if (*params->imodel==3)
        cov += params->V_xi_n[i] + params->V_xi_n[N+(params->hopslist[i*(CONST+diam+N)+j]-1)];
      for (p=0;p<*params->P_e;p++)
        {
        cov += params->V_xi_e[p]*params->XX_e[(i*N + params->hopslist[i*(CONST+diam+N)+j]-1)* *params->P_e + p];
        cov2+= params->V_psi2_e[p]*params->XX_e[(i*N + params->hopslist[i*(CONST+diam+N)+j]-1)* *params->P_e + p];
        }
      for (p=0;p<P_n;p++)
        cov2+= params->V_psi2_n[p];
      log_like += (cov-tmp)*params->Y[i*N+(params->hopslist[i*(CONST+diam+N)+j]-1)] - log(1.0+exp(cov+0.5*cov2-tmp));
      // j->i part
      cov=0.0; cov2=0.0;
      if (*params->imodel==1)
        cov += params->V_xi_n[(params->hopslist[i*(CONST+diam+N)+j]-1)];
      if (*params->imodel==2)
        cov += params->V_xi_n[i];
      if (*params->imodel==3)
        cov += params->V_xi_n[N+i] + params->V_xi_n[(params->hopslist[i*(CONST+diam+N)+j]-1)];
      for (p=0;p<*params->P_e;p++)
        {
        cov += params->V_xi_e[p]*params->XX_e[((params->hopslist[i*(CONST+diam+N)+j]-1)*N+i)* *params->P_e + p];
        cov2+= params->V_psi2_e[p]*params->XX_e[((params->hopslist[i*(CONST+diam+N)+j]-1)*N+i)* *params->P_e + p];
        }
      for (p=0;p<P_n;p++)
        cov2+= params->V_psi2_n[p];
      log_like += (cov-tmp)*params->Y[(params->hopslist[i*(CONST+diam+N)+j]-1)*N+i] - log(1.0+exp(cov+0.5*cov2-tmp));
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
          for (d = 0; d < *params->D; d++)
            tmp += pow (params->V_z[i* *params->D + d] - params->V_z[j* *params->D + d], 2.0);
          tmp = SQRT (tmp + *params->D*(params->V_sigma2[i] + params->V_sigma2[j]));
          //tmp = -1.0/SQRT (tmp + *params->D*(params->V_sigma2[i] + params->V_sigma2[j])); // NEW dists
          // i->j part
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
          log_like += -Nnon/(NC2)*log(1.0+exp(cov+0.5*cov2-tmp));
          // j->i part
          cov=0.0; cov2=0.0;
          if (*params->imodel==1)
            cov += params->V_xi_n[j];
          if (*params->imodel==2)
            cov += params->V_xi_n[i];
          if (*params->imodel==3)
            cov += params->V_xi_n[N+i] + params->V_xi_n[j];
          for (p=0;p<P_n;p++)
            cov2+= params->V_psi2_n[p];
          for (p=0;p<*params->P_e;p++)
            {
            cov += params->V_xi_e[p]*params->XX_e[(j*N + i)* *params->P_e + p];
            cov2+= params->V_psi2_e[p]*params->XX_e[(j*N + i)* *params->P_e + p];
            }
          log_like += -Nnon/(NC2)*log(1.0+exp(cov+0.5*cov2-tmp));
          }
        hsum += Nnon;
        free(sample_nodes);
        }
      }
    }
  else 
    {
    sample_non_edges = calloc(*params->NnonE, sizeof(int));
    for (i=0;i<*params->NE;i++) // loop over all edges
      {
      tmp = 0.0;
      // note that to covert from R indexing to C indexing i becomes i-1 and i+1 becomes i
      for (d = 0; d < *params->D; d++)
        tmp += pow (params->V_z[(params->E[i*2]-1)* *params->D + d] - params->V_z[(params->E[i*2+1]-1)* *params->D + d], 2.0);
      tmp = SQRT (tmp + *params->D*(params->V_sigma2[params->E[i*2]-1] + params->V_sigma2[params->E[i*2+1]-1]));
      //tmp = -1.0/SQRT (tmp + *params->D*(params->V_sigma2[params->E[i*2]-1] + params->V_sigma2[params->E[i*2+1]-1])); // NEW dists
      cov=0.0; cov2=0.0;
      if (*params->imodel==1)
        cov += params->V_xi_n[(params->E[i*2]-1)];
      if (*params->imodel==2)
        cov += params->V_xi_n[(params->E[i*2+1]-1)];
      if (*params->imodel==3)
        cov += params->V_xi_n[(params->E[i*2]-1)] + params->V_xi_n[N+(params->E[i*2+1]-1)];
      for (p=0;p<P_n;p++)
        cov2+= params->V_psi2_n[p];
      for (p=0;p<*params->P_e;p++)
        {
        cov += params->V_xi_e[p]*params->XX_e[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* *params->P_e + p];
        cov2+= params->V_psi2_e[p]*params->XX_e[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* *params->P_e + p];
        }
      log_like += -tmp - log(1.0+exp(cov+0.5*cov2-tmp));
      }
    sample_permutation(*params->NnonE, sample_non_edges, params->seed);
    NC1=MIN(*params->NnonE, (int)(*params->NC* *params->NE));
    for (j=0;j<NC1;j++) // loop over NC1 of the non-edges
      {
      i=sample_non_edges[j];
      tmp = 0.0;
      for (d = 0; d < *params->D; d++)
        tmp += pow (params->V_z[(params->nonE[i*2]-1)* *params->D + d] - 
	            params->V_z[(params->nonE[i*2+1]-1)* *params->D + d], 2.0);
      tmp = SQRT (tmp + *params->D*(params->V_sigma2[params->nonE[i*2]-1] + params->V_sigma2[params->nonE[i*2+1]-1]));
      //tmp = -1.0/SQRT (tmp + *params->D*(params->V_sigma2[params->nonE[i*2]-1] + params->V_sigma2[params->nonE[i*2+1]-1])); // NEW dists
      cov=0.0; cov2=0.0;
      if (*params->imodel==1)
        cov += params->V_xi_n[(params->nonE[i*2]-1)];
      if (*params->imodel==2)
        cov += params->V_xi_n[(params->nonE[i*2+1]-1)];
      if (*params->imodel==3)
        cov += params->V_xi_n[(params->nonE[i*2]-1)] + params->V_xi_n[N+(params->nonE[i*2+1]-1)];
      for (p=0;p<P_n;p++)
        cov2+= params->V_psi2_n[p];
      for (p=0;p<*params->P_e;p++)
        {
        cov += params->V_xi_e[p]*params->XX_e[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)* *params->P_e + p];
        cov2+= params->V_psi2_e[p]*params->XX_e[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)* *params->P_e + p];
        }
      log_like += -*params->NnonE/NC1*log(1.0+exp(cov+0.5*cov2-tmp));
      }
    free(sample_non_edges);
    }
  return log_like;
  }


