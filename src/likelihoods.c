#include "headers.h"
#include <math.h>

double logistic_log_like()
  {
  double log_like = 0.0, tmp, cov, cov2;
  int i=*params->i, j, d, p, k;
  int N=*params->N;
  int diam=*params->diam, h, hsum, Nnon;
  int *sample_non_edges;
  int *sample_nodes;
   if (*params->flag==1 || *params->flag==2)
    {
    for (j = 2+diam; j < params->hopslist[i*(CONST+diam+N)+1]+2+diam; j++) // loop over all edges
        {
        tmp = 0.0;
        for (d = 0; d < *params->D; d++)
          tmp += pow (params->V_z[i* *params->D + d] - params->V_z[(params->hopslist[i*(CONST+diam+N)+j]-1)* *params->D + d], 2.0);
        tmp = SQRT (tmp + *params->D*(params->V_sigma2[i] + params->V_sigma2[params->hopslist[i*(CONST+diam+N)+j]-1]));
        cov=0.0; cov2=0.0;
              for (p=0;p<*params->P;p++)
                {
                cov += params->V_xi[p]*params->XX[(i*N + params->hopslist[i*(CONST+diam+N)+j]-1)* *params->P + p];
                cov2+= params->V_psi2[p]*params->XX[(i*N + params->hopslist[i*(CONST+diam+N)+j]-1)* *params->P + p];
                }
        log_like += (cov-tmp) - log(1.0+exp(cov+0.5*cov2-tmp));
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
          for (d = 0; d < *params->D; d++)
            tmp += pow (params->V_z[i* *params->D + d] - params->V_z[j* *params->D + d], 2.0);
          tmp = SQRT (tmp + *params->D*(params->V_sigma2[i] + params->V_sigma2[j]));
          cov=0.0; cov2=0.0;
          for (p=0;p<*params->P;p++)
            {
            cov += params->V_xi[p]*params->XX[(i*N + j)* *params->P + p];
            cov2+= params->V_psi2[p]*params->XX[(i*N + j)* *params->P + p];
            }
          log_like += -Nnon/(STRATSUB2)*log(1.0+exp(cov+0.5*cov2-tmp));
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
      cov=0.0; cov2=0.0;
      for (p=0;p<*params->P;p++)
        {
        cov += params->V_xi[p]*params->XX[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* *params->P + p];
        cov2+= params->V_psi2[p]*params->XX[((params->E[i*2]-1)*N + params->E[i*2+1]-1)* *params->P + p];
        }
      log_like += -tmp - log(1.0+exp(cov+0.5*cov2-tmp));
      }
    sample_permutation(*params->NnonE, sample_non_edges, params->seed);
    for (j=0;j<STRATSUB1;j++) // loop over STRATSUB1 of the non-edges
      {
      i=sample_non_edges[j];
      tmp = 0.0;
      for (d = 0; d < *params->D; d++)
        tmp += pow (params->V_z[(params->nonE[i*2]-1)* *params->D + d] - params->V_z[(params->nonE[i*2+1]-1)* *params->D + d], 2.0);
      tmp = SQRT (tmp + *params->D*(params->V_sigma2[params->nonE[i*2]-1] + params->V_sigma2[params->nonE[i*2+1]-1]));
      cov=0.0; cov2=0.0;
      for (p=0;p<*params->P;p++)
        {
        cov += params->V_xi[p]*params->XX[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)* *params->P + p];
        cov2+= params->V_psi2[p]*params->XX[((params->nonE[i*2]-1)*N + params->nonE[i*2+1]-1)* *params->P + p];
        }
      log_like += -*params->NnonE/STRATSUB1*log(1.0+exp(cov+0.5*cov2-tmp));
      }
    free(sample_non_edges);
    }
  return log_like;
  }


