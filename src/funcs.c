#include <sys/time.h>
#include "math.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_permutation.h"
#include "headers.h"

void sample_permutation(int N, int *samp, double *seed)
  {
  int i;
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  gsl_rng * r = gsl_rng_alloc (T);
  if (isnan(*seed))
    {
    struct timeval tv;
    gettimeofday(&tv,0);
    gsl_rng_set(r, (tv.tv_sec + tv.tv_usec));
    }
  else 
    gsl_rng_set(r, *seed);
  gsl_permutation * tmp = gsl_permutation_alloc (N);
  gsl_permutation_init (tmp);
  gsl_ran_shuffle (r, tmp->data, N, sizeof(size_t));
  for (i=0; i<N; i++)
    samp[i] = tmp->data[i];
    //samp[i] = i;
  gsl_permutation_free (tmp);
  *seed=gsl_rng_get(r);
  gsl_rng_free(r);
  return;
  }

