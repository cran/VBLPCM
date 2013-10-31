\name{vblpcmfit}
\alias{vblpcmfit}
\title{fit the variational model through EM type iterations
}
\description{Perform optimisation of the variational parameters of the variational approximation to
the posterior for the latent position cluster model for network data. 
}
\usage{vblpcmfit(variational.start, STEPS = 50, maxiter = 100, tol=1e-6, NC=NULL, 
           seed=NaN, d_vector=rep(TRUE,9))
}
\arguments{
  \item{variational.start}{The starting configuration; use vblpcmstart() to generate this. 
}
  \item{STEPS}{Maximum number of iterations in the main VBEM loop. 
}
  \item{maxiter}{Maximum number of iterations for the internal univariate optimisation loops. 
}
  \item{tol}{tolerance of change in variational parameter updates below which the algorithm is deemed 
  to have converged for that parameter. 
}
  \item{NC}{Number of non-links sampled in the case-control type sampler. 
  Results in a speedup but loss of accuracy. 
}
  \item{seed}{Optional seed for the random number generator. Supplying NaN is equivalent to
  not supplying it. Supply a value so that results may be replicated. 
}
  \item{d_vector}{Optional logical vector specifying which sets of variational parameters are to be
  updated. See Details for more information. 
}
}
\value{A v.params list containing the fitted variational parameters for the 
latent positions, clustering membership probabilities, etc. 
conv indicated whether convergence was obtained within the specified number of iterations. 
}
\references{
Michael Salter-Townshend and Thomas Brendan Murphy (2009).
"Variational Bayesian Inference for the Latent Position Cluster Model."
Workshop on Analyzing Networks and Learning with Graphs.
Neural Information Processing Systems.
}
\details{d_vector is a logical vector of length 9 that can be used to select which variational
parameters are held fixed and which are updated. The parameters are in the following order:
z (latent positions), sigma2 (variance of latent positions), lambda (membership probability 
matrix), eta (cluster centres), omega2 (cluster variances), alpha (cluster specific variance of
nodes), nu (Dirichlet parameter for marginal cluster probabilities), xi (likelihood intercept term 
mean), psi2 (likelihood intercept term variance). 
}
\author{Michael Salter-Townshend
}
\seealso{vblpcmstart, latentnet::ergmm
}
