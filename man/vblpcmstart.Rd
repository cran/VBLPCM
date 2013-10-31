\name{vblpcmstart}
\alias{vblpcmstart}
\title{Generate sensible starting configuration for the variational parameter set. 
}
\description{Uses fast methods to generate sensible and coherent values for the parameters of the
variational method. There are returned as a list and that list may be passed directly to
vblpcmfit(). User specification of the configuration is recommended as tweaks to this list only. 
}
\usage{vblpcmstart(g.network,G=1,d=2,LSTEPS=5e3,model="plain", CLUST=0, B=NULL,
           lcc=TRUE, edgecovs=NULL,sendcovs=NULL,receivecovs=NULL,
           socialcovs=NULL,START="FR", seed=NaN)
}
\arguments{
  \item{g.network}{a network object
}
  \item{G}{Desired number of groups
}
  \item{d}{Desired dimensionality of the latent space
}
  \item{LSTEPS}{Number of steps in the log-likelihood forces algorithm
}
  \item{model}{model specified as "plain", "rreceiver", "rsender" or "rsocial". See vblpcmcovs for details. 
}
  \item{CLUST}{degree of push to clustering at the start
}
  \item{B}{default intercept value
}
  \item{lcc}{logical indicator. TRUE => analyze largest connected component 
  of g.network only FALSE => analyze the whole network. 
}
\item{edgecovs}{optional edge covariates.
}
\item{sendcovs}{optional sender node covariates.
}
\item{receivecovs}{optional receiver node covariates.
}
\item{socialcovs}{optional sociality node covariates.
}
\item{START}{what to start the initial positions with. 
                "FR" for Fruchterman-Reingold. 
                "geodesic" for geodesic distances.
		"laplace" for using the Graph Laplacian.
                "random" for random. 
}
\item{seed}{Optional seed for the random number generator in R. Equivalent to using set.seed(seed). 
  The default NaN value does not call set.seed(). 
}
}
\value{A v.params list containing the latent positions, clustering membership probabilities, etc. 
}
\author{Michael Salter-Townshend
}
\seealso{vblpcmfit, vblpcmcovs
}
\examples{
data(sampson)
### plot the mean posterior positions with initial estimations for variational parameters
plot(vblpcmstart(samplike,G=3),main="Sampson's Monks: VB Initial Values")
### plot the mean posterior positions with final estimations for variational parameters
plot(vblpcmfit(vblpcmstart(samplike,G=3)),main="Sampson's Monks: VB Solution")
}
