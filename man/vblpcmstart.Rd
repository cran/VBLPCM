\name{vblpcmstart}
\alias{vblpcmstart}
\title{Generate sensible starting configuration for the variational parameter set. 
}
\description{Uses fast methods to generate sensible and coherent values for the parameters of the
variational method. There are returned as a list and that list may be passed directly to
vblpcmfit(). User specification of the configuration is recommended as tweaks to this list only. 
}
\usage{vblpcmstart(g.network,G=1,d=2,LSTEPS=5e3,model="plain", CLUST=0, 
           lcc=TRUE, edgecovs=NULL,nodecovs=NULL,START="FR", seed=0)
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
  \item{model}{model specified as "plain", "receiver", "sender" or "social"
}
  \item{CLUST}{degree of push to clustering at the start
}
  \item{lcc}{logical indicator. TRUE => analyze largest connected component 
  of g.network only FALSE => analyze the whole network. 
}
\item{edgecovs}{optional edge covariates.
}
\item{nodecovs}{optional node covariates.
}
\item{START}{what to start the initial positions with. 
                "FR" for Fruchterman-Reingold. 
                "geodesic" for geodesic distances.
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
\seealso{vblpcmfit
}
\examples{
data(sampson)
### plot the mean posterior positions
plot(vblpcmfit(vblpcmstart(samplike,G=3)),main="Sampson's Monks: VB Solution")
### plot the mean posterior positions with larger node sizes
plot(vblpcmfit(vblpcmstart(samplike,G=3)),R2=1,main="Sampson's Monks: VB Solution")
}
