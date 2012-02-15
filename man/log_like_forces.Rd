\name{log_like_forces}
\alias{log_like_forces}
\title{create an initial configuration for the latent positions. 
}
\description{This performs an iterative relaxation type algorithm to approximately find the
positions of the nodes in the latent space that maximises the log-likelihood. 
}
\usage{log_like_forces(net, D, X, B, m ,steps)
}
\arguments{
  \item{net}{network object on which to perform layout. 
}
  \item{D}{dimension of the latent space.
}
  \item{X}{the initial guess for X
}
  \item{B}{the intercept term.
}
  \item{m}{usually N will suffice. 
}
  \item{steps}{maximum number of iteration steps. 
}
}
\details{
Usually only used internally in vblpcmstart()
}
\value{Matrix of latent positions X 
}
\author{Michael Salter-Townshend
}
\seealso{igraph::layout.fruchterman.reingold
}
\examples{
data(sampson)
N=network.size(samplike)
X=matrix(runif(N*2,-2,2),ncol=2)
XX=vblpcmcovs(N,"plain",as.sociomatrix(samplike))
out<-log_like_forces(samplike, 2, X, 0, m=N, 1e3)
plot(samplike,coord=out$X)
}
