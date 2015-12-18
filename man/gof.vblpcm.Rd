\name{gof.vblpcm}
\alias{gof.vblpcm}
\title{Goodness of fit based on simulations from the fitted object. 
}
\description{Create a goodness of fit statistics and plots based on the degree distributions of networks simulated fitted 
from a fitted variational approximation. 
}
\usage{
\method{gof}{vblpcm}(object, \dots,
         nsim=100,
         GOF=NULL,
         verbose=FALSE)
}
\arguments{
  \item{object}{fitted VBLPCM object; usually output from vblpcmfit() or vblpcmstart()
}
  \item{\dots}{optional arguments for lower level functions
}
  \item{nsim}{number of networks to simulate
}
 \item{GOF}{formula; an \R formula object, of the form
    \code{~ <model terms>} specifying the
    statistics to use to diagnosis the goodness-of-fit of the model.
    They do not need to be in the model formula specified in
    \code{formula}, and typically are not.
    Examples are the degree distribution ("degree"),
    minimum geodesic distances ("dist"), and shared partner distributions
    ("espartners" and "dspartners").
    For the details on the possible
    \code{<model terms>}, see \code{\link[ergm]{ergm-terms}}.
} 
 \item{verbose}{Provide verbose information on the progress of the
    simulation.
} 
}
\details{A sample of graphs is randomly drawn from the posterior of the vblpcmfit() result. 

     A plot of the summary measures may then be plotted using plot().  
}
\author{Michael Salter-Townshend
}

\seealso{latentnet::gof.ergmm
}
\keyword{degree distribution}
\examples{
data(sampson,package="VBLPCM")
v.start<-vblpcmstart(samplike,G=3,model="rreceiver",LSTEPS=1e3)
v.fit<-vblpcmfit(v.start,STEPS=20)
### plot the mean posterior positions
plot(v.fit, R2=0.05,main="Sampson's Monks: VB with Receiver Effects")
### Look at gof plots
plot(gof(v.fit,GOF=~distance,nsim=50))
}
