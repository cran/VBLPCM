\name{predict.vblpcm}
\alias{predict.vblpcm}
\title{Find all link probabilities
}
\description{generate a matrix of link probabilities based on the fitted VB model.
}
\usage{
\method{predict}{vblpcm}(object, \dots)
}
\arguments{
  \item{object}{The fitted values; output from vblpcmfit()
}
  \item{\dots}{optional additional arguments.
}
}
\value{The posterior predictive link probabilities given the fitted object
}
\author{Michael Salter-Townshend
}
\examples{
data(sampson)
v.fit<-vblpcmfit(vblpcmstart(samplike,G=3))
### create a matrix of link posterior probabilities given the fitted model
probs<-predict.vblpcm(v.fit)
# show this graphically; separation of the boxes implies a good fit to the data
boxplot(split(probs,v.fit$Y),
        ylab=expression(paste("P(",Y[i][j],"=1)")),xlab=expression(paste(Y[i][j])))
}
