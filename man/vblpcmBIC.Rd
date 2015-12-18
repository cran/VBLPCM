\name{vblpcmbic}
\alias{vblpcmbic}
\title{calculate the BIC for the fitted VBLPCM object
}
\usage{vblpcmbic(v.params)
}
\description{
calculate the BIC for the fitted VBLPCM object
}
\arguments{
  \item{v.params}{The fitted values; output from vblpcmfit()
}
}
\details{BIC = BIC(edges | positions) + BIC(positions | clusters)
w/ BIC(edges | positions) = -2 loglikelihood + (P+1)log(number of edges) 
and BIC(positions | clusters) as per mclust
}
\value{The scalar value of the BIC 
}
\references{
Mark S. Handcock, Adrian E. Raftery and Jeremy Tantrum (2007).
     "Model-Based Clustering for Social Networks."
     Journal of the Royal Statistical Society: Series A (Statistics in Society), 170(2), 301-354.
}
\author{Michael Salter-Townshend
}

\seealso{latentnet::summary.ergmm
}
\keyword{BIC}
\examples{
data(sampson)
set.seed(1)
### plot the BIC for G=2,3,4 groups 
gbic<-list(groups=NULL,bic=NULL)
for (g in 2:4)
  {
  v.fit<-vblpcmfit(vblpcmstart(samplike,G=g,LSTEPS=1e3),STEPS=20)
  gbic$groups[g]=v.fit$G
  gbic$bic[g]=v.fit$BIC$overall
  }
plot(gbic$groups, gbic$bic, main="BIC results", xlab="# groups", ylab="BIC", t='b')
}
