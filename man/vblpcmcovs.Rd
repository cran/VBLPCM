\name{vblpcmcovs}
\alias{vblpcmcovs}
\title{create the design matrix for the network analysis
}
\description{
Add intercept (column of ones) and degree-based covariates (if model is for sociality effects)
to a user-supplied (default is NULL) edge covariates matrix of size N^2 rows and C columns where
C is the number of covariates. Node covariates may be converted to difference-between-pairs for
edges. 
}
\usage{vblpcmcovs(N, model, Y, edgecovs=NULL, nodecovs=NULL)
}
\arguments{
  \item{N}{number of nodes
}
  \item{model}{model; may be "plain", "receiver", "sender" or "social"
}
  \item{Y}{adjacency matrix
}
  \item{edgecovs}{optional additional covariate / attribute data on the edges
}
  \item{nodecovs}{optional additional covariate / attribute data on the nodes
}
}
\details{Can be used to construct design matrices with edge covariates, node covariates 
and / or sociality effects. Node covariates will be differences and treated as edge covariates. 
}
\value{A design matrix that is P x N^2
}
\author{Michael Salter-Townshend
}

\seealso{vblpcmstart
}
\keyword{covariates}
\keyword{sociality effects}
