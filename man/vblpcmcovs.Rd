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
\usage{vblpcmcovs(N, model, Y, edgecovs=NULL, sendcovs=NULL, receivecovs=NULL,
                  socialcovs=NULL)
}
\arguments{
  \item{N}{number of nodes
}
  \item{model}{model; may be "plain", "rreceiver", "rsender" or "rsocial". See Details. 
}
  \item{Y}{adjacency matrix
}
  \item{edgecovs}{optional additional covariate / attribute data on the edges
}
  \item{sendcovs}{optional additional covariate / attribute data on the nodes for links out
}
  \item{receivecovs}{optional additional covariate / attribute data on the nodes for links in
} 
  \item{socialcovs}{optional additional covariate / attribute data on the nodes for links in and out
}
}
\details{Can be used to construct design matrices with edge covariates or node covariates 
and / or sociality effects. "rreceiver", "rsender" and "rsocial" model random social effects. 
Node covariates are differenced and treated as edge covariates. 
}
\value{An edge design matrix that is Pe x N^2 and a node design matrix that is Pn x N
where Pe is the number of edge covariates and Pn is the number of node covariates. 
}
\author{Michael Salter-Townshend
}

\seealso{vblpcmstart
}
\keyword{covariates}
\keyword{sociality effects}
