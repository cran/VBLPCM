\name{E_to_Y}
\alias{E_to_Y}
\title{create an adjacency matrix from an edgelist. 
}
\description{uses a call to C to transform edgelist to adjacency matrix. 
}
\usage{E_to_Y(N, NE, directed, E)
}
\arguments{
  \item{N}{number of nodes
}
  \item{NE}{number of edges
}
  \item{directed}{logical indicator of directedness (TRUE=>directed, FALSE=>undirected
}
  \item{E}{the input edgelist
}
}
\value{NxN sociomatrix / adjacency matrix
}
\author{Michael Salter-Townshend
}

\seealso{sociomatrix, Y_to_E
}
\keyword{sociomatrix}
