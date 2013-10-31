\name{Y_to_nonE}
\alias{Y_to_nonE}
\title{calculate a non-edge list from an adjacency matrix
}
\description{uses C code to quickly calculate all non-edges as a two column matrix given an
adjacency matrix. i.e. all zeros in the adjacency matrix will correspond to a row in the
non-edgelist nonE
}
\usage{Y_to_nonE(N, NnonE, directed, Y)
}
\arguments{
  \item{N}{number of nodes
}
  \item{NnonE}{number of non-edges
}
  \item{directed}{logical indicator of directedness; TRUE=>directed FALSE=>undirected
}
  \item{Y}{input adjacency matrix
}
}
\value{A matrix of the non-edges with NnonE rows and 2 columns 
where NnonE is the number of non-edges. 
}
\author{Michael Salter-Townshend
}

\seealso{Y_to_E, Y_to_M, E_to_Y
}
