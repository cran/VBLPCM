\name{Y_to_nonE}
\alias{Y_to_nonE}
\title{calculate a non-edge list from an adjacency matrix
}
\description{uses C code to quickly calculate all non-edges as a two column matrix given an
adjacency matrix. i.e. all zeros in the adjacency matrix will correspond to a row in the
non-edgelist nonE
}
\usage{Y_to_nonE(N, directed, Y)
}
\arguments{
  \item{N}{number of nodes
}
  \item{directed}{logical indicator of directedness; TRUE=>directed FALSE=>undirected
}
  \item{Y}{input adjacency matrix
}
}
\value{A matrix of the non-edges with N^2-NE rows and 2 columns 
where NE is the number of edges. 
}
\author{Michael Salter-Townshend
}

\seealso{Y_to_E, Y_to_M, E_to_Y
}
