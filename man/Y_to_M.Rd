\name{Y_to_M}
\alias{Y_to_M}
\title{calculate the missing edges as an edgelist from an adjacency matrix with NaNs indicating
missing links
}
\description{uses C code to quickly find all pairs of nodes for which we do not know whether there
is a link or not, given an adjacency matrix with NaNs indicating unknown / unobserved linkage
}
\usage{Y_to_M(N, NM, directed, Y)
}
\arguments{
  \item{N}{number of nodes
}
  \item{NM}{number of missing edges
}
  \item{directed}{logical indicator of directedness; TRUE=>directed FALSE=>undirected
}
  \item{Y}{input adjacency matrix
}
}
\value{A matrix of missing edges M
}
\author{Michael Salter-Townshend
}

\seealso{Y_to_E, E_to_Y, Y_to_nonE
}
