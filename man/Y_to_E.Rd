\name{Y_to_E}
\alias{Y_to_E}
\title{calculate the edgelist for a given adjacency matrix
}
\description{calls C code to quickly transform from adjacency to edgelist
}
\usage{Y_to_E(N, NE, directed, Y)
}
\arguments{
  \item{N}{number of nodes
}
  \item{NE}{number of edges
}
  \item{directed}{logical indicator of directedness; TRUE=>directed FALSE=>undirected
}
  \item{Y}{input adjacency matrix
}
}
\value{An edgelist matrix E of size NE x 2
}
\author{Michael Salter-Townshend
}

\seealso{edgelist, E_to_Y
}
