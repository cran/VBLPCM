\name{hops_to_hopslist}
\alias{hops_to_hopslist}
\title{create a handy matrix of vectors to store the hopslist
}
\description{Designed for nternal use only; store the geodesic distances in a handy format Each node gets a vector in the hopslist matrix. 
Each row describes a node and for each row: 
The first diam entries state the number of nodes that are that distance away by shortest path where diam is the 
maximum shortest path between two nodes (the graph diameter). 
eg if entry 3 in row 4 is a 5 then there are exactly 5 nodes that are 4 hops away from node 3. 
This vector is followed by the indices of all the nodes, grouped by the length of the shortest paths. 
}
\usage{hops_to_hopslist(hops, diam, N)
}
\arguments{
  \item{hops}{matrix of geodesic distances
}
  \item{diam}{diameter of the network
}
  \item{N}{total number of nodes in the network
}
}
\author{Michael Salter-Townshend
}

