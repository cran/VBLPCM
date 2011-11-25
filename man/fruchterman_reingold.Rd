\name{fruchterman_reingold}
\alias{fruchterman_reingold}
\title{Perform Fruchterman-Reingold layout of a network in 2 or more dimensions.
}
\description{This was written and incorporated into the VBLPCM package because the
Fruchterman-Reingold routine in the network package only works in two dimensions. 
}
\usage{fruchterman_reingold(net, D=2, steps=1000, repulserad=N^D, m=N*(D-1), 
                            volume=N*(D-1))
}
\arguments{
\item{net}{network object on which to perform Fruchterman-Reingold layout.
}
\item{D}{Desired dimension of the space in which to lay out the network.
}
\item{steps}{Number of desired iterations. 
}
\item{repulserad}{The radius at which repulsion and attraction of linked nodes are equal. 
}
\item{m}{The maximum change in position per iteration. 
}
\item{volume}{The volume of space within which to position the nodes. 
}
}
\value{An N*D matrix of coordinates. 
}
\author{Michael Salter-Townshend
}

\seealso{log_like_forces
}
\examples{
### 2D example
### load the aids blogs dataset
data(aids)
### perform the Fruchterman-Reingold layout
X<-fruchterman_reingold(aids.net, D=2, steps=1e3)
### plot the results
plot(X)

### 3D example
### load the aids blogs dataset
data(aids)
### perform the Fruchterman-Reingold layout
X<-fruchterman_reingold(aids.net, D=3, steps=1e3)
### Not run
### plot the results in 3D
# library(rgl)
# plot3d(X)
}
