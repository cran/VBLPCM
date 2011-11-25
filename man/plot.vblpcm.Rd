\name{plot.vblpcm}
\alias{plot.vblpcm}
\title{plot the posterior latent positions and groupings and network
}
\description{Plot the network using the estimated positions with clustering. The nodes are plotted
as pie-charts to show group membership probabilities. The group means are coloured crosses and the group
standard deviations are shown with coloured circles. 
}
\usage{
\method{plot}{vblpcm}(x, \dots, R2 = 0.2, main = "Variational-Bayes Positions", 
                   alpha = 0.5, colours=1:x$G, RET=FALSE)
}
\arguments{
  \item{x}{The fitted values; output from vblpcmfit()
}
  \item{\dots}{optional arguments to be passed to lower level functions
}
  \item{R2}{scaling factor for the size of each node in the plot
}
  \item{main}{main title for the plot
}
  \item{alpha}{transparency of the links
}
  \item{colours}{colours of the groups
}
  \item{RET}{whether to return the 2D postions of nodes and clusters
}
}
\details{Plots the latent positions and clustering of a network fitted via vblpcmfit() or vblpcmstart()

Each node appears in the latent space as a pie chart with segments size proportional to group memberships. 
The clusters are represented as circles in the latent space centred on the expected position of the group mean 
and with size proportional to the cluster standard deviation. 

If applicable, the size of the pie charts represents the expected sociality effect of the node. 
}
\author{Michael Salter-Townshend
}

\seealso{latentnet::plot.ergmm
}
