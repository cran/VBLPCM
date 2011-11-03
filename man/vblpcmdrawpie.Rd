\name{vblpcmdrawpie}
\alias{vblpcmdrawpie}
\title{add a piechart of group memberships of a node to a network plot; taken mainly from latentnet equivalent
}
\usage{vblpcmdrawpie(center,radius,probs,n=50,colours=1:length(probs))
}
\description{add a piechart of group memberships of a node to a network plot; taken mainly from latentnet equivalent
}
\arguments{
  \item{center}{where to postion the piechart
}
  \item{radius}{radius of the piechart / node
}
  \item{probs}{probability vector of cluster memberships
}
  \item{n}{order of polygon to approximate a circle
}
  \item{colours}{the colours used; default is from palette()
}
}
\author{Michael Salter-Townshend
}
\note{Thanks to Pavel N. Krivitsky of the latentnet package as I copied this from there. 
}

\seealso{plot.vblpcm
}
