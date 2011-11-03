\name{vblpcmKL}
\alias{vblpcmKL}
\title{print and returns the Kullback-Leibler divergence from the fitted vblpcm object to the true
LPCM posterior
}
\description{print and returns the Kullback-Leibler divergence from the fitted vblpcm object to the true
LPCM posterior
}
\usage{vblpcmKL(x)
}
\arguments{
  \item{x}{The fitted values; output from vblpcmfit() or vblpcmstart()
}
}
\author{Michael Salter-Townshend
}
\details{The normalising constant of the posterior is unknown and therefore 
the Kullback-Leibler divergence is missing a constant. 
}
