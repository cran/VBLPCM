\name{vblpcmroc}
\alias{vblpcmroc}
\title{ROC curve plot for vblpcmfit
}
\description{Plot a Receiver Operating Curve to show model fit in terms of link prediction.
}
\usage{vblpcmroc(v.params, NUM=100)
}
\arguments{
  \item{v.params}{The fitted values; output from vblpcmfit()
}
  \item{NUM}{The number of intervals on the roc curve
}
}
\details{A threshold is varied between zero and one. At each point the probability of a link 
between all pairs of nodes is calculated on the v.params argument containing a fitted vblpcm object. 
If greater than the threshold the link is "predicted" present, else it is "predicted" absent. A plot 
of the proportion of true and false positives for each threshold value is thus obtained. 
}
\value{The Area Under the Curve (AUC). The closer to 1 the better the fit. 
}
\author{Michael Salter-Townshend
}
