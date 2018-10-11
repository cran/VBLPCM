\name{VBLPCM-package}
\alias{VBLPCM-package}
\alias{VBLPCM}
\docType{package}
\title{VBLPCM: Variational Bayes for the Latent Position Cluster Model for networks
}
\description{A faster approximate alternative to using latentnet. Interfaces C code to fit a
Variational Bayes approximation to the posterior for the Latent Position Cluster Model for networks. 
}
\details{
\tabular{ll}{
Package: \tab VBLPCM\cr
Type: \tab Package\cr
Version: \tab 2.4.5\cr
Date: \tab 2018-10-09\cr
License: \tab GPL (>=2) \cr
LazyLoad: \tab yes\cr
}
This package is designed to be used as an alternative to the latentnet package when network
size computationally prohibits latentnet. It uses a Variational Bayesian Expectation Maximisation algorithm
to compute a closed-form approximation to the posterior that the ergmm function in latentnet samples
from. It may be thought of as an intermediary approximation that is more accurate than
the two-stage MLE fit provided by latentnet but a faster approximation to the MCMC sampler provided
by latentnet. In fact, the VB iterations also converge quicker than the two-stage MLE. 

VBLPCM can also take advantage of the stratified sampler of Adrian Raftery, Xiaoyue Niu, Peter
Hoff and Ka Yee Yeung. This approximation to the (log)likelihood allows for even larger networks to
be analysed (see tech report below). Rather than using a fixed number of "controls" per geodesic
distance we set a probability of sampling each non-link at each level. 

We also provide four choices of model; these are "plain" and three with random node-specific social effects. 
"rsender" for sender random effects, "rreceiver" for receiver random effects and "rsocial" for both. 
For undirected networks only "plain" or "rsocial" may be chosen. 
}
\references{
Michael Salter-Townshend and Thomas Brendan Murphy (2013).
"Variational Bayesian Inference for the Latent Position Cluster Model."
Computational Statistics and Data Analysis, volume 57, number 1, pages 661-671.
DOI=10.1016/j.csda.2012.08.004


Pavel N. Krivitsky and Mark S. Handcock (2008). 
"Fitting Latent Cluster Models for Social Networks with latentnet."
Journal of Statistical Software, number 5, volume 24, pages 1-23. 


Mark S. Handcock, Adrian E. Raftery and Jeremy Tantrum (2007).
"Model-Based Clustering for Social Networks."  
Journal of the Royal Statistical Society: Series A (Statistics in Society), 170(2), 301-354.


Adrian Raftery, Xiaoyue Niu, Peter Hoff and Ka Yee Yeung (2012). 
"Fast Inference for the Latent Space Network Model Using a Case-Control Approximate Likelihood."
Journal of Computational and Graphical Statistics. 
doi: 10.1080/10618600.2012.679240 


Sucharita Gopal (2007). 
"The Evolving Social Geography of Blogs"
Societies and Cities in the Age of Instant Access Berlin:Springer, 275--294
}
\keyword{graphs}
\keyword{networks}
\keyword{variational}
\keyword{package}
\keyword{models}
\keyword{latentnet}
\seealso{
\code{\link[VBLPCM]{vblpcmstart}}
\code{\link[VBLPCM]{vblpcmfit}}
}
\examples{
### Sampson's monks with sender random effects ###
data(sampson,package="VBLPCM")
v.start<-vblpcmstart(samplike,G=3,model="rreceiver",LSTEPS=1e3)
v.fit<-vblpcmfit(v.start,STEPS=20)
### plot the mean posterior positions ###
plot(v.fit, R2=0.05,main="Sampson's Monks: VB with Receiver Random Effects")
### Who's in each group?  ###
vblpcmgroups(v.fit)

### Look at a goodness-of-fit plot ###
plot(gof(v.fit,GOF=~distance))

### create a matrix of link posterior probabilities given the fitted model ###
probs<-predict.vblpcm(v.fit)
### create a boxplot goodness-of-fit graphic ###
boxplot(split(probs,as.sociomatrix(samplike)))

### run a bigger example, using the likelihood sampler set to 0.1 ###
### ## Not run:
##  data(aids,package="VBLPCM")
##  v.start<-vblpcmstart(aids.net,G=7,model="rsender",d=3)
##  use the case-control sampler with 10 controls per case
##  v.fit<-vblpcmfit(v.start,NC=10)
### plot the mean posterior positions ###
##  plot(v.fit, R2=0.1,main="Aids Blogs with Sender Random Effects")

### Use ROC / AUC to get a measure of model fit to the data ###
##  vblpcmroc(v.fit)
## End(Not run)
}
