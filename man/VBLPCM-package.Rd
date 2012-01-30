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
Version: \tab 1.3\cr
Date: \tab 2012-01-30\cr
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

Note that in this package (unlike in latentnet) we choose to model sociality effects as fixed rather
than random. We use the degree of a node as the effect size and there is therefore a single scalar
parameter controlling each sociality effect modelled for the whole network. This choice resulted
from the observation that the correlation between the social random effects in latentnet and the
degree of the corresponding nodes was very high. 

}
\references{
Michael Salter-Townshend and Thomas Brendan Murphy (2009).
"Variational Bayesian Inference for the Latent Position Cluster Model."
Workshop on Analyzing Networks and Learning with Graphs. 
Neural Information Processing Systems. 


Pavel N. Krivitsky and Mark S. Handcock (2008). 
"Fitting Latent Cluster Models for Social Networks with latentnet."
Journal of Statistical Software, number 5, volume 24, pages 1-23. 


Mark S. Handcock, Adrian E. Raftery and Jeremy Tantrum (2007).
     "Model-Based Clustering for Social Networks."  
     Journal of the Royal Statistical Society: Series A (Statistics in Society), 170(2), 301-354.


Adrian Raftery, Xiaoyue Niu, Peter Hoff and Ka Yee Yeung (2010). 
"Fast Inference for the Latent Space Network Model Using a Case-Control Approximate Likelihood."
University of Washington, Department of Statistics, Technical Report 572. 

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
### simulated data example with missing edges ###
rm(list=ls())
data(simulated.network,package="VBLPCM")
g<-network(Y, directed=FALSE)
v.start<-vblpcmstart(g,G=3)
v.fit<-vblpcmfit(v.start)
### plot the mean posterior positions ###
plot(v.fit, main="simulated network with missing edges: VB Solution") 
### Use ROC / AUC to get a measure of model fit to the data ###
vblpcmroc(v.fit)

### Sampson's monks ###
rm(list=ls())
data(sampson,package="VBLPCM")
v.start<-vblpcmstart(samplike,G=3)
v.fit<-vblpcmfit(v.start)
### plot the mean posterior positions ###
plot(v.fit, main="Sampson's Monks: VB Solution")
### Who's in each group?  ###
vblpcmgroups(v.fit)
### Use ROC / AUC to get a measure of model fit to the data ###
vblpcmroc(v.fit)

### create a matrix of link posterior probabilities given the fitted model ###
probs<-predict.vblpcm(v.fit)
### create a boxplot goodness-of-fit graphic ###
boxplot(split(probs,as.sociomatrix(samplike)))

### Sampson's monks with fixed social effects ###
rm(list=ls())
data(sampson,package="VBLPCM")
v.start<-vblpcmstart(samplike,G=3,model="sender")
v.fit<-vblpcmfit(v.start)
### plot the mean posterior positions ###
plot(v.fit, R2=0.02,main="Sampson's Monks: VB with Sender Effects")

### Use ROC / AUC to get a measure of model fit to the data ###
vblpcmroc(v.fit)
### Look at gof plots ###
plot(gof(v.fit))

### find the most probable link that does not in fact occur ###
monks<-network.vertex.names(samplike)
# save adjacency matrix to Y
Y<-as.sociomatrix(samplike)
probs<-predict.vblpcm(v.fit)
expand.grid(monks,monks)[Y==0,][order(probs[Y==0],decreasing=TRUE),][1,]
### this is a self link; let's find the most probable non-link that's not a self link ###
diag(Y)<-NaN
expand.grid(monks,monks)[Y==0,][order(probs[Y==0],decreasing=TRUE),][1,]

### find the least probable link that does in fact occur ###
expand.grid(monks,monks)[Y==1,][order(probs[Y==1]),][1,]

### run a bigger example, using the likelihood sampler set to 0.2 ###
rm(list=ls())
data(aids,package="VBLPCM")
v.start<-vblpcmstart(aids.net,G=7,model="sender")
v.fit<-vblpcmfit(v.start, STRAT=0.2)
### plot the mean posterior positions ###
plot(v.fit, R2=0.05,main="Aids Blogs with Sender Effects")

### Use ROC / AUC to get a measure of model fit to the data ###
vblpcmroc(v.fit)
}
