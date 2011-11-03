predict.vblpcm<-function(object, ...)
  {
  d<-object$d
  N<-object$N
  XX<-object$XX
  V_xi<-object$V_xi
  #P=object$P;V_xi[2:P]<- -V_xi[2:P]
  V_psi2<-object$V_psi2
  V_z<-object$V_z
  V_sigma2<-object$V_sigma2
  DIST<-c(sqrt(as.matrix(dist(V_z, diag=1, upper=1)^2) + matrix(d*apply(expand.grid(V_sigma2,V_sigma2),1,sum),N))) 
  #tmp<-(XX%*%V_xi-DIST)
  #probs<-matrix(exp(tmp-log(1+exp(tmp+0.5*XX%*%V_psi2))),N)
  # the next two lines are because the exponent of the expected log-likelihood is not the same as the expected likelihood
  tmp<-(XX%*%V_xi-DIST)+0.5*XX%*%V_psi2
  probs<-matrix(exp(tmp-log(1+exp(tmp))),N) 
  return(probs)
  }
