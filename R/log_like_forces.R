log_like_forces<-function(net, D, P, X, B, XX, m=network.size(net), steps=1e3)
  {
  directed<-is.directed(net)
  N<-network.size(net)
  Y<-as.sociomatrix(net)
  C_log_like_forces<- function(directed, N, D, P, steps, Y, X, B, XX, m)
     {
     ans<-.C("log_like_forces", NAOK=TRUE, directed=as.integer(directed), N=as.integer(N), 
             D=as.integer(D), P=as.integer(P), steps=as.integer(steps), 
             Y=as.double(t(Y)), X=as.numeric(t(X)), B=as.numeric(B), XX=as.numeric(t(XX)), m=as.numeric(m),
	     PACKAGE="VBLPCM")
             return(ans)
     }
  
  delete<-seq(from=1, to=N*N, by=(N+1))
  y<-c(Y)[-delete]# logistic regression
  y[is.na(y)]<-0
  loglike<-function(Beta, x, y) 
    {
    covs<-(XX%*%Beta)[-delete]
    sum(y*(covs-x)) - sum(log(1+exp(covs-x)))
    }
  
  if (!exists("doB")) doB<-1
  
  for (i in 1:1)
    {
    # update B
    tmpx<-c(as.matrix(dist(X)))[-delete]
    if (doB==1)
      B<-optim(B, loglike, x=tmpx, y=y, method="BFGS", control=list(fnscale=-1))$par
    # update X
    out<-C_log_like_forces(directed, N, D, P, steps, Y, X, B, XX, m=N)
    out$X<-t(matrix(out$X,ncol=N))
    # centre
    out$X <- out$X - t(matrix(rep(apply(out$X,2,mean),N),nrow=D))

    out$B<-B
    }
  return (out)
  }
