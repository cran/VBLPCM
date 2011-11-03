vblpcmcovs<-function(N, model, Y, edgecovs=NULL,nodecovs=NULL)
  {
  XX<-matrix(rep(1,N^2),ncol=1) # all get the intercept term
  P=ncol(XX)
  if (model=="receiver")
    {
    # receiver random effects
    #P=P+N
    #XX<-cbind(XX,t(matrix(rep(diag(1,N),N),N)))
    P=P+1
    tmp<-apply(Y,2,sum,na.rm=1)
    tmp<-(tmp-mean(tmp))/sd(tmp)
    XX<-cbind(XX,rep(tmp,N))
    }
  
  if (model=="sender")
    {
    # sender random effects
    #XX<-cbind(XX,matrix(0,N^2,N))
    #for (i in 1:N)
    #  XX[((i-1)*N+1):(i*N),i+P]<-1
    #P=P+N
    P=P+1
    tmp<-apply(Y,1,sum,na.rm=1)
    tmp<-(tmp-mean(tmp))/sd(tmp)
    XX<-cbind(XX,c(t(matrix(rep(tmp,N),N))))
    }
  
  if (model=="social")
    {
    # sender random effects
    #XX<-cbind(XX,matrix(0,N^2,N))
    #for (i in 1:N)
    #  XX[((i-1)*N+1):(i*N),i+P]<-1
    #P=P+N
    #XX[,(P-N+1):P]<-XX[,(P-N+1):P]+t(matrix(rep(diag(1,N),N),N))
    #XX[,(P-N+1):P][XX[,(P-N+1):P]>1]<-1
    tmp1<-apply(Y,1,sum,na.rm=1)
    tmp1<-(tmp1-mean(tmp1))/sd(tmp1)
    tmp2<-apply(Y,2,sum,na.rm=1)
    tmp2<-(tmp2-mean(tmp2))/sd(tmp2)
    P=P+2
    XX<-cbind(XX,c(t(matrix(rep(tmp1,N),N))),rep(tmp2,N))
    }
  if (!is.null(nodecovs))
    {
    tmp<-expand.grid(1:N,1:N)
    nodecovs<-as.matrix(nodecovs)
    nodeedgecovs<-nodecovs[tmp[,1],]-nodecovs[tmp[,2],]
    }
  if (!is.null(nodecovs))
    XX<-cbind(XX,nodeedgecovs)
  if (!is.null(edgecovs))
    XX<-cbind(XX,edgecovs)
  return(as.matrix(XX))
  }
  
  
