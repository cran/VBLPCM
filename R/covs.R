vblpcmcovs<-function(N, model, Y, edgecovs=NULL,nodecovs=NULL)
  {
  XX_n<-NULL
  XX_e<-matrix(rep(1,N^2),ncol=1) # all get the intercept term
  if (model=="rreceiver")
    {
    # receiver random effects
    XX_n<-cbind(XX_n,rep(1,N))
    }
  if (model=="rsender")
    {
    # sender random effects
    XX_n<-cbind(XX_n,rep(1,N))
    }
  if (model=="rsocial")
    {
    # sender random effects
    XX_n<-cbind(XX_n,rep(1,N),rep(1,N))
    }
  if (!is.null(nodecovs)) # include option to not model nodecovs as edgecovs?
    {
    tmp<-expand.grid(1:N,1:N)
    nodecovs<-as.matrix(nodecovs)
    nodeedgecovs<-nodecovs[tmp[,1],]-nodecovs[tmp[,2],]
    }
  if (!is.null(nodecovs))
    XX_n<-cbind(XX_n,nodeedgecovs)
  if (!is.null(edgecovs))
    XX_e<-cbind(XX_e,edgecovs)
  return(list("XX_n"=XX_n,"XX_e"=XX_e))
  }
  
  
