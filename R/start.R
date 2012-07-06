# functions for calculating adjacency matrices, edges, non-edges, etc from each other
#source("adjacency_to_edges.R") 
vblpcmstart<-function(g.network, G=1, d=2, LSTEPS=5e3, model="plain", CLUST=0,
                       lcc=TRUE,edgecovs=NULL,nodecovs=NULL,START="FR", seed=0)
  {
  if (!is.nan(seed))
    set.seed(seed) # use this to seed the random number generator in R
  directed=is.directed(g.network)
  if (!(model=="plain" | model=="rsender" | model=="rreceiver" | model=="rsocial"))
    {
    print("Error: unknown model") 
    print("Defaulting to plain model without social effects") 
    print("see vblpcmcovs for details") 
    model="plain"
    }
  if (!directed & (model=="rsender" | model=="rreceiver"))
    {
    print("Error: Model with directed effects cannot be applied to an undirected network") 
    print("Defaulting to plain model without social effects") 
    model="plain"
    }
  if (lcc)
    {
    all_g.network<-g.network
    # only look at largest connected component
    N=all_g.network$gal$n
    # can't have unknown links for the components function to work
    tmp<-as.sociomatrix(g.network)+t(as.sociomatrix(g.network))
    tmp[tmp>1]<-1
    tmp[is.na(tmp)]<-0 # fill NAs with zeros
    tmp<-component.largest(tmp)
    # leads to stack overflow in some cases:
    #g.network<-get.inducedSubgraph(all_g.network,(1:N)[component.largest(tmp)]) 
    # safer:
    g.network<-network(as.sociomatrix(all_g.network)[(1:N)[tmp],(1:N)[tmp]],directed=is.directed(all_g.network))
    for (att in list.vertex.attributes(all_g.network))
      set.vertex.attribute(g.network,att,get.vertex.attribute(all_g.network,att)[tmp])
    rm(all_g.network)
    }
  
  N=g.network$gal$n
  NE=network.edgecount(g.network)
  Y<-as.sociomatrix(g.network)
  E<-Y_to_E(N,NE,g.network$gal$directed,Y)
  nonE<-Y_to_nonE(N, directed, Y)
  NnonE=nrow(nonE)
  M<-Y_to_M(N, directed, Y)
  NM=nrow(M)
  # create edge/non-edge matrix for all nodes
  EnonE<-matrix(NaN,N,N)
  numedges<-matrix(NaN,N,2) # columns are #edges and #non-links
  numedges[,1]<-apply(Y,1,sum,na.rm=1)
  numedges[,2]<-apply(!is.na(Y),1,sum)-numedges[,1]
  for (i in 1:N)
    {
    if (numedges[i,1]!=0)
      EnonE[i,1:numedges[i,1]]<-((1:N)[Y[i,] == 1])[!is.na((1:N)[Y[i,] == 1])]     # edges
    EnonE[i,(numedges[i,1]+1):(N-sum(is.na(Y[i,])))]<-((1:N)[Y[i,]==0])[!is.na((1:N)[Y[i,] == 0])] # non-edges
    }
  tmpY<-Y+t(Y);tmpY[tmpY>1]<-1
  hops<-geodist(tmpY)$gdist
  hops[is.infinite(hops)]<-10+max(hops[is.finite(hops)]) # 1 further than the furthest actually observed
  #hops<-breadth_first_hops(Y,2,5)
  diam<-max(hops,na.rm=1)
  #hops[is.na(Y)]<-1+hops[is.na(Y)] # for the missing links
  # each row contains the numbers of nodes that are #hops away from the node indexed by the row
  # followed by all nodes grouped and ordered by this.
  hopslist<-hops_to_hopslist(hops,diam,N) 
  # create covariates design matrix
  XX<-vblpcmcovs(N,model,Y,edgecovs,nodecovs)
  XX_n<-XX$XX_n
  XX_e<-XX$XX_e
  P_n=ncol(XX_n)
  if (is.null(P_n))  
    {
    P_n=0
    XX_n=0
    }
  P_e=ncol(XX_e)
  
  # variational parameters are: 
  #1. V_z 
  #2. S
  #3. V_lambda
  #4. V_eta
  #5. V_omega2
  #6. V_alpha
  #7. V_nu
  #8. V_xi
  #9. V_psi2
  # initialise variational parameters
  
  sigma02=0.125
  inv_sigma02=1/sigma02
  
  if (START=="random")
    X<-matrix(runif(N*d,-2*d,2*d),ncol=d) # may lead to local minima
  if (START=="geodesic")
    X<-jitter(cmdscale(hops,k=d)) # can't rely on this for data with large social effects
  if (START=="laplace")
    {
    Tinv<-diag(degree(g.network)^-0.5)
    Tinv[is.infinite(Tinv)]<-0
    tmpY<-Y
    tmpY[is.nan(Y)]<-0
    L<-Tinv%*%(diag(degree(g.network))-tmpY%*%t(tmpY))%*%Tinv
    X<-jitter(cmdscale(L,k=d)) 
    }
  if (START=="FR")
    {
    #require(igraph)
    #X<-layout.fruchterman.reingold(graph.adjacency(Y))
    #detach("package:igraph")
    #X<-network.layout.fruchtermanreingold(g.network,layout.par=list(area=N)) # only works in 2D
    X<-fruchterman_reingold(g.network, d)
    }
  if (!exists("B")) 
    {
    delete<-seq(from=1, to=N*N, by=(N+1))
    y<-c(Y)[-delete]# logistic regression
    tmpx<- c(as.matrix(dist(X)))[-delete]
    p=mean(y,na.rm=1)
    B<-mean(tmpx)+log(p)-log(1-p) # "average distance" + log-odds(p) 
    }
  out<-log_like_forces(g.network, d, X, B, m=N, LSTEPS)
  initial_V_z<-out$X
  initial_V_xi_e<-out$B
  if (model=="rreceiver") 
    {
    initial_V_xi_n<-apply(Y,2,sum,na.rm=1)
    initial_V_xi_n<-(initial_V_xi_n-mean(initial_V_xi_n))/sd(initial_V_xi_n)
    }
  if (model=="rsender") 
    {
    initial_V_xi_n<-apply(Y,1,sum,na.rm=1)
    initial_V_xi_n<-(initial_V_xi_n-mean(initial_V_xi_n))/sd(initial_V_xi_n)
    }
  if (model =="rsocial") 
    {
    tmp1<-apply(Y,2,sum,na.rm=1)
    tmp1<-(tmp1-mean(tmp1))/sd(tmp1)
    tmp2<-apply(Y,1,sum,na.rm=1)
    tmp2<-(tmp2-mean(tmp2))/sd(tmp2)
    initial_V_xi_n<-c(t(matrix(c(tmp1, tmp2),N)))
    }
  if (d>1) 
    {
    fitmc<-summary(mclustBIC(initial_V_z,G=G,modelNames="VII"),initial_V_z)
    } else fitmc<-summary(mclustBIC(initial_V_z,G=G,modelNames="V"),initial_V_z)
  if (is.null(fitmc$bic))
    {
    print("Couldn't fit initial clustering using mclust")
    print("Try using a different initialisation or different number of clusters or latent space dimension")
    print("Using no groups for now")
    G=1
    if (d>1) 
      {
      fitmc<-summary(mclustBIC(initial_V_z,G=G,modelNames="VII"),initial_V_z)
      } else fitmc<-summary(mclustBIC(initial_V_z,G=G,modelNames="V"),initial_V_z)
    }
  initial_V_eta<-t(fitmc$parameter$mean)
  initial_V_omega2<-c(t(fitmc$parameter$variance$sigmasq)) # for use with EII
  if (G > 1) initial_V_lambda<-t(fitmc$z)
  if (G==1) initial_V_lambda<-matrix(1,1,N)
  # centering: 1 is use clusters, 0 is use original positions; anything in between is ok
  initial_V_z<-t(t(initial_V_eta)%*%initial_V_lambda)+
               (1-CLUST)*(initial_V_z-t(t(initial_V_eta)%*%initial_V_lambda)) 
  #if (G>1) initial_V_sigma2<-(fitmc$uncertainty+0.5*mean(fitmc$uncertainty)) 
  if (G>1) initial_V_sigma2<-rep(0.2,N)
  if (G==1) initial_V_sigma2<-rep(0.2,N)
  if (G>1) initial_V_nu<-fitmc$parameters$pro
  if (G==1) initial_V_nu<-0.5
  initial_V_psi2<-2e0
  initial_V_alpha<-1/fitmc$parameter$variance$scale
  
  conv_check=1
  
  V_z<-initial_V_z
  V_eta<-initial_V_eta
  V_lambda<-initial_V_lambda
  V_sigma2<-initial_V_sigma2
  V_omega2<-initial_V_omega2
  V_nu<-initial_V_nu
  V_alpha<-initial_V_alpha
  initial_V_xi_e<-c(initial_V_xi_e,rep(0,P_e-1))
  initial_V_psi2_e<-rep(initial_V_psi2,P_e) 
  if (P_n>0)
    {
    V_xi_n<-initial_V_xi_n
    V_psi2_n<-rep(initial_V_psi2,P_n)
    } else {V_xi_n<-NaN; V_psi2_n<-NaN}
  V_xi_e<-initial_V_xi_e
  V_psi2_e<-initial_V_psi2_e
  
  # priors
  xi=0
  mu_nought=rep(0,d)
  # priors from latentnet package...
  psi2=9.0e0
  omega2=0.75
  alpha=rep(4.5,G)*sigma02
  nu=rep(2.5,G)
  ###################################
  variational.start<-list()
  seed->variational.start$seed
  g.network->variational.start$net
  P_n->variational.start$P_n
  P_e->variational.start$P_e
  model->variational.start$model
  d->variational.start$d
  N->variational.start$N
  NE->variational.start$NE
  NnonE->variational.start$NnonE
  NM->variational.start$NM
  G->variational.start$G
  Y->variational.start$Y
  E->variational.start$E
  nonE->variational.start$nonE
  M->variational.start$M
  numedges->variational.start$numedges
  EnonE->variational.start$EnonE
  diam->variational.start$diam
  hopslist->variational.start$hopslist
  XX_n->variational.start$XX_n
  XX_e->variational.start$XX_e
  V_xi_n->variational.start$V_xi_n
  V_xi_e->variational.start$V_xi_e
  V_psi2_n->variational.start$V_psi2_n
  V_psi2_e->variational.start$V_psi2_e
  V_z->variational.start$V_z
  V_sigma2->variational.start$V_sigma2
  V_eta->variational.start$V_eta
  V_lambda->variational.start$V_lambda
  V_omega2->variational.start$V_omega2
  V_nu->variational.start$V_nu
  V_alpha->variational.start$V_alpha
  xi->variational.start$xi
  psi2->variational.start$psi2
  sigma02->variational.start$sigma02
  omega2->variational.start$omega2
  nu->variational.start$nu
  alpha->variational.start$alpha
  inv_sigma02->variational.start$inv_sigma02
  BIC<-vblpcmbic(variational.start)
  BIC->variational.start$BIC
  class(variational.start)<-"vblpcm"
  return(variational.start)
  }
