vblpcmKL<-function(x)
  {
  P<-x$P
  d<-x$d
  N<-x$N
  NE<-x$NE
  NnonE<-x$NnonE
  NM<-x$NM
  G<-x$G
  Y<-x$Y
  E<-x$E
  nonE<-x$nonE
  M<-x$M
  numedges<-x$numedges
  EnonE<-x$EnonE
  diam<-x$diam
  hopslist<-x$hopslist
  XX<-x$XX
  V_xi<-x$V_xi
  V_psi2<-x$V_psi2
  V_z<-x$V_z
  V_sigma2<-x$V_sigma2
  V_eta<-x$V_eta
  V_lambda<-x$V_lambda
  V_omega2<-x$V_omega2
  V_nu<-x$V_nu
  V_alpha<-x$V_alpha
  xi<-x$xi
  psi2<-x$psi2
  sigma02<-x$sigma02
  omega2<-x$omega2
  nu<-x$nu
  alpha<-x$alpha
  inv_sigma02<-x$inv_sigma02
  if (!exists("x$STRAT")) # required but not provided in values from vblpcmstart
    x$STRAT=1
  STRAT=x$STRAT
  KL=0
  total_KL<-function(P, D, N, NE, NnonE, NM, G, Y, E, nonE, M, numedges, EnonE, diam, hopslist, XX, V_xi, V_psi2, V_z, V_sigma2, V_eta, V_lambda, 
                     V_omega2, V_nu, V_alpha, xi, psi2, sigma02, omega2, nu, alpha, inv_sigma02, STRAT, KL) 
  		   {
                     ans<-.C("KL_total", NAOK=TRUE, P=as.integer(P),D=as.integer(d), N=as.integer(N), 
  		     NE=as.integer(NE), NnonE=as.integer(NnonE), NM=as.integer(NM),
                     G=as.integer(G), Y=as.numeric(t(Y)), E=as.integer(t(E)), nonE=as.integer(t(nonE)), M=as.integer(t(M)),
                     numedges=as.integer(t(numedges)), EnonE=as.integer(t(EnonE)),
		     diam=as.integer(diam), hopslist=as.integer(t(hopslist)),
                     XX=as.double(t(XX)), V_xi=as.double(V_xi), V_psi2=as.double(V_psi2), V_z=as.double(t(V_z)),
                     V_sigma2=as.double(V_sigma2), V_eta=as.double(t(V_eta)),
                     V_lambda=as.double(t(V_lambda)),
                     V_omega2=as.double(V_omega2), V_nu=as.double(V_nu), V_alpha=as.double(V_alpha),
                     xi=as.double(xi), psi2=as.double(psi2), sigma02=as.double(sigma02),
                     omega2=as.double(omega2), nu=as.double(nu), alpha=as.double(alpha),
                     inv_sigma02=as.double(inv_sigma02), dists=as.double(t(as.matrix(dist(V_z)))),
  		     STRAT=as.double(STRAT), KL=as.double(KL), PACKAGE="VBLPCM")
                     return(ans)
                     }
  final_KL<-total_KL(P, d, N, NE, NnonE, NM, G, Y, E, nonE, M, numedges, EnonE, diam, hopslist, XX, V_xi, V_psi2, V_z, V_sigma2, V_eta, 
                     V_lambda, V_omega2, V_nu, V_alpha, xi, psi2, sigma02, omega2, nu, alpha,
  		   inv_sigma02, STRAT, KL)
  cat("KL distance to true posterior is ", final_KL$KL, "+ constant \n")
  final_KL$KL
  }

