vblpcmfit<-function(variational.start, STEPS=10, maxiter=100, tol=1e-6, STRAT=1, seed=NaN, d_vector=rep(TRUE,9))
  {
  if (length(d_vector)!=9)
    stop("You must supply a d_vector of length 9. Please refer to the help file for vblpcmfit\n")
  P<-variational.start$P
  model<-variational.start$model
  d<-variational.start$d
  N<-variational.start$N
  NE<-variational.start$NE
  NnonE<-variational.start$NnonE
  NM<-variational.start$NM
  G<-variational.start$G
  Y<-variational.start$Y
  E<-variational.start$E
  nonE<-variational.start$nonE
  M<-variational.start$M
  numedges<-variational.start$numedges
  EnonE<-variational.start$EnonE
  diam<-variational.start$diam
  hopslist<-variational.start$hopslist
  XX<-variational.start$XX
  V_xi<-variational.start$V_xi
  V_psi2<-variational.start$V_psi2
  V_z<-variational.start$V_z
  V_sigma2<-variational.start$V_sigma2
  V_eta<-variational.start$V_eta
  V_lambda<-variational.start$V_lambda
  V_omega2<-variational.start$V_omega2
  V_nu<-variational.start$V_nu
  V_alpha<-variational.start$V_alpha
  xi<-variational.start$xi
  psi2<-variational.start$psi2
  sigma02<-variational.start$sigma02
  omega2<-variational.start$omega2
  nu<-variational.start$nu
  alpha<-variational.start$alpha
  inv_sigma02<-variational.start$inv_sigma02
  conv=0 # not converged to start with
  out<-.C("Rf_VB_bbs", NAOK=TRUE, steps=as.integer(STEPS), max_iter=as.integer(maxiter), P=as.integer(P), 
           D=as.integer(d), N=as.integer(N), NE=as.integer(NE), NnonE=as.integer(NnonE), NM=as.integer(NM),
           G=as.integer(G), Y=as.numeric(t(Y)), E=as.integer(t(E)), nonE=as.integer(t(nonE)), M=as.integer(t(M)),
  	   numedges=as.integer(t(numedges)), EnonE=as.integer(t(EnonE)), diam=as.integer(diam),
	   hopslist=as.integer(t(hopslist)), XX=as.double(t(XX)),
           V_xi=as.double(V_xi), V_psi2=as.double(V_psi2), V_z=as.double(t(V_z)),
           V_sigma2=as.double(V_sigma2), V_eta=as.double(t(V_eta)),
           V_lambda=as.double(t(V_lambda)),
           V_omega2=as.double(V_omega2), V_nu=as.double(V_nu), V_alpha=as.double(V_alpha),
           xi=as.double(xi), psi2=as.double(psi2), sigma02=as.double(sigma02),
           omega2=as.double(omega2), nu=as.double(nu), alpha=as.double(alpha),
           inv_sigma02=as.double(inv_sigma02), tol=as.double(tol), STRAT=as.double(STRAT), 
	   seed=as.double(seed), d_vector=as.double(d_vector), conv=as.integer(conv),
	   PACKAGE="VBLPCM")
  
  V_xi<-out$V_xi
  V_z<-t(matrix(out$V_z,ncol=N))
  V_sigma2<-out$V_sigma2
  V_eta<-t(matrix(out$V_eta,ncol=G))
  V_omega2<-out$V_omega2
  V_lambda<-t(matrix(out$V_lambda,ncol=G))
  V_nu<-out$V_nu
  V_alpha<-out$V_alpha
  V_psi2<-out$V_psi2
  
  V_eta<-t(apply(V_eta, 1, "-", apply(V_z, 2, mean)))
  V_z<-t(apply(V_z, 1, "-", apply(V_z, 2, mean)))
  
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
  variational.params<-list()
  variational.params$net<-variational.start$net
  P->variational.params$P
  model->variational.params$model
  d->variational.params$d
  N->variational.params$N
  NE->variational.params$NE
  NnonE->variational.params$NnonE
  NM->variational.params$NM
  G->variational.params$G
  Y->variational.params$Y
  E->variational.params$E
  nonE->variational.params$nonE
  M->variational.params$M
  numedges->variational.params$numedges
  EnonE->variational.params$EnonE
  diam->variational.params$diam
  hopslist->variational.params$hopslist
  XX->variational.params$XX
  V_xi->variational.params$V_xi
  V_psi2->variational.params$V_psi2
  V_z->variational.params$V_z
  V_sigma2->variational.params$V_sigma2
  V_eta->variational.params$V_eta
  V_lambda->variational.params$V_lambda
  V_omega2->variational.params$V_omega2
  V_nu->variational.params$V_nu
  V_alpha->variational.params$V_alpha
  xi->variational.params$xi
  psi2->variational.params$psi2
  sigma02->variational.params$sigma02
  omega2->variational.params$omega2
  nu->variational.params$nu
  alpha->variational.params$alpha
  inv_sigma02->variational.params$inv_sigma02
  STRAT->variational.params$STRAT
  as.logical(out$conv)->variational.params$conv
  seed->variational.params$seed # this is the original seed value, not the value the RNG is now using
  BIC<-vblpcmbic(variational.params)
  BIC->variational.params$BIC
  class(variational.params)<-"vblpcm"
  return(variational.params)
  }
