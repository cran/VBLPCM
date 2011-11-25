# lower BIC is preferred
vblpcmbic<-function(v.params)
  {
  N<-v.params$N
  NE<-v.params$NE
  Y<-v.params$Y
  P<-v.params$P
  XX<-v.params$XX
  V_xi<-v.params$V_xi
  V_psi2<-v.params$V_psi2
  V_z<-v.params$V_z
  V_sigma2<-v.params$V_sigma2
  V_omega2<-v.params$V_omega2
  V_eta<-v.params$V_eta
  V_lambda<-v.params$V_lambda
  G<-v.params$G
  d<-v.params$d
  cov1<-XX%*%V_xi
  cov2<-XX%*%V_psi2
  LL<-0
  for (i in 1:N)
    for (j in (1:N)[-i])
    if (!is.na(Y[i,j]))
      {
      tmp<-cov1[(i-1)*N+j]-sqrt(t(V_z[i,]-V_z[j,])%*%(V_z[i,]-V_z[j,]) + d*(V_sigma2[i]+V_sigma2[j]))
      LL = LL + Y[i,j]*tmp - log(1+exp(tmp + 0.5*cov2[(i-1)*N+j]))
      }
  BIC<-list(
           Y=-2*LL + (N*d+(P-1)*N)*log(N),
           #Y=-2*LL + (P+1)*log(NE),
	   MBC =
                                if(G>0){
				  tmp=0
				  for (g in 1:G)
				  for (i in 1:N)
				    tmp =
				    tmp-2*sum(V_lambda[g,i]*dnorm(V_z[i,],V_eta[g,],sqrt(V_sigma2[i]+V_omega2[g]),log=TRUE))
				    tmp=tmp+(G+G*d+G-1)*log(N)
				    tmp
                                } else {
                                  -2*sum(dnorm(V_z,0,sqrt(mean(V_z^2)*d),log=TRUE))+1*log(N*d)
                                }
           )
  BIC[["overall"]]<-sum(unlist(BIC))
  BIC
  }
