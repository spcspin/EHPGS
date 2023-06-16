#' The Bayesian Gibbs sampling algorithm
#'
#' @param pheno   A numeric vector of phenotype which include all the lines in relationship matrix.
#' @param KA      Additive relationship matrix.
#' @param KD      Dominance relationship matrix.
#' @param mu.ini  Initial values for the constant term.
#' @param ga.ini  Initial values for the additive genotypic value.
#' @param gd.ini  Initial values for the dominance genotypic value
#' @param vE.ini  Initial values for the error variance.
#' @param vA.ini  Initial values for the additive variance.
#' @param vD.ini  Initial values for the dominance variance.
#' @param iter    Number of iteration in each chain.
#' @param m       Number of independent chain.
#' @param S_star  The shape of freedom of the variance components for the prior distribution of a scaled inverse Chi square.
#' @param nu_star The degree of freedom of the variance components for the prior distribution of a scaled inverse Chi square.
#' @param cat.itr A boolean variable determines whether to print out the iteration number.
#'
#' @return This function will return the MCMC value of additive genotypic value (g_A;ga), dominance genotypic value (g_D;gd), constant term (\eqn{mu};mu), and variance components(additive variance;vA, dominance variance;vD, error variance;vE)
#' @export
#' @import MASS stats utils
#' @examples
#'
#' KA <- kinship(train.geno)
#' Xd <- replace(train.geno,train.geno==0,99)
#' Xd <- replace(Xd,Xd==1,0)
#' Xd <- replace(Xd,Xd==-1,0)
#' Xd <- replace(Xd,Xd==99,1)
#'
#' KD <- kinship(Xd)
#' KA <- KA +diag(10^(-9),nrow = nrow(KA), ncol = ncol(KA))
#' KD <- KD +diag(10^(-9),nrow = nrow(KD), ncol = ncol(KD))
#'
#' result <- BGS(pheno =train.pheno$F1.weight,KA = KA,KD = KD,mu.ini = 0,ga.ini = 0,gd.ini = 0,
#' vE.ini = 1,vA.ini = 0.5,vD.ini = 0.5,iter = 100,m = 1,S_star = 0.5*var(train.pheno$F1.weight),nu_star = 5)
BGS <- function(pheno,KA,KD,mu.ini=0,ga.ini=0,gd.ini=0,vE.ini=1,vA.ini=0,vD.ini=0,iter=5000,m=1,S_star,nu_star=5,cat.itr=TRUE){
  y=pheno
  n=length(pheno)
  n_fix=1

  if (is.null(mu.ini)) {
    mu.ini = mean(train.pheno,na.rm =TRUE)
  }

  ga_dat <- array(NA,c(m,iter,n))
  gd_dat <- array(NA,c(m,iter,n))
  mu_dat <- matrix(NA,m,iter)
  vA_dat <- matrix(NA,m,iter)
  vD_dat <- matrix(NA,m,iter)
  vE_dat <- matrix(NA,m,iter)


  for (i in 1:m) {
    if (cat.itr) {
      cat(paste("\n","#chain=",i,"\n",sep = ""))
    }
    for (j in 0:iter) {
      if (j==0) {
        if (cat.itr) {
          cat("iter=")
        }
        vA=vA.ini;vD=vD.ini;vE=vE.ini
        lamA=vE/vA;lamD=vE/vD
        X=rep(1,n)
        Z=diag(1,n,n)
        W= cbind(X,Z,Z)


        bl11=matrix(0,n_fix,n_fix);
        bl12=matrix(0,n_fix,n);
        bl21=matrix(0,n,n_fix);
        bl13=matrix(0,n_fix,n);
        bl31=matrix(0,n,n_fix);
        bl23=matrix(0,n,n);
        bl32=matrix(0,n,n);

        E=rbind(cbind(bl11,bl12,bl13),cbind(bl21,solve(KA)*lamA,bl23),cbind(bl31,bl32,solve(KD)*lamD))

        C= t(W)%*%W+E

        r1=t(rep(1,n))%*%y
        r2=y
        r3=y

        g1=as.matrix(rep(mu.ini,n_fix))
        g2=rep(ga.ini,n)
        g3=rep(gd.ini,n)

      }

      if (cat.itr) {
        cat(".")
      }

      g1_star=solve(C[1:n_fix,1:n_fix])%*%(r1-C[n_fix,seq(n_fix+1,n_fix+n)]%*%g2-C[n_fix,seq(n_fix+n+1,n_fix+n+n)]%*%g3)
      g1=MASS::mvrnorm(n=1,mu=g1_star,Sigma = vE*solve(C[n_fix,n_fix]))
      g2_star=solve(C[seq(n_fix+1,n_fix+n),seq(n_fix+1,n_fix+n)])%*%(r2-C[seq(n_fix+1,n_fix+n),n_fix]%*%as.matrix(g1)-C[seq(n_fix+1,n_fix+n),seq(n_fix+n+1,n_fix+n+n)]%*%g3)
      g2=MASS::mvrnorm(n=1,mu=g2_star,Sigma = vE*solve(C[seq(n_fix+1,n_fix+n),seq(n_fix+1,n_fix+n)]))
      g3_star=solve(C[seq(n_fix+n+1,n_fix+n+n),seq(n_fix+n+1,n_fix+n+n)])%*%(r3-C[seq(n_fix+n+1,n_fix+n+n),n_fix]%*%as.matrix(g1)-C[seq(n_fix+n+1,n_fix+n+n),seq(n_fix+1,n_fix+n)]%*%g2)
      g3=MASS::mvrnorm(n=1,mu=g3_star,Sigma = vE*solve(C[seq(n_fix+n+1,n_fix+n+n),seq(n_fix+n+1,n_fix+n+n)]))

      e=y-W%*%c(g1,g2,g3)

      S_star=0.5*stats::var(y);nu_star=5

      vE=(t(e)%*%e+S_star*nu_star)/stats::rchisq(1,n+nu_star)

      vA=(t(g2)%*%solve(KA)%*%g2+S_star*nu_star)/stats::rchisq(1,n+nu_star)

      vD=(t(g3)%*%solve(KD)%*%g3+S_star*nu_star)/stats::rchisq(1,n+nu_star)

      vE <-as.numeric(vE)

      lamA=vE/vA;lamD=vE/vD

      E=rbind(cbind(bl11,bl12,bl13),cbind(bl21,solve(KA)*as.numeric(lamA),bl23),cbind(bl31,bl32,solve(KD)*as.numeric(lamD)))

      C= t(W)%*%W+E

      ga_dat[i,j,] <- g2
      gd_dat[i,j,] <- g3
      mu_dat[i,j] <- g1
      vA_dat[i,j] <- vA
      vD_dat[i,j] <- vD
      vE_dat[i,j] <- vE

      if (cat.itr) {
        cat(j,sep="")
      }
    }

  }

  result=list(ga=ga_dat,gd=gd_dat,mu=mu_dat,vA=vA_dat,vD=vD_dat,vE=vE_dat)
  return(result)

}
