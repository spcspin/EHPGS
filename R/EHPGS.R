#' Evaluation of hybrid performance in plant breeding via genomic selection
#'
#' @param num.F1        An integer number of superior hybrid combinations are select.
#' @param num.P         An integer number of potential parental lines are select.
#' @param train.pheno   A numeric vector of training set's phenotype.
#' @param train.geno    A numerical matrix of training set's genotype. Matrix is coded as -1,0 and 1 for Minor, Hetero, and Major, respectively.
#' @param parent.geno   A numerical matrix of parental lines' genotype. Matrix is coded as -1,0 and 1 for Minor, Hetero, and Major, respectively.
#' @param train.in.test A character vector of lines' names which are both in training set and half-diallel design.
#' @param mu.ini        Initial values for the constant term.
#' @param ga.ini        Initial values for the additive genotypic value.
#' @param gd.ini        Initial values for the dominance genotypic value
#' @param vE.ini        Initial values for the error variance.
#' @param vA.ini        Initial values for the additive variance.
#' @param vD.ini        Initial values for the dominance variance.
#' @param iter          Number of iteration in each chain.
#' @param m             Number of independent chain.
#'
#' @return This function will return 3 data frames, F1.out includes the top n GEBV F1s (n=num.F1) and its GEBV,SCA,MPH, and BPH. P.out includes the top m GEBV parents (m=num.P) and its GEBV,GCA. V.out includes the estimate of constant term and variance components (additive variance,dominance variance,and error variance)
#' @export
#' @import MASS stats utils
#' @examples
#' EHPGS(num.F1=15,num.P=5,train.pheno=BLUP.pheno,train.geno,parent.geno,
#' train.in.test=rownames(train.geno),mu.ini=NULL,
#' ga.ini = 0,gd.ini = 0,vE.ini = 1,vA.ini = 0.5,vD.ini = 0.5,iter = 100,m = 1)
EHPGS <- function(num.F1,num.P,train.pheno,train.geno,parent.geno,train.in.test=NULL,mu.ini=NULL,ga.ini = 0,gd.ini = 0,vE.ini = 1,vA.ini = 0.5,vD.ini = 0.5,iter = 10000,m = 5){
  ### set prior parameters

  S_star = 0.5*stats::var(train.pheno,na.rm =TRUE)
  nu_star = 5

  if (is.null(mu.ini)) {
    mu.ini = mean(train.pheno,na.rm =TRUE)
  }

  ###generate half-diallel design hybrid genotypes

  hybrid.geno <- c()
  F1 <- c()
  P1 <- c()
  P2 <- c()

  for (i in seq(1,nrow(parent.geno)-1)) {
    for (j in seq(i+1,nrow(parent.geno))) {
      P1 <- c(P1,rownames(parent.geno)[i])
      P2 <- c(P2,rownames(parent.geno)[j])
      F1 <- c(F1,paste(rownames(parent.geno)[i],":",rownames(parent.geno)[j],sep = ""))
      hybrid <- (parent.geno[i,]+parent.geno[j,])/2
      hybrid.geno <- rbind(hybrid.geno,hybrid)

    }

  }

  if (is.null(train.in.test)==FALSE) {
    dup.dat <- data.frame(do.call("rbind", strsplit(train.in.test, ":", fixed = TRUE)))
    colnames(dup.dat) <- c("P1","P2")
    re.train.in.test <- paste(dup.dat$P2,":",dup.dat$P1,sep = "")
    for (i in 1:length(train.in.test)) {
      if (re.train.in.test[i]%in%F1) {
        F1[match(re.train.in.test[i],F1)] <- train.in.test[i]
      }
    }
  }

  rownames(hybrid.geno) <- F1

 
  index = match(train.in.test,F1)
  de.hybrid.geno <- hybrid.geno[-index,]
  Xa <- rbind(train.geno,de.hybrid.geno)
  Xd <- replace(Xa,Xa==0,99)
  Xd <- replace(Xd,Xd==1,0)
  Xd <- replace(Xd,Xd==-1,0)
  Xd <- replace(Xd,Xd==99,1)
  KA <- kinship(Xa)
  KD <- kinship(Xd)

  rownames(KA) = rownames(KD) = colnames(KA) = colnames(KD) = c(rownames(train.geno),rownames(de.hybrid.geno))
  
  Ka <- KA[rownames(train.geno),rownames(train.geno)]
  Kd <- KD[rownames(train.geno),rownames(train.geno)]


  matrix.inverse=FALSE
  k=-9
  Ka.ori <- Ka
  Kd.ori <- Kd
  while (matrix.inverse==FALSE) {
    Ka.check <- class(tryCatch(solve(Ka), error = function(e) e))
    Kd.check <- class(tryCatch(solve(Kd), error = function(e) e))
    if (any(Ka.check=="error")|any(Kd.check=="error")) {
      if (k<=-4) {

        Ka <- Ka.ori +diag(10^(k),nrow = nrow(Ka), ncol = ncol(Ka))
        Kd <- Kd.ori +diag(10^(k),nrow = nrow(Kd), ncol = ncol(Kd))
        k=k+1

      }else{matrix.inverse=TRUE}
    }else{matrix.inverse=TRUE}

  }
  invisible(utils::capture.output(test <- BGS(pheno = train.pheno,KA = Ka,KD = Kd,mu.ini ,ga.ini,gd.ini ,vE.ini,vA.ini,vD.ini,iter = 1000,m = 1,S_star,nu_star)))
  vA <- mean(apply(test[[4]],1,function(x) mean(x[91:100])))
  vD <- mean(apply(test[[5]],1,function(x) mean(x[91:100])))
  vAvDratio=FALSE
  while (vAvDratio==FALSE) {
    if (vA/vD>10|vA/vD<0.1) {
      if (k<=-4) {
        Ka <- Ka.ori +diag(10^(k),nrow = nrow(Ka), ncol = ncol(Ka))
        Kd <- Kd.ori +diag(10^(k),nrow = nrow(Kd), ncol = ncol(Kd))
        invisible(capture.output(test <- BGS(pheno = train.pheno,KA = Ka,KD = Kd,mu.ini ,ga.ini,gd.ini ,vE.ini,vA.ini,vD.ini,iter = 100,m = 1,S_star,nu_star)))
        vA <- mean(apply(test[[4]],1,function(x) mean(x[91:100])))
        vD <- mean(apply(test[[5]],1,function(x) mean(x[91:100])))
        k=k+1
      }else{vAvDratio=TRUE}
    }else{vAvDratio=TRUE}

  }
  result <- BGS(pheno = train.pheno,KA = Ka,KD = Kd,mu.ini ,ga.ini,gd.ini ,vE.ini,vA.ini,vD.ini,iter ,m ,S_star,nu_star)
  vA <- mean(apply(result[[4]],1,function(x) mean(x[seq(iter-iter/10+1,iter)])))
  vD <- mean(apply(result[[5]],1,function(x) mean(x[seq(iter-iter/10+1,iter)])))
  vE <- mean(apply(result[[6]],1,function(x) mean(x[seq(iter-iter/10+1,iter)])))
  mu <-mean(apply(result[[3]],1,function(x) mean(x[seq(iter-iter/10+1,iter)])))

  ga <- apply(result[[1]],3,function(x) mean(x[seq(iter-iter/10+1,iter)]))
  gd <- apply(result[[2]],3,function(x) mean(x[seq(iter-iter/10+1,iter)]))

  if (length(train.in.test)<length(F1)) {
    ga <- KA[,match(rownames(Ka),rownames(KA))] %*% solve(Ka) %*% ga
    gd <- KD[,match(rownames(Kd),rownames(KD))] %*% solve(Kd) %*% gd
  }else{
    F1 <- rownames(train.geno)
    P1 <- data.frame(do.call("rbind", strsplit(F1, ":", fixed = TRUE)))[,1]
    P2 <- data.frame(do.call("rbind", strsplit(F1, ":", fixed = TRUE)))[,2]
  }

  BGS.allga <- as.data.frame(cbind(ga=ga,P1=P1,P2=P2,F1=F1))
  BGS.allgd <- as.data.frame(cbind(gd=gd,P1=P1,P2=P2,F1=F1))

  colnames(BGS.allga) <-c("ga","P1" ,"P2" ,"F1")
  colnames(BGS.allgd) <-c("gd","P1" ,"P2" ,"F1")
  BGS.allga$ga <-as.numeric(BGS.allga$ga)
  BGS.allgd$gd <-as.numeric(BGS.allgd$gd)

  #calculate GCA

  P <- union(BGS.allga$P1,BGS.allga$P2)
  GCA <- rep(NA,length(P))
  Ga.bar <- mean(BGS.allga$ga)
  for (i in 1:length(P)) {
    N0 <- length(which(BGS.allga$P1==P[i]))+length(which(BGS.allga$P2==P[i]))
    ga.i <- BGS.allga$ga[which(BGS.allga$P1==P[i]|BGS.allga$P2==P[i])]
    Ga.bar.i <- sum(ga.i)/(N0-1)
    GCA[i] <- ((N0-1)*Ga.bar.i)/(N0-2)-N0*Ga.bar/(2*(N0-2))
  }
  names(GCA) <- P

  #calculate SCA
  F1 <- BGS.allgd$F1
  SCA <- BGS.allgd$gd
  names(SCA) <- F1

  #calculate BPH
  BPH <- rep(NA,length(F1))
  for (i in 1:length(F1)) {
    GCA.i <- GCA[which(names(GCA)==P1[i])]
    GCA.j <- GCA[which(names(GCA)==P2[i])]
    SCA.ij <- SCA[which(names(SCA)==F1[i])]
    BPH[i] <- SCA.ij-abs(GCA.i-GCA.j)

  }
  names(BPH) <- F1

  #calculate GEBV.F1
  GEBV.F1 <-rep(NA,length(SCA))
  for (i in 1:length(F1)) {
    GEBV.F1[i] <- mu+ BGS.allga$ga[i]+BGS.allgd$gd[i]
  }
  names(GEBV.F1) <- F1

  #calculate GEBV.P
  P <- names(GCA)
  GEBV.P <-mu+2*GCA
  names(GEBV.P) <- P

  F1.result <- as.data.frame(cbind(GEBV=round(GEBV.F1,3),SCA=round(SCA,3),MPH=round(SCA,3),BPH=round(BPH,3)))
  P.result <- as.data.frame(cbind(GEBV=round(GEBV.P,3),GCA=round(GCA,3)))

  F1.out <- F1.result[order(F1.result$GEBV, decreasing = TRUE)[1:num.F1], ]
  P.out <- P.result[order(P.result$GEBV, decreasing = TRUE)[1:num.P],]
  V.out <- as.data.frame(cbind(mu=round(mu,3),vA=round(vA,3),vD=round(vD,3),vE=round(vE,3)))

  output <- list(F1.out=F1.out,P.out=P.out,V.out=V.out )
  return(output)

}
