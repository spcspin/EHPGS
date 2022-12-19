#' Relationship matrix
#'
#'This function calculate additive and dominance relationship matrix.
#'
#' @param X  A numerical matrix of genotype. If is additive relationship matrix can be coded as -1,0,and 1.If is dominance relationship matrix can be coded as 0,1,and 0 (row: sample; column: marker).
#'
#' @return K A numerical matrix (relationship matrix)
#' @export
#' @import MASS stats utils
#' @examples
#' KA <- kinship(train.geno)
kinship <- function(X){
  X <-as.matrix(X[,which(apply(X,2,stats::sd)!=0)])
  sta <- apply(X,2,scale)
  sta <- t(stats::na.omit(t(sta)))
  K <- sta%*%t(sta)/ncol(sta)
  return(K)
}
