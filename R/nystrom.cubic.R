#' @title Nystrom Approximation for Cubic Smoothing Spline Estimates
#' @description Computes a Nystrom approximation based on randomly selected columns of $\Sigma$ matrix, where estimation utilizes functions for linear mixed effect model (LME).
#' @param x The values of independent variable. It should be a vector.
#' @param y The values of dependent variable. It should be a vector.
#' @param xg The gridpoints used for approximation of
#' @param e A list of two elements of "values" and "vectors", which refer respectively, the eigenvalues and eigenfunctions of reproducing kernel for cubic splines at the pre-selected gridpoints (must agree with xg).
#' @param p An integer value. The selection parameter indicates the number of columns for random selection and approximation. The default value is 30.
#' @param method A character string. If "REML" the LME model is fit by maximizing the restricted log-likelihood. If "ML" the log-likelihood is maximized. Defaults to "REML".
#' @param pstd An indicator of whether standard deviation is desired. The default value is FALSE.
#' @export
#' @import mgcv
#' @return A vector(s) of following component(s):
#' \item{fit}{The Nystrom approximation of cubic smoothing spline estimate.}
#' \item{pstd}{The corrsponding posterior standard deviation.}
#' @examples
#' \dontrun{
#' data(eigenM)
#' x <- runif(1000)
#' y <- sin(32*pi*x)-8*(x-.5)^2 + rnorm(1000)
#' nystrom.cubic(x,y,xg,e,K,method="REML",pstd=FALSE)
#' }

nystrom.cubic <- function(x, y, p=30, method="REML", pstd=FALSE){
  n <- length(y)
  selcol <- sort(sample(1:n,p))
  W <- cubic(x[selcol])
  w <- eigen(W)
  C <-cubic(x,x[selcol])
  Z <- t(1/sqrt(w$values)*t(C%*%w$vectors))
  all1 <- rep(1,n)
  tmp <- lme(y~x, random=list(all1=pdIdent(~Z-1)),method = method)
  if (pstd){
    T <-matrix(c(rep(1,n),x),nr=n)
    #sigmaM <- cubic(x)
    v <- as.numeric(VarCorr(tmp)[p:(p+1),1])
    A <- solve(crossprod(Z)+v[2]/v[1]*diag(p))
    M.inv <-v[1]/v[2]*(diag(n)-Z%*%( tcrossprod( A ,Z) ) )
    TMT <- solve(t(T)%*%M.inv%*%T)
    #H <- M.inv-(M.inv%*%T)%*%TMT%*%(crossprod(T,M.inv))
    var1 <- v[1]*T%*%(tcrossprod(TMT,T))
    #var2 <- v[1]*(tcrossprod(Z)-Z%*%tcrossprod(crossprod(Z,H)%*%Z,Z))
    var2 <- v[2]*Z%*%(tcrossprod(A,Z))+v[1]*(Z%*%A)%*%crossprod(Z,T)%*%TMT%*%crossprod(T,Z)%*%tcrossprod(A,Z)
    cov <- -v[1]*(T%*%TMT)%*%(crossprod(T,M.inv)%*%tcrossprod(Z))
    VarCov <- var1+var2+2*cov
    pred <- list(fit=as.vector(tmp$fit[,2]),pstd=as.vector(sqrt(diag(VarCov))))
  } else{
    pred <- as.vector(tmp$fit[,2])
  }
  return(pred)
}

